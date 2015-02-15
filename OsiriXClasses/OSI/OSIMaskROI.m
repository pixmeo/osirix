//
//  OSIMaskROI.m
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 9/26/12.
//  Copyright (c) 2012 OsiriX Team. All rights reserved.
//

#import "OSIMaskROI.h"
#import "CPRGenerator.h"
#import "CPRGeneratorRequest.h"
#import "OSIFloatVolumeData.h"
#import "Notifications.h"
#include <OpenGL/CGLMacro.h>
#include <Accelerate/Accelerate.h>

@interface OSIMaskROI ()
@property (nonatomic, readwrite, retain) NSString *name;
@property (nonatomic, readwrite, retain) NSColor *fillColor;
- (NSData *)_maskRunsDataForSlab:(OSISlab)slab dicomToPixTransform:(N3AffineTransform)dicomToPixTransform minCorner:(N3VectorPointer)minCornerPtr;
- (N3Vector)_cachedBitmapMaskCorner;
- (BOOL)_isMaskInCachedBitmapMask:(OSIROIMask *)mask;
@end

@implementation OSIMaskROI

@synthesize mask = _mask;
@synthesize volumeTransform = _volumeTransform;
@synthesize name = _name;
@synthesize fillColor = _fillColor;

- (instancetype)initWithROI:(OSIROI *)roi sampledOnVolumeData:(OSIFloatVolumeData *)floatVolumeData // initializing by copying the values and getting the mask from the passed in ROI.
{
    if ( (self = [super init])) {
        _mask = [[roi ROIMaskForFloatVolumeData:floatVolumeData] retain];
        _volumeTransform = [floatVolumeData volumeTransform];
        _name = [[roi name] copy];
    }
    return self;
}

- (instancetype)initWithROIMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform name:(NSString *)name
{
	if ( (self = [super init]) ) {
        _mask = [mask retain];
        _volumeTransform = volumeTransform;
        _name = [name copy];
	}
	return self;
}

- (instancetype)initWithROIMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform sampledOnVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name reinterpolateMask:(BOOL)reinterpolateMask
{
	if ( (self = [super init]) ) {
        if (reinterpolateMask = NO) {
            _mask = [mask retain];
            _volumeTransform = volumeTransform;
            _name = [name retain];
        } else {
            OSIROIMask *mappedMask = [self.mask ROIMaskByResamplingFromVolumeTransform:volumeTransform toVolumeTransform:floatVolumeData.volumeTransform interpolationMode:CPRInterpolationModeNearestNeighbor];
            _mask = [[mappedMask ROIMaskCroppedToWidth:floatVolumeData.pixelsWide height:floatVolumeData.pixelsHigh depth:floatVolumeData.pixelsDeep] retain];
            _volumeTransform = floatVolumeData.volumeTransform;
            _name = [name retain];
        }
	}
    return self;
}

- (instancetype)initWithCoder:(NSCoder *)decoder
{
    if ( (self = [super init]) ) {
        _mask = [[decoder decodeObjectOfClass:[OSIROIMask class] forKey:@"mask"] retain];
        _volumeTransform = [decoder decodeN3AffineTransformForKey:@"volumeTransform"];
        _name = [[decoder decodeObjectOfClass:[NSString class] forKey:@"name"] retain];
        _fillColor = [[decoder decodeObjectOfClass:[NSColor class] forKey:@"fillColor"] retain];
    }
    return self;
}

- (void)encodeWithCoder:(NSCoder *)encoder
{
    [encoder encodeObject:_mask forKey:@"mask"];
    [encoder encodeN3AffineTransform:_volumeTransform forKey:@"volumeTransform"];
    [encoder encodeObject:_name forKey:@"name"];
    [encoder encodeObject:_fillColor forKey:@"fillColor"];
}

- (OSIMaskROI *)maskROIRepresention
{
    return self;
}

- (void)dealloc
{
    self.mask = nil;
    self.name = nil;
    self.fillColor = nil;
    [_cachedBitmapMask release];
    _cachedBitmapMask = nil;

    [super dealloc];
}

- (void)setVolumeTransform:(N3AffineTransform)volumeTransform
{
    [self setVolumeTransform:volumeTransform reinterpolateMask:NO];
}

- (void)setVolumeTransform:(N3AffineTransform)volumeTransform reinterpolateMask:(BOOL)reinterpolateMask
{
    if (N3AffineTransformEqualToTransform(volumeTransform, _volumeTransform) == NO) {
        [_cachedBitmapMask release];
        _cachedBitmapMask = nil;

        if (reinterpolateMask = NO) {
            [self willChangeValueForKey:@"volumeTransform"];
            [self willChangeValueForKey:@"mask"]; // ok, so the mask itself didn't change, but what the mask logically represented did
            _volumeTransform = volumeTransform;
            [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self userInfo:nil];
            [self didChangeValueForKey:@"mask"];
            [self didChangeValueForKey:@"volumeTransform"];
        } else {
            [self willChangeValueForKey:@"volumeTransform"];
            [self willChangeValueForKey:@"mask"]; // ok, so the mask itself didn't change, but what the mask logically represented did
            OSIROIMask *mappedMask = [[self.mask ROIMaskByResamplingFromVolumeTransform:_volumeTransform toVolumeTransform:volumeTransform interpolationMode:CPRInterpolationModeNearestNeighbor] retain];
            [_mask release];
            _mask = mappedMask;
            _volumeTransform = volumeTransform;
            [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self userInfo:nil];
            [self didChangeValueForKey:@"mask"];
            [self didChangeValueForKey:@"volumeTransform"];
        }
    }
}

- (void)setMask:(OSIROIMask *)mask
{
    if (_mask != mask) {
        [_mask release];
        _mask = [mask retain];
        [_cachedBitmapMask release];
        _cachedBitmapMask = nil;
        [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self userInfo:nil];
    }
}

- (void)unionWithMask:(OSIROIMask *)mask
{
    [self unionWithMask:mask volumeTransform:_volumeTransform];
}

- (void)unionWithMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform
{
    if ([mask maskRunCount] == 0) {
        return;
    }

    mask = [mask ROIMaskByResamplingFromVolumeTransform:volumeTransform toVolumeTransform:_volumeTransform interpolationMode:CPRInterpolationModeNearestNeighbor];

    if ([self _isMaskInCachedBitmapMask:mask]) {
        NSAssert(_cachedBitmapMask != nil, @"this can't happen because _isMaskInCachedBitmapMask would have returned no, but assert to quiet the static analyzer");
        N3Vector corner = [self _cachedBitmapMaskCorner];

        OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[mask maskRunsData] bytes];
        NSInteger maskRunCount = [mask maskRunCount];
        NSInteger i;

        CPRVolumeDataInlineBuffer inlineBuffer;

        [_cachedBitmapMask aquireInlineBuffer:&inlineBuffer];
        float *bitmapMaskBytes = (float *)CPRVolumeDataFloatBytes(&inlineBuffer);
        NSInteger width = _cachedBitmapMask.pixelsWide;
        NSInteger height = _cachedBitmapMask.pixelsHigh;

        // draw in the runs
        for (i = 0; i < maskRunCount; i++) {
            NSInteger x = maskRuns[i].widthRange.location - corner.x;
            NSInteger y = maskRuns[i].heightIndex - corner.y;
            NSInteger z = maskRuns[i].depthIndex - corner.z;

            vDSP_vfill(&(maskRuns[i].intensity), &(bitmapMaskBytes[x + y*width + z*width*height]), 1, maskRuns[i].widthRange.length);
        }
        OSIROIMask *newMask = [self.mask ROIMaskByUnioningWithMask:mask];
        [self willChangeValueForKey:@"mask"];
        [_mask release];
        _mask = [newMask retain];
        [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self userInfo:nil];
        [self didChangeValueForKey:@"mask"];
    } else {
        self.mask = [self.mask ROIMaskByUnioningWithMask:mask];
    }
}

- (void)subtractMask:(OSIROIMask *)mask
{
    [self subtractMask:mask volumeTransform:_volumeTransform];
}

- (void)subtractMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform
{
    if ([mask maskRunCount] == 0) {
        return;
    }

    mask = [mask ROIMaskByResamplingFromVolumeTransform:volumeTransform toVolumeTransform:_volumeTransform interpolationMode:CPRInterpolationModeNearestNeighbor];

    if ([self _isMaskInCachedBitmapMask:mask]) {
        NSAssert(_cachedBitmapMask != nil, @"this can't happen because _isMaskInCachedBitmapMask would have returned no, but assert to quiet the static analyzer");

        N3Vector corner = [self _cachedBitmapMaskCorner];

        OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[mask maskRunsData] bytes];
        NSInteger maskRunCount = [mask maskRunCount];
        NSInteger i;

        CPRVolumeDataInlineBuffer inlineBuffer;

        [_cachedBitmapMask aquireInlineBuffer:&inlineBuffer];
        float *bitmapMaskBytes = (float *)CPRVolumeDataFloatBytes(&inlineBuffer);
        NSInteger width = _cachedBitmapMask.pixelsWide;
        NSInteger height = _cachedBitmapMask.pixelsHigh;

        // draw in the runs
        for (i = 0; i < maskRunCount; i++) {
            NSInteger x = maskRuns[i].widthRange.location - corner.x;
            NSInteger y = maskRuns[i].heightIndex - corner.y;
            NSInteger z = maskRuns[i].depthIndex - corner.z;

            memset(&(bitmapMaskBytes[x + y*width + z*width*height]), 0, maskRuns[i].widthRange.length * sizeof(float));
        }
        OSIROIMask *newMask = [self.mask ROIMaskBySubtractingMask:mask];
        [self willChangeValueForKey:@"mask"];
        [_mask release];
        _mask = [newMask retain];
        [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self userInfo:nil];
        [self didChangeValueForKey:@"mask"];
    } else {
        self.mask = [self.mask ROIMaskBySubtractingMask:mask];
    }
}

- (NSArray *)convexHull
{

    return [self.mask convexHull]; // TOTALLY BOGUS, this is in the wrong coordinate space!!!!!!!
}

- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume
{
    OSIROIMask *mappedMask = [self.mask ROIMaskByResamplingFromVolumeTransform:_volumeTransform toVolumeTransform:floatVolume.volumeTransform interpolationMode:CPRInterpolationModeNearestNeighbor];
    mappedMask = [mappedMask ROIMaskCroppedToWidth:floatVolume.pixelsWide height:floatVolume.pixelsHigh depth:floatVolume.pixelsDeep];

    return mappedMask;
}

- (OSIFloatVolumeData *)bitmapMask
{
    if (_cachedBitmapMask != nil) {
        return _cachedBitmapMask;
    }

    NSInteger maskRunCount = [self.mask maskRunCount];
    if (maskRunCount == 0) {
        return nil;
    }

    // figure out how bit the mask needs to be
    NSInteger maxWidth = NSIntegerMin;
    NSInteger minWidth = NSIntegerMax;
    NSInteger maxHeight = NSIntegerMin;
    NSInteger minHeight = NSIntegerMax;
    NSInteger maxDepth = NSIntegerMin;
    NSInteger minDepth = NSIntegerMax;

    OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[self.mask maskRunsData] bytes];
    NSInteger i;

    for (i = 0; i < maskRunCount; i++) {
        maxWidth = MAX(maxWidth, (NSInteger)OSIROIMaskRunLastWidthIndex(maskRuns[i]));
        minWidth = MIN(minWidth, (NSInteger)OSIROIMaskRunFirstWidthIndex(maskRuns[i]));

        maxHeight = MAX(maxHeight, (NSInteger)maskRuns[i].heightIndex);
        minHeight = MIN(minHeight, (NSInteger)maskRuns[i].heightIndex);

        maxDepth = MAX(maxDepth, (NSInteger)maskRuns[i].depthIndex);
        minDepth = MIN(minDepth, (NSInteger)maskRuns[i].depthIndex);
    }

    NSInteger width = (maxWidth - minWidth) + 1;
    NSInteger height = (maxHeight - minHeight) + 1;
    NSInteger depth = (maxDepth - minDepth) + 1;

    N3AffineTransform bitmapMaskVolumeTransform = N3AffineTransformConcat(_volumeTransform, N3AffineTransformMakeTranslation(-1.0*(CGFloat)minWidth, -1.0*(CGFloat)minHeight, -1.0*(CGFloat)minDepth));

    // create the FloatVolumeData
    float *bitmapMaskBytes = malloc(width * height * depth * sizeof(float));
    memset(bitmapMaskBytes, 0, width * height * depth * sizeof(float));
    OSIFloatVolumeData *bitmapMask = [[OSIFloatVolumeData alloc] initWithFloatBytesNoCopy:bitmapMaskBytes pixelsWide:width pixelsHigh:height pixelsDeep:depth volumeTransform:bitmapMaskVolumeTransform outOfBoundsValue:0 freeWhenDone:YES];

    // draw in the runs
    for (i = 0; i < maskRunCount; i++) {
        NSInteger x = maskRuns[i].widthRange.location - minWidth;
        NSInteger y = maskRuns[i].heightIndex - minHeight;
        NSInteger z = maskRuns[i].depthIndex - minDepth;

        vDSP_vfill(&(maskRuns[i].intensity), &(bitmapMaskBytes[x + y*width + z*width*height]), 1, maskRuns[i].widthRange.length);
    }

    _cachedBitmapMask = bitmapMask;
    return _cachedBitmapMask;
}


- (void)drawSlab:(OSISlab)slab inCGLContext:(CGLContextObj)cgl_ctx pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform
{
    OSIROIMaskRun maskRun;
    NSData *maskRunsData;
    N3Vector minCorner;
    NSInteger i;
    NSInteger runsCount;
    const OSIROIMaskRun *maskRunsBytes;
    double widthIndex;
    double maxWidthIndex;
    double heightIndex;
    double depthIndex;

    if (self.fillColor == nil) {
        return;
    }

    NSColor *deviceColor = [self.fillColor colorUsingColorSpaceName:NSDeviceRGBColorSpace];

    maskRunsData = [self _maskRunsDataForSlab:slab dicomToPixTransform:dicomToPixTransform minCorner:&minCorner];

    if ([maskRunsData length] == 0) {
        return;
    }

    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glColor4f((float)[deviceColor redComponent], (float)[deviceColor greenComponent], (float)[deviceColor blueComponent], (float)[deviceColor alphaComponent]);
    glBegin(GL_QUADS);
    runsCount = [maskRunsData length] / sizeof(OSIROIMaskRun);
    maskRunsBytes = [maskRunsData bytes];
    for (i = 0; i < runsCount; i++) {
        maskRun = maskRunsBytes[i];
        widthIndex = (double)maskRun.widthRange.location + minCorner.x;
        maxWidthIndex = widthIndex + (double)maskRun.widthRange.length;
        heightIndex = (double)maskRun.heightIndex + minCorner.y;
        depthIndex = maskRun.depthIndex;

        glVertex3d(widthIndex, heightIndex, depthIndex);
        glVertex3d(maxWidthIndex, heightIndex, depthIndex);
        glVertex3d(maxWidthIndex, heightIndex + 1.0, depthIndex);
        glVertex3d(widthIndex, heightIndex + 1.0, depthIndex);
    }
    glEnd();

    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_POLYGON_SMOOTH);
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);

}

- (CGFloat)volume
{
    N3AffineTransform maskToDicom = N3AffineTransformInvert(_volumeTransform);
    return [_mask maskIndexCount] * N3VectorLength(N3VectorApplyTransform(N3VectorXBasis, maskToDicom)) *
                                    N3VectorLength(N3VectorApplyTransform(N3VectorYBasis, maskToDicom)) *
                                    N3VectorLength(N3VectorApplyTransform(N3VectorZBasis, maskToDicom));
}

- (NSData *)_maskRunsDataForSlab:(OSISlab)slab dicomToPixTransform:(N3AffineTransform)dicomToPixTransform minCorner:(N3VectorPointer)minCornerPtr;
{
    CPRVolumeData *floatVolumeData;
    N3Vector corner;
    N3Vector minCorner;
    N3Vector maxCorner;
    N3AffineTransform coalescedVolumeMaskToPixTransform;
    NSInteger width;
    NSInteger height;
    CPRObliqueSliceGeneratorRequest *sliceRequest;
    OSIROIMask *sliceMask;

    minCorner = N3VectorMake(CGFLOAT_MAX, CGFLOAT_MAX, 0);
    maxCorner = N3VectorMake(-CGFLOAT_MAX, -CGFLOAT_MAX, 0);

    OSIFloatVolumeData *bitmapMask = [self bitmapMask];
    if (bitmapMask == nil) {
        return nil;
    }

    coalescedVolumeMaskToPixTransform = N3AffineTransformConcat(N3AffineTransformInvert(bitmapMask.volumeTransform), dicomToPixTransform);

    // first off figure out where this float volume needs to be
    corner = N3VectorApplyTransform(N3VectorMake(0,                     0,                     0),                     coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(bitmapMask.pixelsWide, 0,                     0),                     coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                     bitmapMask.pixelsHigh, 0),                     coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(bitmapMask.pixelsWide, bitmapMask.pixelsHigh, 0),                     coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                     0,                     bitmapMask.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(bitmapMask.pixelsWide, 0,                     bitmapMask.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                     bitmapMask.pixelsDeep, bitmapMask.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(bitmapMask.pixelsWide, bitmapMask.pixelsDeep, bitmapMask.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);

    minCorner.x = floor(minCorner.x) - 1;
    minCorner.y = floor(minCorner.y) - 1;
    maxCorner.x = ceil(maxCorner.x) + 1;
    maxCorner.y = ceil(maxCorner.y) + 1;


    width = maxCorner.x - minCorner.x;
    height = maxCorner.y - minCorner.y;

    sliceRequest = [[[CPRObliqueSliceGeneratorRequest alloc] init] autorelease];
    sliceRequest.pixelsWide = width;
    sliceRequest.pixelsHigh = height;
    sliceRequest.slabWidth = slab.thickness;
    sliceRequest.projectionMode = CPRProjectionModeMIP;
    sliceRequest.interpolationMode = CPRInterpolationModeNearestNeighbor;

    sliceRequest.sliceToDicomTransform = N3AffineTransformConcat(N3AffineTransformMakeTranslation(minCorner.x + 0.5, minCorner.y + 0.5, 0), N3AffineTransformInvert(dicomToPixTransform));

    floatVolumeData = [CPRGenerator synchronousRequestVolume:sliceRequest volumeData:bitmapMask];
    sliceMask = [OSIROIMask ROIMaskFromVolumeData:(OSIFloatVolumeData *)floatVolumeData volumeTransform:NULL];

    if (minCornerPtr) {
        *minCornerPtr = minCorner;
    }

    return [sliceMask maskRunsData];
}

- (N3Vector)_cachedBitmapMaskCorner
{
    if (_cachedBitmapMask == nil) {
        return N3VectorZero;
    }

    N3Vector corner = N3VectorApplyTransform(N3VectorZero, N3AffineTransformConcat(N3AffineTransformInvert(_cachedBitmapMask.volumeTransform), _volumeTransform));
    return N3VectorRound(corner);
}

- (BOOL)_isMaskInCachedBitmapMask:(OSIROIMask *)mask
{
    if (_cachedBitmapMask == nil) {
        return NO;
    }

    if ([mask maskRunCount] == 0) {
        return YES;
    }

    // find the bounds of the cached bitmap
    // first find the corner of the cached bitmap mask
    N3Vector corner = [self _cachedBitmapMaskCorner];
    NSInteger cachedMinWidth = (NSInteger)corner.x;
    NSInteger cachedMinHeight = (NSInteger)corner.y;
    NSInteger cachedMinDepth = (NSInteger)corner.z;
    NSInteger cachedMaxWidth = (cachedMinWidth + _cachedBitmapMask.pixelsWide) - 1;
    NSInteger cachedMaxHeight = (cachedMinHeight + _cachedBitmapMask.pixelsHigh) - 1;
    NSInteger cachedMaxDepth = (cachedMinDepth + _cachedBitmapMask.pixelsDeep) - 1;

    NSInteger maxWidth = NSIntegerMin;
    NSInteger minWidth = NSIntegerMax;
    NSInteger maxHeight = NSIntegerMin;
    NSInteger minHeight = NSIntegerMax;
    NSInteger maxDepth = NSIntegerMin;
    NSInteger minDepth = NSIntegerMax;

    OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[mask maskRunsData] bytes];
    NSInteger maskRunCount = [mask maskRunCount];
    NSInteger i;

    for (i = 0; i < maskRunCount; i++) {
        maxWidth = MAX(maxWidth, (NSInteger)OSIROIMaskRunLastWidthIndex(maskRuns[i]));
        minWidth = MIN(minWidth, (NSInteger)OSIROIMaskRunFirstWidthIndex(maskRuns[i]));

        maxHeight = MAX(maxHeight, (NSInteger)maskRuns[i].heightIndex);
        minHeight = MIN(minHeight, (NSInteger)maskRuns[i].heightIndex);

        maxDepth = MAX(maxDepth, (NSInteger)maskRuns[i].depthIndex);
        minDepth = MIN(minDepth, (NSInteger)maskRuns[i].depthIndex);
    }

    if (cachedMinWidth <= minWidth && cachedMinHeight <= minHeight && cachedMinDepth <= minDepth &&
        cachedMaxWidth >= maxWidth && cachedMaxHeight >= maxHeight && cachedMaxDepth >= maxDepth) {
        return YES;
    } else {
        return NO;
    }
}

- (OSIROIMask *)ROIMaskWithSphereDiameter:(CGFloat)diameter atDicomVector:(N3Vector)position
{
    NSAssert(diameter > 0,  @"%s called with negative diameter", __PRETTY_FUNCTION__);

    N3AffineTransform inverseTransform = N3AffineTransformInvert([self volumeTransform]);

    // Make round to the nearest odd number
#if CGFLOAT_IS_DOUBLE
    NSUInteger width = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorXBasis, inverseTransform)))/2.0));
    NSUInteger height = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorYBasis, inverseTransform)))/2.0));
    NSUInteger depth = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorZBasis, inverseTransform)))/2.0));
#else
    NSUInteger width = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorXBasis, inverseTransform)))/2.0));
    NSUInteger height = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorYBasis, inverseTransform)))/2.0));
    NSUInteger depth = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorZBasis, inverseTransform)))/2.0));
#endif

    OSIROIMask *sphereMask = [OSIROIMask ROIMaskWithElipsoidWidth:width height:height depth:depth];

    N3Vector maskPosition = N3VectorApplyTransform(position, [self volumeTransform]);
    maskPosition = N3VectorRound(maskPosition);
    maskPosition = N3VectorSubtract(maskPosition, N3VectorMake((width-1)/2, (height-1)/2, (depth-1)/2));

    return [sphereMask ROIMaskByTranslatingByX:maskPosition.x Y:maskPosition.y Z:maskPosition.z];
}


- (OSIROIMask *)ROIMaskBrushWithSphereDiameter:(CGFloat)diameter fromDicomVector:(N3Vector)fromDicomVector toDicomVector:(N3Vector)toDicomVector
{
    NSAssert(diameter > 0,  @"%s called with negative diameter", __PRETTY_FUNCTION__);

    N3AffineTransform inverseTransform = N3AffineTransformInvert([self volumeTransform]);

    // Make round to the nearest odd number
#if CGFLOAT_IS_DOUBLE
    NSUInteger width = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorXBasis, inverseTransform)))/2.0));
    NSUInteger height = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorYBasis, inverseTransform)))/2.0));
    NSUInteger depth = 1+(2*floor((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorZBasis, inverseTransform)))/2.0));
#else
    NSUInteger width = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorXBasis, inverseTransform)))/2.0));
    NSUInteger height = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorYBasis, inverseTransform)))/2.0));
    NSUInteger depth = 1+(2*floorf((diameter / N3VectorLength(N3VectorApplyTransformToDirectionalVector(N3VectorZBasis, inverseTransform)))/2.0));
#endif

    OSIROIMask *sphereMask = [OSIROIMask ROIMaskWithElipsoidWidth:width height:height depth:depth];

    N3Vector fromMaskVector = N3VectorApplyTransform(fromDicomVector, [self volumeTransform]);
    N3Vector toMaskVector = N3VectorApplyTransform(toDicomVector, [self volumeTransform]);

    NSInteger stampCount = N3VectorDistance(fromMaskVector, toMaskVector) + 1;
    OSIROIMask *mask = [OSIROIMask ROIMask];

    NSInteger i;
    const CGFloat sphereSize = 2;
    for (i = 0; i < stampCount; i++) {
        N3Vector brushPosition = N3VectorRound(N3VectorLerp(fromMaskVector, toMaskVector, (CGFloat)i/(CGFloat)stampCount));
        OSIROIMask *brushMask = [sphereMask ROIMaskByTranslatingByX:brushPosition.x - (width-1)/2 Y:brushPosition.y - (height-1)/2 Z:brushPosition.z - (depth-1)/2];
        mask = [mask ROIMaskByUnioningWithMask:brushMask];
    }

    return mask;
}


@end
