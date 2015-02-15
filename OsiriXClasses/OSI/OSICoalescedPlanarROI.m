/*=========================================================================
 Program:   OsiriX
 
 Copyright (c) OsiriX Team
 All rights reserved.
 Distributed under GNU - LGPL
 
 See http://www.osirix-viewer.com/copyright.html for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.
 =========================================================================*/

#import "OSICoalescedPlanarROI.h"
#import "OSIFloatVolumeData.h"
#import "OSIROIMask.h"
#import "CPRGenerator.h"
#import "CPRGeneratorRequest.h"
#import "OSIGeometry.h"
#include <OpenGL/CGLMacro.h>

@interface OSICoalescedPlanarROI ()

@property (nonatomic, readonly, retain) OSIFloatVolumeData *coalescedROIMaskVolumeData;

- (NSData *)_maskRunsDataForSlab:(OSISlab)slab dicomToPixTransform:(N3AffineTransform)dicomToPixTransform minCorner:(N3VectorPointer)minCornerPtr;

@end

@implementation OSICoalescedPlanarROI

@synthesize sourceROIs = _sourceROIs;
@synthesize volumeTransform = _volumeTransform;

- (id)initWithSourceROIs:(NSArray *)rois
{
    if ( (self = [super init]) ) {
        _sourceROIs = [rois copy];
    }
    
    return self;
}

- (void)dealloc
{
    [_sourceROIs release];
    _sourceROIs = nil;
    [_coalescedROIMaskVolumeData release];
    _coalescedROIMaskVolumeData = nil;
    [_cachedMaskRunsData release];
    _cachedMaskRunsData = nil;
    
    [super dealloc];
}

- (NSString *)name
{
    if ([_sourceROIs count] > 0) {
        return [[_sourceROIs objectAtIndex:0] name];
    }
    return nil;
}

- (NSArray *)convexHull
{
    NSMutableArray *hull;
    OSIROI *roi;
    
    hull = [NSMutableArray array];
    
    for (roi in _sourceROIs) {
        [hull addObjectsFromArray:[roi convexHull]];
    }
    
    return hull;
}

- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume
{
    NSMutableArray *maskRuns;
    OSIFloatVolumeData *coalescedROIMaskVolume;
    CPRObliqueSliceGeneratorRequest *sliceRequest;
    CPRVolumeData *resampledVolume;
    OSIROIMask *resampledMask;
    OSIROI *roi;

    maskRuns = [NSMutableArray array];
    
    if (N3AffineTransformEqualToTransform(floatVolume.volumeTransform, self.volumeTransform)) {
        for (roi in _sourceROIs) {
            [maskRuns addObjectsFromArray:[[roi ROIMaskForFloatVolumeData:floatVolume] maskRuns]];
        }
        
        return [[[OSIROIMask alloc] initWithMaskRuns:maskRuns] autorelease];
    } else { // we need to create a 
        coalescedROIMaskVolume = self.coalescedROIMaskVolumeData;
        
        sliceRequest = [[[CPRObliqueSliceGeneratorRequest alloc] init] autorelease];
        sliceRequest.pixelsWide = floatVolume.pixelsWide;
        sliceRequest.pixelsHigh = floatVolume.pixelsHigh;
        sliceRequest.slabSampleDistance = floatVolume.pixelSpacingZ;
        sliceRequest.slabWidth = floatVolume.pixelSpacingZ * floatVolume.pixelsDeep;
        sliceRequest.sliceToDicomTransform = N3AffineTransformInvert(N3AffineTransformConcat(floatVolume.volumeTransform, N3AffineTransformMakeTranslation(0, 0, (float)floatVolume.pixelsDeep/2.0)));
        
        resampledVolume = [CPRGenerator synchronousRequestVolume:sliceRequest volumeData:coalescedROIMaskVolume];
        
        assert(floatVolume.pixelsWide == resampledVolume.pixelsWide);
        assert(floatVolume.pixelsHigh == resampledVolume.pixelsHigh);
        assert(floatVolume.pixelsDeep == resampledVolume.pixelsDeep);
        
        resampledMask = [OSIROIMask ROIMaskFromVolumeData:(OSIFloatVolumeData *)resampledVolume volumeTransform:NULL];
        
        return resampledMask;
    }
}

- (NSSet *)osiriXROIs
{
    OSIROI *roi;
    NSMutableSet *osirixROIs;
    
    osirixROIs = [NSMutableSet setWithCapacity:[_sourceROIs count]];
    
    for (roi in _sourceROIs) {
        [osirixROIs unionSet:[roi osiriXROIs]];
    }
    
    return osirixROIs;
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
    
    if (_cachedMaskRunsData && OSISlabEqualTo(slab, _cachedSlab) && N3AffineTransformEqualToTransform(dicomToPixTransform, _cachedDicomToPixTransform)) {
        maskRunsData = _cachedMaskRunsData;
        minCorner = _cachedMinCorner;
    } else {
        [_cachedMaskRunsData release];
        _cachedMaskRunsData = [[self _maskRunsDataForSlab:slab dicomToPixTransform:dicomToPixTransform minCorner:&_cachedMinCorner] retain];
        _cachedSlab = slab;
        _cachedDicomToPixTransform = dicomToPixTransform;
        maskRunsData = _cachedMaskRunsData;
        minCorner = _cachedMinCorner;
    }


    glLineWidth(3.0);    

    glColor4f(1, 0, 0, .4);
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
}


- (OSIFloatVolumeData *)coalescedROIMaskVolumeData
{
    N3Vector minCorner;
    N3Vector maxCorner;
    OSIROI *roi;
    NSValue *hullPointValue;
    N3Vector hullPoint;
    NSInteger width;
    NSInteger height;
    NSInteger depth;
    NSInteger i;
    N3AffineTransform coalescedROIMaskVolumeTransform;
    float *coalescedROIMaskVolumeBytes;
    float *bytesPtr;
    OSIROIMask *coalescedMask;
    OSIROIMaskRun maskRun;
    NSValue *maskRunValue;


    if (_coalescedROIMaskVolumeData == nil) {
        minCorner = N3VectorMake(CGFLOAT_MAX, CGFLOAT_MAX, CGFLOAT_MAX);
        maxCorner = N3VectorMake(-CGFLOAT_MAX, -CGFLOAT_MAX, -CGFLOAT_MAX);
        
        for (roi in _sourceROIs) {
            for (hullPointValue in [roi convexHull]) {
                hullPoint = N3VectorApplyTransform([hullPointValue N3VectorValue], self.volumeTransform);
                
                minCorner.x = MIN(minCorner.x, hullPoint.x);
                minCorner.y = MIN(minCorner.y, hullPoint.y);
                minCorner.z = MIN(minCorner.z, hullPoint.z);
                maxCorner.x = MAX(maxCorner.x, hullPoint.x);
                maxCorner.y = MAX(maxCorner.y, hullPoint.y);
                maxCorner.z = MAX(maxCorner.z, hullPoint.z);
            }
        }
        
        minCorner.x = floor(minCorner.x) - 1;
        minCorner.y = floor(minCorner.y) - 1;
        minCorner.z = floor(minCorner.z) - 1;
        maxCorner.x = ceil(maxCorner.x) + 1;
        maxCorner.y = ceil(maxCorner.y) + 1;
        maxCorner.z = ceil(maxCorner.z) + 1;
        
        width = maxCorner.x - minCorner.x;
        height = maxCorner.y - minCorner.y;
        depth = maxCorner.z - minCorner.z;
        
        coalescedROIMaskVolumeTransform = N3AffineTransformConcat(self.volumeTransform, N3AffineTransformMakeTranslation(-minCorner.x, -minCorner.y, -minCorner.z));
        
        coalescedROIMaskVolumeBytes = malloc(width * height * depth * sizeof(float));
        memset(coalescedROIMaskVolumeBytes, 0, width * height * depth * sizeof(float));
        _coalescedROIMaskVolumeData = [[OSIFloatVolumeData alloc] initWithFloatBytesNoCopy:coalescedROIMaskVolumeBytes pixelsWide:width pixelsHigh:height pixelsDeep:depth volumeTransform:coalescedROIMaskVolumeTransform outOfBoundsValue:0 freeWhenDone:YES];
        
        for (roi in _sourceROIs) {
            coalescedMask = [roi ROIMaskForFloatVolumeData:_coalescedROIMaskVolumeData];
            assert([_coalescedROIMaskVolumeData checkDebugROIMask:coalescedMask]);
            
            for (maskRunValue in [coalescedMask maskRuns]) {
                maskRun = [maskRunValue OSIROIMaskRunValue];
                
                for (i = maskRun.widthRange.location; i < NSMaxRange(maskRun.widthRange); i++) {
                    bytesPtr = &(coalescedROIMaskVolumeBytes[maskRun.depthIndex * width * height + maskRun.heightIndex * width + i]);
                    *bytesPtr = MAX(*bytesPtr, maskRun.intensity);
                }
            }
        }

    }
    
    return _coalescedROIMaskVolumeData;
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
    coalescedVolumeMaskToPixTransform = N3AffineTransformConcat(N3AffineTransformInvert(self.coalescedROIMaskVolumeData.volumeTransform), dicomToPixTransform);
    
    // first off figure out where this float volume needs to be
    corner = N3VectorApplyTransform(N3VectorMake(0,                                          0,                                          0), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(self.coalescedROIMaskVolumeData.pixelsWide, 0,                                          0), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                                          self.coalescedROIMaskVolumeData.pixelsHigh, 0), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(self.coalescedROIMaskVolumeData.pixelsWide, self.coalescedROIMaskVolumeData.pixelsHigh, 0), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                                          0,                                          self.coalescedROIMaskVolumeData.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(self.coalescedROIMaskVolumeData.pixelsWide, 0,                                          self.coalescedROIMaskVolumeData.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(0,                                          self.coalescedROIMaskVolumeData.pixelsDeep, self.coalescedROIMaskVolumeData.pixelsDeep), coalescedVolumeMaskToPixTransform);
    minCorner.x = MIN(minCorner.x, corner.x); minCorner.y = MIN(minCorner.y, corner.y);
    maxCorner.x = MAX(maxCorner.x, corner.x); maxCorner.y = MAX(maxCorner.y, corner.y);
    corner = N3VectorApplyTransform(N3VectorMake(self.coalescedROIMaskVolumeData.pixelsWide, self.coalescedROIMaskVolumeData.pixelsDeep, self.coalescedROIMaskVolumeData.pixelsDeep), coalescedVolumeMaskToPixTransform);
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
    
    sliceRequest.sliceToDicomTransform = N3AffineTransformConcat(N3AffineTransformMakeTranslation(minCorner.x, minCorner.y, 0), N3AffineTransformInvert(dicomToPixTransform));
    
    floatVolumeData = [CPRGenerator synchronousRequestVolume:sliceRequest volumeData:self.coalescedROIMaskVolumeData];
    sliceMask = [OSIROIMask ROIMaskFromVolumeData:(OSIFloatVolumeData *)floatVolumeData volumeTransform:NULL];
    
    if (minCornerPtr) {
        *minCornerPtr = minCorner;
    }
    
    return [sliceMask maskRunsData];
 }

@end
























