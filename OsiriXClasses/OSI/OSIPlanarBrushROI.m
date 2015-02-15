//
//  OSIPlanarBrushROI.m
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 9/26/11.
//  Copyright 2011 OsiriX Team. All rights reserved.
//

#import "OSIPlanarBrushROI.h"
#import "OSIROI+Private.h"
#import "OSIROIMask.h"
#import "ROI.h"
#import "DCMView.h"
#import "OSIGeometry.h"
#import "OSIFloatVolumeData.h"
#import "CPRGenerator.h"
#import "CPRGeneratorRequest.h"
#import "CPRVolumeData.h"

@interface OSIPlanarBrushROI ()

@end

@implementation OSIPlanarBrushROI (Private)

- (id)initWithOsiriXROI:(ROI *)roi pixToDICOMTransfrom:(N3AffineTransform)pixToDICOMTransfrom
{
	NSMutableArray *hullPoints;
    NSInteger i;
    NSInteger j;
    N3AffineTransform volumeTransform;
    float* mask;
	
	if ( (self = [super init]) ) {
		_osiriXROI = [roi retain];
		
		_plane = N3PlaneApplyTransform(N3PlaneZZero, pixToDICOMTransfrom);

		if ([roi type] == tPlain) {
            hullPoints = [[NSMutableArray alloc] init];
            
            [hullPoints addObject:[NSValue valueWithN3Vector:N3VectorApplyTransform(N3VectorMake(roi.textureUpLeftCornerX, roi.textureUpLeftCornerY, 0), pixToDICOMTransfrom)]];
            [hullPoints addObject:[NSValue valueWithN3Vector:N3VectorApplyTransform(N3VectorMake(roi.textureUpLeftCornerX, roi.textureDownRightCornerY, 0), pixToDICOMTransfrom)]];
            [hullPoints addObject:[NSValue valueWithN3Vector:N3VectorApplyTransform(N3VectorMake(roi.textureDownRightCornerX, roi.textureUpLeftCornerY, 0), pixToDICOMTransfrom)]];
            [hullPoints addObject:[NSValue valueWithN3Vector:N3VectorApplyTransform(N3VectorMake(roi.textureDownRightCornerX, roi.textureDownRightCornerY, 0), pixToDICOMTransfrom)]];
            _convexHull = hullPoints;
            
            mask = malloc(roi.textureWidth * roi.textureHeight * sizeof(float));
            memset(mask, 0, roi.textureWidth * roi.textureHeight * sizeof(float));
            
            for (j = 0; j < roi.textureHeight; j++) {
                for (i = 0; i < roi.textureWidth; i++) {
                    mask[j*roi.textureWidth + i] = ((float)roi.textureBuffer[j*roi.textureWidth + i])/255.0;
                }
            }
            volumeTransform = N3AffineTransformConcat(N3AffineTransformInvert(pixToDICOMTransfrom), N3AffineTransformMakeTranslation(-roi.textureUpLeftCornerX, -roi.textureUpLeftCornerY, 0));
            _brushMask = [[OSIFloatVolumeData alloc] initWithFloatBytesNoCopy:mask pixelsWide:roi.textureWidth pixelsHigh:roi.textureHeight pixelsDeep:1 volumeTransform:volumeTransform outOfBoundsValue:0 freeWhenDone:YES];
        } else {
			[self autorelease];
			self = nil;
		}
	}
	return self;
}


@end


@implementation OSIPlanarBrushROI

- (void)dealloc
{
    [_osiriXROI release];
    _osiriXROI = nil;
    [_brushMask release];
    _brushMask = nil;
    [_convexHull release];
    _convexHull = nil;
    
    [super dealloc];
}

- (NSString *)name
{
	return [_osiriXROI name];
}

- (NSArray *)convexHull
{
    return _convexHull;
}

- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume
{
    CPRObliqueSliceGeneratorRequest *request;
    CPRVolumeData *volume;
    OSIROIMask *roiMask;
    N3AffineTransform sliceToDicomTransform;
    N3Vector planePixelPoint;
    
    // make sure floatVolume's z direction is perpendicular to the plane
    assert(N3VectorLength(N3VectorCrossProduct(N3VectorApplyTransformToDirectionalVector(_plane.normal, floatVolume.volumeTransform), N3VectorMake(0, 0, 1))) < 0.0001);
    
    planePixelPoint = N3VectorApplyTransform(_plane.point, floatVolume.volumeTransform);
    sliceToDicomTransform = N3AffineTransformInvert(N3AffineTransformConcat([floatVolume volumeTransform], N3AffineTransformMakeTranslation(0, 0, -planePixelPoint.z)));
    
    request = [[[CPRObliqueSliceGeneratorRequest alloc] init] autorelease];
    request.sliceToDicomTransform = sliceToDicomTransform;
    request.pixelsWide = floatVolume.pixelsWide;
    request.pixelsHigh = floatVolume.pixelsHigh;
    request.interpolationMode = CPRInterpolationModeNearestNeighbor;
    
    volume = [CPRGenerator synchronousRequestVolume:request volumeData:_brushMask];    
    roiMask = [OSIROIMask ROIMaskFromVolumeData:(OSIFloatVolumeData *)volume volumeTransform:NULL];
    
#if CGFLOAT_IS_DOUBLE
    return [roiMask ROIMaskByTranslatingByX:0 Y:0 Z:round(planePixelPoint.z)];
#else
    return [roiMask ROIMaskByTranslatingByX:0 Y:0 Z:roundf(planePixelPoint.z)];
#endif

}

- (void)drawRect:(NSRect)rect inSlab:(OSISlab)slab inCGLContext:(CGLContextObj)cgl_ctx pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform;
{
	double dicomToPixGLTransform[16];
	
	if (OSISlabContainsPlane(slab, _plane) == NO) {
		return; // this ROI does not live on this slice
	}
    
    N3AffineTransformGetOpenGLMatrixd(dicomToPixTransform, dicomToPixGLTransform);
	
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glMultMatrixd(dicomToPixGLTransform);
    
    glLineWidth(3.0);    
    
    // let's try drawing some the mask
    OSIROIMask *mask;
    NSArray *maskRuns;
    OSIROIMaskRun maskRun;
    NSValue *maskRunValue;
    N3AffineTransform inverseVolumeTransform;
    N3Vector lineStart;
    N3Vector lineEnd;

    // dicom to pix transform gives us the actual pixel spacing to use
    // it might be interesting to make a float volume data that only initializes its memory
    // when it is really needed, that might allow a function like this to
    // but for now make a float volume data that is the size of what is being visualized.

    // rect will give the width and height, but we will use slab to find the depth
    // to slab gives us the height in cm, and the dicomToPixTransform will give us the pixelspacing

    // all this, but for now, ignore the thickness
//    N3AffineTransform pixToDicomTransform = N3AffineTransformInvert(dicomToPixTransform);
//    CGFloat cmPerPix = N3VectorLength(N3VectorApplyTransform(N3VectorZBasis, pixToDicomTransform));
//    CGFloat slabThicknessPix = ceil(slab.thickness / cmPerPix);
    N3AffineTransform dicomToFloatDataTransform = N3AffineTransformConcat(dicomToPixTransform, N3AffineTransformMakeTranslation(-rect.origin.x, -rect.origin.y, 0));
    NSData *floatData = [NSData dataWithBytes:NULL length:rect.size.width * rect.size.height];
    OSIFloatVolumeData *volumeData = [[OSIFloatVolumeData alloc] initWithData:floatData pixelsWide:rect.size.width pixelsHigh:rect.size.height pixelsDeep:1
                                                              volumeTransform:dicomToFloatDataTransform outOfBoundsValue:0];
    
    inverseVolumeTransform = N3AffineTransformInvert([volumeData volumeTransform]);
    mask = [self ROIMaskForFloatVolumeData:volumeData];
    maskRuns = [mask maskRuns];
    
    glColor3f(1, 0, 1);
    glBegin(GL_LINES);
    for (maskRunValue in maskRuns) {
        maskRun = [maskRunValue OSIROIMaskRunValue];
        
        lineStart = N3VectorMake(maskRun.widthRange.location, maskRun.heightIndex + 0.5, maskRun.depthIndex);
        lineEnd = N3VectorMake(NSMaxRange(maskRun.widthRange), maskRun.heightIndex + 0.5, maskRun.depthIndex);
        
        lineStart = N3VectorApplyTransform(lineStart, inverseVolumeTransform);
        lineEnd = N3VectorApplyTransform(lineEnd, inverseVolumeTransform);
        
        glVertex3d(lineStart.x, lineStart.y, lineStart.z);
        glVertex3d(lineEnd.x, lineEnd.y, lineEnd.z);
    }
    glEnd();
    
    glPopMatrix();
}


- (NSSet *)osiriXROIs
{
	return [NSSet setWithObject:_osiriXROI];
}



@end




















