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

#import "CPRVolumeData.h"
#import "DCMPix.h"
#import "CPRUnsignedInt16ImageRep.h"
#import "CPRGenerator.h"
#import "CPRGeneratorRequest.h"

@interface CPRVolumeData ()

- (BOOL)_testOrientationMatrix:(double[9])orientation; // returns YES if the orientation matrix's determinant is non-zero

@property (nonatomic, readonly, assign) float* floatBytes;

@end


@implementation CPRVolumeData

@synthesize outOfBoundsValue = _outOfBoundsValue;
@synthesize pixelsWide = _pixelsWide;
@synthesize pixelsHigh = _pixelsHigh;
@synthesize pixelsDeep = _pixelsDeep;
@synthesize volumeTransform = _volumeTransform;

- (instancetype)initWithFloatBytesNoCopy:(const float *)floatBytes pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
                         volumeTransform:(N3AffineTransform)volumeTransform outOfBoundsValue:(float)outOfBoundsValue freeWhenDone:(BOOL)freeWhenDone // volumeTransform is the transform from Dicom (patient) space to pixel data
{
    return [self initWithData:[NSData dataWithBytesNoCopy:(void *)floatBytes length:sizeof(float) * pixelsWide * pixelsHigh * pixelsDeep freeWhenDone:freeWhenDone]
                   pixelsWide:pixelsWide pixelsHigh:pixelsHigh pixelsDeep:pixelsDeep volumeTransform:volumeTransform outOfBoundsValue:outOfBoundsValue];
}

- (instancetype)initWithData:(NSData *)data pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
             volumeTransform:(N3AffineTransform)volumeTransform outOfBoundsValue:(float)outOfBoundsValue // volumeTransform is the transform from Dicom (patient) space to pixel data
{
    if ( (self = [super init]) ) {
        _floatData = [data retain];
        _outOfBoundsValue = outOfBoundsValue;
        _pixelsWide = pixelsWide;
        _pixelsHigh = pixelsHigh;
        _pixelsDeep = pixelsDeep;
        _volumeTransform = volumeTransform;
    }
    return self;
}

- (instancetype)initWithVolumeData:(CPRVolumeData *)volumeData
{
    if ( (self = [super init]) ) {
        _floatData = [volumeData->_floatData retain];
        _outOfBoundsValue = volumeData->_outOfBoundsValue;
        _pixelsWide = volumeData->_pixelsWide;
        _pixelsHigh = volumeData->_pixelsHigh;
        _pixelsDeep = volumeData->_pixelsDeep;
        _volumeTransform = volumeData->_volumeTransform;
    }
    return self;
}

- (void)dealloc
{
    [_floatData release];
    _floatData = nil;
    [super dealloc];
}

- (BOOL)isRectilinear
{
    return N3AffineTransformIsRectilinear(_volumeTransform);
}

- (CGFloat)minPixelSpacing
{
    return MIN(MIN(self.pixelSpacingX, self.pixelSpacingY), self.pixelSpacingZ);
}

- (CGFloat)pixelSpacingX
{
    N3Vector zero;
    N3AffineTransform inverseTransform;

    if (self.rectilinear) {
        return 1.0/_volumeTransform.m11;
    } else {
        inverseTransform = N3AffineTransformInvert(_volumeTransform);
        zero = N3VectorApplyTransform(N3VectorZero, inverseTransform);
        return N3VectorDistance(zero, N3VectorApplyTransform(N3VectorMake(1.0, 0.0, 0.0), inverseTransform));
    }
}

- (CGFloat)pixelSpacingY
{
    N3Vector zero;
    N3AffineTransform inverseTransform;

    if (self.rectilinear) {
        return 1.0/_volumeTransform.m22;
    } else {
        inverseTransform = N3AffineTransformInvert(_volumeTransform);
        zero = N3VectorApplyTransform(N3VectorZero, inverseTransform);
        return N3VectorDistance(zero, N3VectorApplyTransform(N3VectorMake(0.0, 1.0, 0.0), inverseTransform));
    }
}

- (CGFloat)pixelSpacingZ
{
    N3Vector zero;
    N3AffineTransform inverseTransform;

    if (self.rectilinear) {
        return 1.0/_volumeTransform.m33;
    } else {
        inverseTransform = N3AffineTransformInvert(_volumeTransform);
        zero = N3VectorApplyTransform(N3VectorZero, inverseTransform);
        return N3VectorDistance(zero, N3VectorApplyTransform(N3VectorMake(0.0, 0.0, 1.0), inverseTransform));
    }
}

- (BOOL)isDataValid __deprecated
{
    return YES;
}

- (void)invalidateData __deprecated {}

- (float *)floatBytes
{
    return (float *)[_floatData bytes];
}

// will copy fill length*sizeof(float) bytes
- (BOOL)getFloatRun:(float *)buffer atPixelCoordinateX:(NSUInteger)x y:(NSUInteger)y z:(NSUInteger)z length:(NSUInteger)length
{
    assert(x < _pixelsWide);
    assert(y < _pixelsHigh);
    assert(z < _pixelsDeep);
    assert(x + length < _pixelsWide);
    
    memcpy(buffer, &(self.floatBytes[x + y*_pixelsWide + z*_pixelsWide*_pixelsHigh]), length * sizeof(float));
    return YES;
}

- (CPRUnsignedInt16ImageRep *)unsignedInt16ImageRepForSliceAtIndex:(NSUInteger)z
{
    CPRUnsignedInt16ImageRep *imageRep;
    uint16_t *unsignedInt16Data;
    vImage_Buffer floatBuffer;
    vImage_Buffer unsignedInt16Buffer;
    
    imageRep = [[CPRUnsignedInt16ImageRep alloc] initWithData:NULL pixelsWide:_pixelsWide pixelsHigh:_pixelsHigh];
    imageRep.pixelSpacingX = [self pixelSpacingX];
    imageRep.pixelSpacingY = [self pixelSpacingY];
    imageRep.sliceThickness = [self pixelSpacingZ];
    imageRep.imageToDicomTransform = N3AffineTransformConcat(N3AffineTransformMakeTranslation(0.0, 0.0, (CGFloat)z), N3AffineTransformInvert(_volumeTransform));
    
    unsignedInt16Data = [imageRep unsignedInt16Data];
    
    floatBuffer.data = (void *)self.floatBytes + (_pixelsWide * _pixelsHigh * sizeof(float) * z);
    floatBuffer.height = _pixelsHigh;
    floatBuffer.width = _pixelsWide;
    floatBuffer.rowBytes = sizeof(float) * _pixelsWide;
    
    unsignedInt16Buffer.data = unsignedInt16Data;
    unsignedInt16Buffer.height = _pixelsHigh;
    unsignedInt16Buffer.width = _pixelsWide;
    unsignedInt16Buffer.rowBytes = sizeof(uint16_t) * _pixelsWide;
    
    vImageConvert_FTo16U(&floatBuffer, &unsignedInt16Buffer, -1024, 1, 0);
    imageRep.slope = 1;
    imageRep.offset = -1024;
    
    return [imageRep autorelease];
}

- (CPRVolumeData *)volumeDataForSliceAtIndex:(NSUInteger)z
{
    CPRVolumeData *sliceVolume;
    N3AffineTransform sliceVolumeTransform;
    NSData *sliceData;

    sliceVolumeTransform = N3AffineTransformConcat(_volumeTransform, N3AffineTransformMakeTranslation(0, 0, -z));

    sliceData = [NSData dataWithBytes:self.floatBytes + (_pixelsWide*_pixelsHigh*z) length:_pixelsWide * _pixelsHigh * sizeof(float)];
    sliceVolume = [[[self class] alloc] initWithData:sliceData pixelsWide:_pixelsWide pixelsHigh:_pixelsHigh pixelsDeep:1 volumeTransform:sliceVolumeTransform outOfBoundsValue:_outOfBoundsValue];

    return [sliceVolume autorelease];
}

- (BOOL)getFloat:(float *)floatPtr atPixelCoordinateX:(NSUInteger)x y:(NSUInteger)y z:(NSUInteger)z
{
    CPRVolumeDataInlineBuffer inlineBuffer;
    
    [self aquireInlineBuffer:&inlineBuffer];
    *floatPtr = CPRVolumeDataGetFloatAtPixelCoordinate(&inlineBuffer, x, y, z);
    return YES;
}

- (BOOL)getLinearInterpolatedFloat:(float *)floatPtr atDicomVector:(N3Vector)vector
{
    CPRVolumeDataInlineBuffer inlineBuffer;
    [self aquireInlineBuffer:&inlineBuffer];
    *floatPtr = CPRVolumeDataLinearInterpolatedFloatAtDicomVector(&inlineBuffer, vector);
    return YES;
}

- (BOOL)getNearestNeighborInterpolatedFloat:(float *)floatPtr atDicomVector:(N3Vector)vector
{
    CPRVolumeDataInlineBuffer inlineBuffer;
    
    [self aquireInlineBuffer:&inlineBuffer];
    *floatPtr = CPRVolumeDataNearestNeighborInterpolatedFloatAtDicomVector(&inlineBuffer, vector);
    return YES;
}

- (BOOL)getCubicInterpolatedFloat:(float *)floatPtr atDicomVector:(N3Vector)vector
{
    CPRVolumeDataInlineBuffer inlineBuffer;

    [self aquireInlineBuffer:&inlineBuffer];
    *floatPtr = CPRVolumeDataCubicInterpolatedFloatAtDicomVector(&inlineBuffer, vector);
    return YES;
}


- (instancetype)volumeDataByApplyingTransform:(N3AffineTransform)transform;
{
    return [[[CPRVolumeData alloc] initWithData:_floatData pixelsWide:_pixelsWide pixelsHigh:_pixelsHigh pixelsDeep:_pixelsDeep
                                                                         volumeTransform:N3AffineTransformConcat(_volumeTransform, transform) outOfBoundsValue:_outOfBoundsValue] autorelease];
}

- (instancetype)volumeDataResampledWithVolumeTransform:(N3AffineTransform)transform interpolationMode:(CPRInterpolationMode)interpolationsMode
{
    if (N3AffineTransformEqualToTransform(_volumeTransform, transform)) {
        return self;
    }

    N3AffineTransform oringinalVoxelToDicomTransform = N3AffineTransformInvert(_volumeTransform);
    N3AffineTransform originalVoxelToNewVoxelTransform = N3AffineTransformConcat(oringinalVoxelToDicomTransform, transform);

    N3Vector minCorner = N3VectorZero;
    N3Vector maxCorner = N3VectorZero;

    N3Vector corner1 = N3VectorApplyTransform(N3VectorMake(0,           0,           0), originalVoxelToNewVoxelTransform);
    N3Vector corner2 = N3VectorApplyTransform(N3VectorMake(_pixelsWide, 0,           0), originalVoxelToNewVoxelTransform);
    N3Vector corner3 = N3VectorApplyTransform(N3VectorMake(0,           _pixelsHigh, 0), originalVoxelToNewVoxelTransform);
    N3Vector corner4 = N3VectorApplyTransform(N3VectorMake(_pixelsWide, _pixelsHigh, 0), originalVoxelToNewVoxelTransform);
    N3Vector corner5 = N3VectorApplyTransform(N3VectorMake(0,           0,           _pixelsDeep), originalVoxelToNewVoxelTransform);
    N3Vector corner6 = N3VectorApplyTransform(N3VectorMake(_pixelsWide, 0,           _pixelsDeep), originalVoxelToNewVoxelTransform);
    N3Vector corner7 = N3VectorApplyTransform(N3VectorMake(0,           _pixelsHigh, _pixelsDeep), originalVoxelToNewVoxelTransform);
    N3Vector corner8 = N3VectorApplyTransform(N3VectorMake(_pixelsWide, _pixelsHigh, _pixelsDeep), originalVoxelToNewVoxelTransform);

    minCorner.x = MIN(MIN(MIN(MIN(MIN(MIN(MIN(corner1.x, corner2.x), corner3.x), corner4.x), corner5.x), corner6.x), corner7.x), corner8.x);
    minCorner.y = MIN(MIN(MIN(MIN(MIN(MIN(MIN(corner1.y, corner2.y), corner3.y), corner4.y), corner5.y), corner6.y), corner7.y), corner8.y);
    minCorner.z = MIN(MIN(MIN(MIN(MIN(MIN(MIN(corner1.z, corner2.z), corner3.z), corner4.z), corner5.z), corner6.z), corner7.z), corner8.z);
    maxCorner.x = MAX(MAX(MAX(MAX(MAX(MAX(MAX(corner1.x, corner2.x), corner3.x), corner4.x), corner5.x), corner6.x), corner7.x), corner8.x);
    maxCorner.y = MAX(MAX(MAX(MAX(MAX(MAX(MAX(corner1.y, corner2.y), corner3.y), corner4.y), corner5.y), corner6.y), corner7.y), corner8.y);
    maxCorner.z = MAX(MAX(MAX(MAX(MAX(MAX(MAX(corner1.z, corner2.z), corner3.z), corner4.z), corner5.z), corner6.z), corner7.z), corner8.z);

#if CGFLOAT_IS_DOUBLE
    minCorner.x = floor(minCorner.x);
    minCorner.y = floor(minCorner.y);
    minCorner.z = floor(minCorner.z);
    maxCorner.x = ceil(maxCorner.x);
    maxCorner.y = ceil(maxCorner.y);
    maxCorner.z = ceil(maxCorner.z);
#else
    minCorner.x = floorf(minCorner.x);
    minCorner.y = floorf(minCorner.y);
    minCorner.z = floorf(minCorner.z);
    maxCorner.x = ceilf(maxCorner.x);
    maxCorner.y = ceilf(maxCorner.y);
    maxCorner.z = ceilf(maxCorner.z);
#endif

    NSUInteger width = (NSUInteger)(maxCorner.x - minCorner.x) + 1;
    NSUInteger height = (NSUInteger)(maxCorner.y - minCorner.y) + 1;
    NSUInteger depth = (NSUInteger)(maxCorner.z - minCorner.z) + 1;

    N3AffineTransform shiftedTransform = N3AffineTransformConcat(transform, N3AffineTransformMakeTranslation(-1.0*(CGFloat)minCorner.x, -1.0*(CGFloat)minCorner.y, -1.0*(CGFloat)minCorner.z));

    return [self volumeDataResampledWithVolumeTransform:shiftedTransform pixelsWide:width pixelsHigh:height pixelsDeep:depth interpolationMode:interpolationsMode];
}

- (instancetype)volumeDataResampledWithVolumeTransform:(N3AffineTransform)transform pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
                                        interpolationMode:(CPRInterpolationMode)interpolationsMode;
{
    N3AffineTransform newVoxelToDicomTransform = N3AffineTransformInvert(transform);

    CPRObliqueSliceGeneratorRequest *request = [[[CPRObliqueSliceGeneratorRequest alloc] init] autorelease];

    N3Vector xBasis = N3VectorApplyTransformToDirectionalVector(N3VectorXBasis, newVoxelToDicomTransform);
    N3Vector yBasis = N3VectorApplyTransformToDirectionalVector(N3VectorYBasis, newVoxelToDicomTransform);
    N3Vector zBasis = N3VectorApplyTransformToDirectionalVector(N3VectorZBasis, newVoxelToDicomTransform);

    request.origin = N3VectorApplyTransform(N3VectorZero, newVoxelToDicomTransform);
    request.directionX = N3VectorNormalize(xBasis);
    request.directionY = N3VectorNormalize(yBasis);
    request.directionZ = N3VectorNormalize(zBasis);

    request.pixelSpacingX = N3VectorLength(xBasis);
    request.pixelSpacingY = N3VectorLength(yBasis);
    request.pixelSpacingZ = N3VectorLength(zBasis);

    request.pixelsWide = pixelsWide;
    request.pixelsHigh = pixelsHigh;
    request.slabWidth = N3VectorLength(zBasis) * (CGFloat)pixelsDeep;

    request.interpolationMode = interpolationsMode;

    return [[[[self class] alloc] initWithVolumeData:[CPRGenerator synchronousRequestVolume:request volumeData:self]] autorelease]; // The little extra init is to make sure we return the correct class (as in OSIFloatVolumeData)
}

- (NSUInteger)tempBufferSizeForNumVectors:(NSUInteger)numVectors
{
    return numVectors * sizeof(float) * 11;
}

- (BOOL)isEqual:(id)object
{
    BOOL isEqual = NO;
    if ([object isKindOfClass:[CPRVolumeData class]]) {
        CPRVolumeDataInlineBuffer inlineBuffer1;
        CPRVolumeDataInlineBuffer inlineBuffer2;
        CPRVolumeData *otherVolumeData = (CPRVolumeData *)object;

        [self aquireInlineBuffer:&inlineBuffer1];
        [otherVolumeData aquireInlineBuffer:&inlineBuffer2];

        if (inlineBuffer1.floatBytes == inlineBuffer2.floatBytes &&
            inlineBuffer1.outOfBoundsValue == inlineBuffer2.outOfBoundsValue &&
            inlineBuffer1.pixelsWide == inlineBuffer2.pixelsWide &&
            inlineBuffer1.pixelsHigh == inlineBuffer2.pixelsHigh &&
            inlineBuffer1.pixelsDeep == inlineBuffer2.pixelsDeep &&
            N3AffineTransformEqualToTransform(inlineBuffer1.volumeTransform, inlineBuffer2.volumeTransform)) {
            isEqual = YES;
        }
    }

    return isEqual;
}

// not done yet, will crash if given vectors that are outside of the volume
- (void)linearInterpolateVolumeVectors:(N3VectorArray)volumeVectors outputValues:(float *)outputValues numVectors:(NSUInteger)numVectors tempBuffer:(void *)tempBuffer
{
    float *interpolateBuffer = (float *)tempBuffer;
    
    float *scrap;
    
    float *yFloors;
    float *zFloors;
    float *yFrac;
    float *zFrac;
    float *yNegFrac;
    float *zNegFrac;
    float negOne;
    
    float *i1Positions;
    float *i2Positions;
    float *j1Positions;
    float *j2Positions;
    
    float *yFloorPosition;
    float *yCielPosition;
    float *zFloorPosition;
    float *zCielPosition;

    float width;
    float widthTimesHeight;
    float widthTimesHeightTimeDepth;
    
    float *i1;
    float *i2;
    float *j1;
    float *j2;
    
    float *w1;
    float *w2;

    negOne = -1.0;
    width = _pixelsWide;
    widthTimesHeight = _pixelsWide*_pixelsHigh;
    widthTimesHeightTimeDepth = widthTimesHeight*_pixelsDeep;
    
    yFrac = interpolateBuffer + (numVectors * 0);
    zFrac = interpolateBuffer + (numVectors * 1);
    yFloors = interpolateBuffer + (numVectors * 2);
    zFloors = interpolateBuffer + (numVectors * 3);
    
    vDSP_vfrac(((float *)volumeVectors) + 1, 3, yFrac, 1, numVectors);
    vDSP_vsub(yFrac, 1, ((float *)volumeVectors) + 1, 3, yFloors, 1, numVectors);

    vDSP_vfrac(((float *)volumeVectors) + 2, 3, zFrac, 1, numVectors);
    vDSP_vsub(zFrac, 1, ((float *)volumeVectors) + 2, 3, zFloors, 1, numVectors);
    
    yFloorPosition = interpolateBuffer + (numVectors * 6);
    yCielPosition = interpolateBuffer + (numVectors * 7);
    zFloorPosition = interpolateBuffer + (numVectors * 8);
    zCielPosition = interpolateBuffer + (numVectors * 9);
        
    vDSP_vsmul(yFloors, 1, &width, yFloorPosition, 1, numVectors);
    vDSP_vsadd(yFloorPosition, 1, &width, yCielPosition, 1, numVectors);
    vDSP_vsmul(zFloors, 1, &widthTimesHeight, zFloorPosition, 1, numVectors);
    vDSP_vsadd(zFloorPosition, 1, &widthTimesHeight, zCielPosition, 1, numVectors);
    
    i1Positions = interpolateBuffer + (numVectors * 2);
    i2Positions = interpolateBuffer + (numVectors * 3);
    j1Positions = interpolateBuffer + (numVectors * 4);
    j2Positions = interpolateBuffer + (numVectors * 5);
    
    // i1 yFloor zFloor
    // i2 yFloor zCiel
    // j1 yCiel zFloor
    // j2 yCiel zCiel
    
    vDSP_vadd((float *)volumeVectors, 3, yFloorPosition, 1, i1Positions, 1, numVectors);
    
    vDSP_vadd(i1Positions, 1, zCielPosition, 1, i2Positions, 1, numVectors);
    vDSP_vadd(i1Positions, 1, zFloorPosition, 1, i1Positions, 1, numVectors);
    
    vDSP_vadd((float *)volumeVectors, 3, yCielPosition, 1, j1Positions, 1, numVectors);
    
    vDSP_vadd(j1Positions, 1, zCielPosition, 1, j2Positions, 1, numVectors);
    vDSP_vadd(j1Positions, 1, zFloorPosition, 1, j1Positions, 1, numVectors);
    
    
    i1 = interpolateBuffer + (numVectors * 6);
    i2 = interpolateBuffer + (numVectors * 7);
    j1 = interpolateBuffer + (numVectors * 8);
    j2 = interpolateBuffer + (numVectors * 9);
    
    vDSP_vlint((float *)self.floatBytes, i1Positions, 1, i1, 1, numVectors * 4, widthTimesHeightTimeDepth);
                
    yNegFrac = interpolateBuffer + (numVectors * 2);
    zNegFrac = interpolateBuffer + (numVectors * 3);
    
    vDSP_vsadd(yFrac, 1, &negOne, yNegFrac, 1, numVectors);
    vDSP_vneg(yNegFrac, 1, yNegFrac, 1, numVectors);
    
    vDSP_vsadd(zFrac, 1, &negOne, zNegFrac, 1, numVectors);
    vDSP_vneg(zNegFrac, 1, zNegFrac, 1, numVectors);
        
    w1 = interpolateBuffer + (numVectors * 4);
    w2 = interpolateBuffer + (numVectors * 5);
    
    scrap = interpolateBuffer + (numVectors * 10);
    
    vDSP_vmul(i1, 1, zNegFrac, 1, w1, 1, numVectors);
    vDSP_vmul(i2, 1, zFrac, 1, scrap, 1, numVectors);
    vDSP_vadd(w1, 1, scrap, 1, w1, 1, numVectors);
    
    vDSP_vmul(j1, 1, zNegFrac, 1, w2, 1, numVectors);
    vDSP_vmul(j2, 1, zFrac, 1, scrap, 1, numVectors);
    vDSP_vadd(w2, 1, scrap, 1, w2, 1, numVectors);
    
    
    vDSP_vmul(w1, 1, yNegFrac, 1, outputValues, 1, numVectors);
    vDSP_vmul(w2, 1, yFrac, 1, scrap, 1, numVectors);
    vDSP_vadd(outputValues, 1, scrap, 1, outputValues, 1, numVectors);
}

        

- (BOOL)aquireInlineBuffer:(CPRVolumeDataInlineBuffer *)inlineBuffer
{
    memset(inlineBuffer, 0, sizeof(CPRVolumeDataInlineBuffer));
    inlineBuffer->floatBytes = self.floatBytes;
    inlineBuffer->outOfBoundsValue = _outOfBoundsValue;
    inlineBuffer->pixelsWide = _pixelsWide;
    inlineBuffer->pixelsHigh = _pixelsHigh;
    inlineBuffer->pixelsDeep = _pixelsDeep;
    inlineBuffer->pixelsWideTimesPixelsHigh = _pixelsWide*_pixelsHigh;
    inlineBuffer->volumeTransform = _volumeTransform;
    return YES;
}

- (void)releaseInlineBuffer:(CPRVolumeDataInlineBuffer *)inlineBuffer __deprecated {}

- (BOOL)_testOrientationMatrix:(double[9])orientation // returns YES if the orientation matrix's determinant is non-zero
{
    N3AffineTransform transform;
    
    transform = N3AffineTransformIdentity;
    transform.m11 = orientation[0];
    transform.m12 = orientation[1];
    transform.m13 = orientation[2];
    transform.m21 = orientation[3];
    transform.m22 = orientation[4];
    transform.m23 = orientation[5];
    transform.m31 = orientation[6];
    transform.m32 = orientation[7];
    transform.m33 = orientation[8];
    
    return N3AffineTransformDeterminant(transform) != 0.0;
}

- (NSString *)description
{
	NSMutableString *description = [NSMutableString string];
	[description appendString:NSStringFromClass([self class])];
	[description appendString:[NSString stringWithFormat: @"\nPixels Wide: %lld\n", (long long)_pixelsWide]];
	[description appendString:[NSString stringWithFormat: @"Pixels High: %lld\n", (long long)_pixelsHigh]];
	[description appendString:[NSString stringWithFormat: @"Pixels Deep: %lld\n", (long long)_pixelsDeep]];
	[description appendString:[NSString stringWithFormat: @"Out of Bounds Value: %f\n", _outOfBoundsValue]];
	[description appendString:[NSString stringWithFormat: @"Volume Transform:\n%@\n", NSStringFromN3AffineTransform(_volumeTransform)]];

	return description;
}

@end

@implementation CPRVolumeData (DCMPixAndVolume)

- (id)initWithWithPixList:(NSArray *)pixList volume:(NSData *)volume
{
    DCMPix *firstPix;
    float sliceThickness;
    N3AffineTransform pixToDicomTransform;
    double spacingX;
    double spacingY;
    double spacingZ;
    double orientation[9];
	
    firstPix = [pixList objectAtIndex:0];
    
    sliceThickness = [firstPix sliceInterval];
	if(sliceThickness == 0)
	{
		NSLog(@"slice interval = slice thickness!");
		sliceThickness = [firstPix sliceThickness];
	}
        
    memset(orientation, 0, sizeof(double) * 9);
    [firstPix orientationDouble:orientation];
    spacingX = firstPix.pixelSpacingX;
    spacingY = firstPix.pixelSpacingY;
    if(sliceThickness == 0) { // if the slice thickness is still 0, make it the same as the average of the spacingX and spacingY
        sliceThickness = (spacingX + spacingY)/2.0;
    }
    spacingZ = sliceThickness;
    
    // test to make sure that orientation is initialized, when the volume is curved or something, it doesn't make sense to talk about orientation, and
    // so the orientation is really bogus
    // the test we will do is to make sure that orientation is 3 non-degenerate vectors
    if ([self _testOrientationMatrix:orientation] == NO) {
        memset(orientation, 0, sizeof(double)*9);
        orientation[0] = orientation[4] = orientation[8] = 1;
    }

    // This is not really the pixToDicom because for the CPRVolumeData uses Center rule whereas DCMPix uses Top-Left rule.
    pixToDicomTransform = N3AffineTransformIdentity;
    pixToDicomTransform.m41 = firstPix.originX;
    pixToDicomTransform.m42 = firstPix.originY;
    pixToDicomTransform.m43 = firstPix.originZ;
    pixToDicomTransform.m11 = orientation[0]*spacingX;
    pixToDicomTransform.m12 = orientation[1]*spacingX;
    pixToDicomTransform.m13 = orientation[2]*spacingX;
    pixToDicomTransform.m21 = orientation[3]*spacingY;
    pixToDicomTransform.m22 = orientation[4]*spacingY;
    pixToDicomTransform.m23 = orientation[5]*spacingY;
    pixToDicomTransform.m31 = orientation[6]*spacingZ;
    pixToDicomTransform.m32 = orientation[7]*spacingZ;
    pixToDicomTransform.m33 = orientation[8]*spacingZ;

    self = [self initWithData:volume pixelsWide:[firstPix pwidth] pixelsHigh:[firstPix pheight] pixelsDeep:[pixList count]
                          volumeTransform:N3AffineTransformInvert(pixToDicomTransform) outOfBoundsValue:-1000];
    return self;
}

- (void)getOrientation:(float[6])orientation
{
    double doubleOrientation[6];
    NSInteger i;
    
    [self getOrientationDouble:doubleOrientation];
    
    for (i = 0; i < 6; i++) {
        orientation[i] = doubleOrientation[i];
    }
}

- (void)getOrientationDouble:(double[6])orientation
{
    N3AffineTransform pixelToDicomTransform;
    N3Vector xBasis;
    N3Vector yBasis;
    
    pixelToDicomTransform = N3AffineTransformInvert(_volumeTransform);
    
    xBasis = N3VectorNormalize(N3VectorMake(pixelToDicomTransform.m11, pixelToDicomTransform.m12, pixelToDicomTransform.m13));
    yBasis = N3VectorNormalize(N3VectorMake(pixelToDicomTransform.m21, pixelToDicomTransform.m22, pixelToDicomTransform.m23));
    
    orientation[0] = xBasis.x; orientation[1] = xBasis.y; orientation[2] = xBasis.z;
    orientation[3] = yBasis.x; orientation[4] = yBasis.y; orientation[5] = yBasis.z; 
}

- (float)originX
{
    N3AffineTransform pixelToDicomTransform;
    
    pixelToDicomTransform = N3AffineTransformInvert(_volumeTransform);
    
    return pixelToDicomTransform.m41;
}

- (float)originY
{
    N3AffineTransform pixelToDicomTransform;
    
    pixelToDicomTransform = N3AffineTransformInvert(_volumeTransform);
    
    return pixelToDicomTransform.m42;
}

- (float)originZ
{
    N3AffineTransform pixelToDicomTransform;
    
    pixelToDicomTransform = N3AffineTransformInvert(_volumeTransform);
    
    return pixelToDicomTransform.m43;
}


@end














