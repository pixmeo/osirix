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

#import <Cocoa/Cocoa.h>
#import "N3Geometry.h"

@class CPRUnsignedInt16ImageRep;

CF_EXTERN_C_BEGIN

typedef NS_ENUM(NSInteger, CPRInterpolationMode) {
    CPRInterpolationModeLinear,
    CPRInterpolationModeNearestNeighbor,
	
	CPRInterpolationModeNone = 0xFFFFFF,
};

typedef struct { // build one of these on the stack and then use -[CPRVolumeData aquireInlineBuffer:] to initialize it.
    const float *floatBytes;
    
    float outOfBoundsValue;
    
    NSUInteger pixelsWide;
    NSUInteger pixelsHigh;
    NSUInteger pixelsDeep;
    
    NSUInteger pixelsWideTimesPixelsHigh; // just in the interest of not calculating this a million times...
    
    N3AffineTransform volumeTransform;
} CPRVolumeDataInlineBuffer;

// Interface to the data
@interface CPRVolumeData : NSObject {
    NSData *_floatData;
    float _outOfBoundsValue;
    
    NSUInteger _pixelsWide;
    NSUInteger _pixelsHigh;
    NSUInteger _pixelsDeep;
    
    N3AffineTransform _volumeTransform; // volumeTransform is the transform from Dicom (patient) space to pixel data
}


- (instancetype)initWithFloatBytesNoCopy:(const float *)floatBytes pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
               volumeTransform:(N3AffineTransform)volumeTransform outOfBoundsValue:(float)outOfBoundsValue freeWhenDone:(BOOL)freeWhenDone; // volumeTransform is the transform from Dicom (patient) space to pixel data

- (instancetype)initWithData:(NSData *)data pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
   volumeTransform:(N3AffineTransform)volumeTransform outOfBoundsValue:(float)outOfBoundsValue; // volumeTransform is the transform from Dicom (patient) space to pixel data

- (instancetype)initWithVolumeData:(CPRVolumeData *)volumeData;

@property (readonly) NSUInteger pixelsWide;
@property (readonly) NSUInteger pixelsHigh;
@property (readonly) NSUInteger pixelsDeep;

@property (readonly, getter=isRectilinear) BOOL rectilinear;

@property (readonly) CGFloat minPixelSpacing; // the smallest pixel spacing in any direction;
@property (readonly) CGFloat pixelSpacingX;// mm/pixel
@property (readonly) CGFloat pixelSpacingY;
@property (readonly) CGFloat pixelSpacingZ;

@property (readonly) float outOfBoundsValue;

@property (readonly) N3AffineTransform volumeTransform; // volumeTransform is the transform from Dicom (patient) space to pixel data

// will copy fill length*sizeof(float) bytes, the coordinates better be within the volume!!!
// a run a is a series of pixels in the x direction
- (BOOL)getFloatRun:(float *)buffer atPixelCoordinateX:(NSUInteger)x y:(NSUInteger)y z:(NSUInteger)z length:(NSUInteger)length; 

- (CPRUnsignedInt16ImageRep *)unsignedInt16ImageRepForSliceAtIndex:(NSUInteger)z;
- (CPRVolumeData *)volumeDataForSliceAtIndex:(NSUInteger)z;

- (CPRVolumeData *)volumeDataByApplyingTransform:(N3AffineTransform)transform;

// the first version of this function figures out the dimensions needed to fit the whole volume. Note that with the first version of this function the passed in transform may not
// be equal to the volumeTransform of the returned volumeData because a minimum cube of data needed to fit the was calculated. Any shift in the data is guaranteed to be a multiple
// of the basis vectors of the transform though.
- (instancetype)volumeDataResampledWithVolumeTransform:(N3AffineTransform)transform interpolationMode:(CPRInterpolationMode)interpolationsMode;
- (instancetype)volumeDataResampledWithVolumeTransform:(N3AffineTransform)transform pixelsWide:(NSUInteger)pixelsWide pixelsHigh:(NSUInteger)pixelsHigh pixelsDeep:(NSUInteger)pixelsDeep
                                     interpolationMode:(CPRInterpolationMode)interpolationsMode;

- (BOOL)getFloat:(float *)floatPtr atPixelCoordinateX:(NSUInteger)x y:(NSUInteger)y z:(NSUInteger)z; // returns YES if the float was sucessfully gotten
- (BOOL)getLinearInterpolatedFloat:(float *)floatPtr atDicomVector:(N3Vector)vector; // these are slower, use the inline buffer if you care about speed
- (BOOL)getNearestNeighborInterpolatedFloat:(float *)floatPtr atDicomVector:(N3Vector)vector; // these are slower, use the inline buffer if you care about speed

- (BOOL)aquireInlineBuffer:(CPRVolumeDataInlineBuffer *)inlineBuffer; // always return YES

// not done yet, will crash if given vectors that are outside of the volume
- (NSUInteger)tempBufferSizeForNumVectors:(NSUInteger)numVectors;
- (void)linearInterpolateVolumeVectors:(N3VectorArray)volumeVectors outputValues:(float *)outputValues numVectors:(NSUInteger)numVectors tempBuffer:(void *)tempBuffer;
// end not done

- (BOOL)isEqual:(id)object;

- (BOOL)isDataValid __deprecated;
- (void)invalidateData __deprecated;
- (void)releaseInlineBuffer:(CPRVolumeDataInlineBuffer *)inlineBuffer __deprecated; // CPRVolumeTransform can no longer be invalided so the this is no longer needed

@end


@interface CPRVolumeData (DCMPixAndVolume) // make a nice clean interface between the rest of of OsiriX that deals with pixlist and all their complications, and fill out our convenient data structure.

- (id) initWithWithPixList:(NSArray *)pixList volume:(NSData *)volume;

- (void)getOrientation:(float[6])orientation;
- (void)getOrientationDouble:(double[6])orientation;

@property (readonly) float originX;
@property (readonly) float originY;
@property (readonly) float originZ;

@end

CF_INLINE const float* CPRVolumeDataFloatBytes(CPRVolumeDataInlineBuffer *inlineBuffer)
{
	return inlineBuffer->floatBytes;
}


CF_INLINE float CPRVolumeDataGetFloatAtPixelCoordinate(CPRVolumeDataInlineBuffer *inlineBuffer, NSInteger x, NSInteger y, NSInteger z)
{
    bool outside;
    
    if (inlineBuffer->floatBytes) {
        outside = false;
        
        outside |= x < 0;
        outside |= y < 0;
        outside |= z < 0;
        outside |= x >= inlineBuffer->pixelsWide;
        outside |= y >= inlineBuffer->pixelsHigh;
        outside |= z >= inlineBuffer->pixelsDeep;
        
        if (!outside) {
            return (inlineBuffer->floatBytes)[x + y*inlineBuffer->pixelsWide + z*inlineBuffer->pixelsWideTimesPixelsHigh];
        } else {
            return inlineBuffer->outOfBoundsValue;
        }
    } else {
        return 0;
    }
}

CF_INLINE float CPRVolumeDataLinearInterpolatedFloatAtVolumeCoordinate(CPRVolumeDataInlineBuffer *inlineBuffer, CGFloat x, CGFloat y, CGFloat z) // coordinate in the pixel space
{
    float returnValue;
    
    NSInteger floorX = (x);
    NSInteger ceilX = floorX+1.0;
    NSInteger floorY = (y);
    NSInteger ceilY = floorY+1.0;
    NSInteger floorZ = (z);
    NSInteger ceilZ = floorZ+1.0;
    
    bool outside = false;
    outside |= floorX < 0;
    outside |= floorY < 0;
    outside |= floorZ < 0;
    outside |= ceilX >= inlineBuffer->pixelsWide;
    outside |= ceilY >= inlineBuffer->pixelsHigh;
    outside |= ceilZ >= inlineBuffer->pixelsDeep;
    
    if (outside || !inlineBuffer->floatBytes) {
        returnValue = inlineBuffer->outOfBoundsValue;
    } else {
        float xd = x - floorX;
        float yd = y - floorY;
        float zd = z - floorZ;
//        
//        float xda = 1.0f - xd;
//        float yda = 1.0f - yd;
//        float zda = 1.0f - zd;
//        
//        float i1 = CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, floorX, floorY, floorZ)*zda + CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, floorX, floorY, ceilZ)*zd;
//        float i2 = CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, floorX, ceilY, floorZ)*zda + CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, floorX, ceilY, ceilZ)*zd;
//        
//        float w1 = i1*yda + i2*yd;
//        
//        float j1 = CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, ceilX, floorY, floorZ)*zda + CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, ceilX, floorY, ceilZ)*zd;
//        float j2 = CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, ceilX, ceilY, floorZ)*zda + CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, ceilX, ceilY, ceilZ)*zd;
//        
//        float w2 = j1*yda + j2*yd;
//        
//        returnValue = w1*xda + w2*xd;
        
        
#define trilinFuncMacro(v,x,y,z,a,b,c,d,e,f,g,h)         \
t00 =   a + (x)*(b-a);      \
t01 =   c + (x)*(d-c);      \
t10 =   e + (x)*(f-e);      \
t11 =   g + (x)*(h-g);      \
t0  = t00 + (y)*(t01-t00);  \
t1  = t10 + (y)*(t11-t10);  \
v   =  t0 + (z)*(t1-t0);

        float A, B, C, D, E, F, G, H;
        float t00, t01, t10, t11, t0, t1;
        int Binc, Cinc, Dinc, Einc, Finc, Ginc, Hinc;
        int xinc, yinc, zinc;
        
        xinc = 1;
        yinc = (int)inlineBuffer->pixelsWide;
        zinc = (int)inlineBuffer->pixelsWideTimesPixelsHigh;
        
        // Compute the increments to get to the other 7 voxel vertices from A
        Binc = xinc;
        Cinc = yinc;
        Dinc = xinc + yinc;
        Einc = zinc;
        Finc = zinc + xinc;
        Ginc = zinc + yinc;
        Hinc = zinc + xinc + yinc;
        
        // Set values for the first pass through the loop
        const float *dptr = inlineBuffer->floatBytes + floorZ * zinc + floorY * yinc + floorX;
        A = *(dptr);
        B = *(dptr + Binc);
        C = *(dptr + Cinc);
        D = *(dptr + Dinc);
        E = *(dptr + Einc);
        F = *(dptr + Finc);
        G = *(dptr + Ginc);
        H = *(dptr + Hinc);
        
        trilinFuncMacro( returnValue, xd, yd, zd, A, B, C, D, E, F, G, H );
    }
    
    return returnValue;
}

CF_INLINE float CPRVolumeDataNearestNeighborInterpolatedFloatAtVolumeCoordinate(CPRVolumeDataInlineBuffer *inlineBuffer, CGFloat x, CGFloat y, CGFloat z) // coordinate in the pixel space
{
#if CGFLOAT_IS_DOUBLE
    NSInteger roundX = round(x);
    NSInteger roundY = round(y);
    NSInteger roundZ = round(z);
#else
    NSInteger roundX = roundf(x);
    NSInteger roundY = roundf(y);
    NSInteger roundZ = roundf(z);
#endif

    return CPRVolumeDataGetFloatAtPixelCoordinate(inlineBuffer, roundX, roundY, roundZ);
}

CF_INLINE float CPRVolumeDataLinearInterpolatedFloatAtDicomVector(CPRVolumeDataInlineBuffer *inlineBuffer, N3Vector vector) // coordinate in mm dicom space
{
    vector = N3VectorApplyTransform(vector, inlineBuffer->volumeTransform);
    return CPRVolumeDataLinearInterpolatedFloatAtVolumeCoordinate(inlineBuffer, vector.x, vector.y, vector.z);
}

CF_INLINE float CPRVolumeDataNearestNeighborInterpolatedFloatAtDicomVector(CPRVolumeDataInlineBuffer *inlineBuffer, N3Vector vector) // coordinate in mm dicom space
{
    vector = N3VectorApplyTransform(vector, inlineBuffer->volumeTransform);
    return CPRVolumeDataNearestNeighborInterpolatedFloatAtVolumeCoordinate(inlineBuffer, vector.x, vector.y, vector.z);
}

CF_INLINE float CPRVolumeDataLinearInterpolatedFloatAtVolumeVector(CPRVolumeDataInlineBuffer *inlineBuffer, N3Vector vector)
{
    return CPRVolumeDataLinearInterpolatedFloatAtVolumeCoordinate(inlineBuffer, vector.x, vector.y, vector.z);
}

CF_INLINE float CPRVolumeDataNearestNeighborInterpolatedFloatAtVolumeVector(CPRVolumeDataInlineBuffer *inlineBuffer, N3Vector vector)
{
    return CPRVolumeDataNearestNeighborInterpolatedFloatAtVolumeCoordinate(inlineBuffer, vector.x, vector.y, vector.z);
}

CF_EXTERN_C_END

