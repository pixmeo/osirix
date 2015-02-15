//
//  OSIMaskROI.h
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 9/26/12.
//  Copyright (c) 2012 OsiriX Team. All rights reserved.
//

#import "OSIROI.h"
#import "OSIROIMask.h"
#import "N3Geometry.h"

@interface OSIMaskROI : OSIROI <NSCoding>
{
    OSIROIMask *_mask;
    N3AffineTransform _volumeTransform; // volumeTransform is the transform from Dicom (patient) space to the mask
    NSString *_name;
    NSColor *_fillColor;
    
    OSIFloatVolumeData *_cachedBitmapMask;
}

- (instancetype)initWithROI:(OSIROI *)roi sampledOnVolumeData:(OSIFloatVolumeData *)floatVolumeData; // initializing by copying the values and getting the mask from the passed in ROI.

- (instancetype)initWithROIMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform name:(NSString *)name;

// reinterpolateMask will take the mask with it's volume transform, and reinterpolate so that the mask is aligned on the floatVolumeData
- (instancetype)initWithROIMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform sampledOnVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name reinterpolateMask:(BOOL)reinterpolateMask;


- (void)unionWithMask:(OSIROIMask *)mask;
- (void)unionWithMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform; // volumeTransform of the mask we are unioning with
- (void)subtractMask:(OSIROIMask *)mask;
- (void)subtractMask:(OSIROIMask *)mask volumeTransform:(N3AffineTransform)volumeTransform; // volumeTransform of the mask we are subtracting

- (OSIFloatVolumeData *)bitmapMask;

@property (nonatomic, readwrite, retain) OSIROIMask *mask;
@property (nonatomic, readwrite, assign) N3AffineTransform volumeTransform;


- (void)setVolumeTransform:(N3AffineTransform)volumeTransform reinterpolateMask:(BOOL)reinterpolateMask;

/** Creates and returns an OSIROIMask that is an apporimation of a sphere with the given diameter in DICOM coordinates at the position in DICOM coordinate.

 @return A new OSIROIMask object that contains the value of volumeRun.
 @param diameter The diameter of the sphere in mm.
 @param atDicomVector The center of the sphere in DICOM coordinates.
 */
- (OSIROIMask *)ROIMaskWithSphereDiameter:(CGFloat)diameter atDicomVector:(N3Vector)atDicomVector;

/** Creates and returns an OSIROIMask formed by dragging a sphere from fromDicomVector to toDicomVector.

 @return A new OSIROIMask object that contains the value of volumeRun.
 @param diameter The diameter of the sphere in mm.
 @param fromDicomVector The start of the brush stroke in DICOM coordinates.
 @param toDicomVector The end of the brush stroke in DICOM coordinates.
 */
- (OSIROIMask *)ROIMaskBrushWithSphereDiameter:(CGFloat)diameter fromDicomVector:(N3Vector)fromDicomVector toDicomVector:(N3Vector)toDicomVector;

@end
