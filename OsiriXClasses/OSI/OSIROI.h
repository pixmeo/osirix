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
#import "OSIGeometry.h"
#import <OpenGL/CGLTypes.h>

@class OSIStudy;
@class ROI;
@class OSIROIMask;
@class OSIFloatVolumeData;
@class OSIROIFloatPixelData;
@class N3BezierPath;
@class OSIMaskROI;

/**
 
 OSIROI is an abstract ROI class. It is the super class of all OsiriX SDK plugin ROI types and provides basic support for accessing pixel data and properties common to all ROIs.
 
 Subclasses must implement convexHull, name, and ROIMaskForFloatVolumeData:.
 
 @warning *Important:* For now the Plugin SDK classes only works with intensity data and does not work with RGB data.

 
 */

extern NSString* const OSIPasteboardTypeMaskROI;
extern NSString* const OSIPasteboardTypeCodingROI; // an ROI that implements NSCoding

@interface OSIROI : NSObject <NSPasteboardWriting, NSPasteboardReading>  {
}

//@property (nonatomic, readwrite, assign) void *context;

//- (id)initWithDictionaryRepresentation:(NSDictionary *)dict;

///-----------------------------------
/// @name Getting ROI Attributes
///-----------------------------------

/** Returns the name of the receiver.
 
 @return The name of the receiver
 */
- (NSString *)name;

/** Returns the fill color of the receiver.
 
 This value is equal to nil if the ROI should not be drawn, or if this OSI ROI is backed by multiple osirix ROIs with different colors.
 
 @return The fill color of the receiver
 */
- (NSColor *)fillColor;

/** Set the fill color of the receiver.
 
 Set this value to nil if the ROI should not be drawn.
 
 */
- (void)setFillColor:(NSColor *)color;

/** Returns the stroke color of the receiver.
 
 This value is equal to nil if the outline of the ROI should not be drawn, or if this OSI ROI is backed by multiple osirix ROIs with different colors.
 
 @return The stroke color of the receiver
 */
- (NSColor *)strokeColor;

/** Set the stroke color of the receiver.
 
 Set this value to nil if the ROI should not draw it's outline drawn.
 
 */
- (void)setStrokeColor:(NSColor *)strokeColor;

/** Returns the stroke thickness of the receiver.
 
 This value is 0 if this OSI ROI is backed by multiple osirix ROIs with stroke thicknesses.
 
 @return The stroke color of the receiver
 */
- (CGFloat)strokeThickness;

/** Set the stroke thickness of the receiver.
  
 */
- (void)setStrokeThickness:(CGFloat)strokeThickness;

/** Returns an array of `NSString` objects that represent the names of all the metrics that can be recovered from this ROI.
  
 Concrete subclasses may override this method to return additional metrics.
 
 @return An array of `NSString` objects that represent the names of all the metrics that can be recovered from this ROI.
 @see label
 @see labelForMetric:
 @see unitForMetric:
 @see valueForMetric:floatVolumeData:
*/
- (NSArray *)metricNames;

/** Returns a reasonable label to print for the given metric.
 
 Concrete subclasses may override this method to return a reasonable label for specific metrics.
 
 @return A reasonable label to print for the given metric.
 @param metric The metric name for which to return a label.
 @see label
 @see metricNames
 @see unitForMetric:
 @see valueForMetric:floatVolumeData:
 */
- (NSString *)localizedNameForMetric:(NSString *)metric;

/** Returns the unit for a given metric.
 
 Concrete subclasses may override this method to return units for any additional metric they may define.
 
 @return The unit for a given metric.
 @param metric The metric name for which to return the unit.
 @see label
 @see metricNames
 @see labelForMetric:
 @see valueForMetric:floatVolumeData:
 */
- (NSString *)unitForMetric:(NSString *)metric;

/** Returns the value for a given metric.
 
 Concrete subclasses may override this method to return values for any additional metric they may define.
  
 @return The value for a given metric.
 @param metric The metric name for which to return a value.
 @param floatVolumeData The OSIFloatVolumeData to use to generate the values.
 @see label
 @see metricNames
 @see labelForMetric:
 @see unitForMetric:
 */
- (id)valueForMetric:(NSString *)metric floatVolumeData:(OSIFloatVolumeData *)floatVolumeData;

/** Returns the mean intesity of this ROI under the given float volume data.
 
 @return The mean intesity of this ROI under the given float volume data.
 @param metric The metric name for which to return a value.
 @see label
 @see metricNames
 @see labelForMetric:
 @see unitForMetric:
 */
- (CGFloat)intensityMeanWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;

/** Returns the maximum intesity of this ROI under the given float volume data.
 
 @return The maximum intesity of this ROI under the given float volume data.
 @param metric The metric name for which to return a value.
 @see label
 @see metricNames
 @see labelForMetric:
 @see unitForMetric:
 */
- (CGFloat)intensityMaxWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;

/** Returns the minimum intesity of this ROI under the given float volume data.
 
 @return The minimum intesity of this ROI under the given float volume data.
 @param metric The metric name for which to return a value.
 @see label
 @see metricNames
 @see labelForMetric:
 @see unitForMetric:
 */
- (CGFloat)intensityMinWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;

/** Returns the volume of the ROI in cumbic centimeters
 
 @return The maximum intesity of this ROI under the given float volume data.
 @param metric The metric name for which to return a value.
 @see label
 @see metricNames
 @see labelForMetric:
 @see unitForMetric:
 */
- (CGFloat)volume;

///-----------------------------------
/// @name Getting Pixel/Voxel Data values
///-----------------------------------


//- (OSIStudy *)study;

/** Returns the a OSIROIFloatPixelData object that can be used to access pixels under this ROI in the given Float Volume data.

 This is a convenience for getting the ROI Mask for the given Float Volume Data and generating a Float Pixel Data object from that mask and the given Float Volume Data.

 @return The a OSIROIFloatPixelData object that can be used to access pixels under this ROI in the given Float Volume data.
 @param floatVolume the Float Volume Data for which to generate a Float Pixel Data object.
 @see ROIFloatPixelData
 @see ROIMaskForFloatVolumeData:
 */
- (OSIROIFloatPixelData *)ROIFloatPixelDataForFloatVolumeData:(OSIFloatVolumeData *)floatVolume; // convenience method

/** Returns the a OSIROIMask object that represents the volume the receiver covers in the given Float Volume data.
 
 Concrete subclasses must override this method.
 
 @return The a OSIROIFloatPixelData object that can be used to access pixels under this ROI in the given Float Volume data.
 @param floatVolume the Float Volume Data for which to generate a mask.
 @see ROIFloatPixelData
 @see ROIFloatPixelDataForFloatVolumeData:
 */
- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume;

/** Returns returns a OSIMaskROI that covers the pixels within the ROI, or nil if this ROI can't be represented as a mask.

 Concrete subclasses can implement this function if it is possible to represent this ROI as a mask.

 @return The a OSIROIMask object that includes the pixels under the mask.
 @see ROIFloatPixelDataForFloatVolumeData:
 */
- (OSIMaskROI *)maskROIRepresention;

//- (BOOL)containsVector:(OSIVector)vector;

//- (NSDictionary *)dictionaryRepresentation; // make sure this is a plist serializable dictionary;

///-----------------------------------
/// @name Representing the general position of the ROI
///-----------------------------------


/** Returns an array of points that represent the outside bounds of the ROI.
 
 Note: For now these points don't actually need to lie on the convex hull, but the ROI *must* be within the convex hull of these points.
 
 These points are in patient space.
 
 Concrete subclasses must implement this method.
 
 @return An array of N3Vectors stored at NSValues that represent the outside bounds of the ROI.
 */
- (NSArray *)convexHull; // N3Vectors stored in NSValue objects. The ROI promises to live inside of these points

/** If it makes sense for this ROI, returns a bezierPath that represents this ROI.
  
 Concrete subclasses can implement this method if it makes sense.
 
 @return A N3BezierPath representation of this ROI in DICOM space. Or nil if it does not make sense to represent this ROI as a bezier path.
 */
- (N3BezierPath *)bezierPath;

/** Returns the center of mass of the ROI.
 
 @return The center of mass of the ROI.
 */
- (N3Vector)centerOfMass;

///-----------------------------------
/// @name Drawing
///-----------------------------------

/** Overridden by subclasses to draw the receiver’s image within the passed-in rectangle..
 
 The receiver is expected to draw into OpenGL. The current OpenGL model matrix is set up so that rendering is in pix space.
  
 @param dirtyRect A rectangle defining the dirty area of the view that requires redrawing.
 @param dicomToPixTransform A matrix that converts points in Patient Space (Dicom space in mm) into pix space.
 
 @return An array of points that represent the outside bounds of the ROI.
 */
- (void)drawRect:(NSRect)rect inSlab:(OSISlab)slab inCGLContext:(CGLContextObj)glContext pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform;

- (void)drawSlab:(OSISlab)slab inCGLContext:(CGLContextObj)glContext pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform __deprecated;


// for drawing in 3D what we really want is for the ROI to return a VTK actor, and then it will be the actor and VTK that will decide how to draw

///-----------------------------------
/// @name Breaking out of the SDK
///-----------------------------------

/** Returns an array of OsiriX `ROI` objects that are the basis of this OSIROI.
 
 Concrete subclasses need to implement this method if the receiver depends on OsiriX `ROI` objects.
 
 @return A set of OsiriX `ROI` objects that are the basis of this OSIROI.
 */
- (NSSet *)osiriXROIs; // returns the primitive ROIs that are represented by this object

// at some point I would love to support drawing new ROI types...

@end
