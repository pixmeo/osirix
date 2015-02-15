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

#import "OSIROI.h"
#import "OSIROI+Private.h"
#import "OSIPlanarPathROI.h"
#import "OSIPlanarBrushROI.h"
#import "OSIMaskROI.h"
#import "OSICoalescedPlanarROI.h"
#import "OSIROIFloatPixelData.h"
#import "OSIFloatVolumeData.h"
#import "DCMView.h"
#import "N3Geometry.h"
#import "ROI.h"

NSString* const OSIPasteboardTypeMaskROI = @"com.rossetantoine.osirix.maskROI";
NSString* const OSIPasteboardTypeCodingROI = @"com.rossetantoine.osirix.codingROI";

@implementation OSIROI

- (void)dealloc
{
    [super dealloc];
}

- (NSString *)name
{
	assert(0);
	return nil;
}

- (NSColor *)fillColor
{
    NSSet *osiriXROIs = [self osiriXROIs];
    NSColor *color = nil;
    
    for (ROI *roi in osiriXROIs) {
        if (color == nil) {
            if ([roi type] == tPlain) {
                color = [[roi NSColor] colorWithAlphaComponent:[roi opacity]];
            }
        } else if ([color isEqual:[roi NSColor]] == NO) {
            return nil;
        }
    }
    
    return color;
}


- (void)setFillColor:(NSColor *)color
{
    NSSet *osiriXROIs = [self osiriXROIs];
    
    for (ROI *roi in osiriXROIs) {
        if ([roi type] == tPlain) {
            [roi setNSColor:[color colorWithAlphaComponent:1]];
            [roi setOpacity:[color alphaComponent]];
        }
    }
}

- (NSColor *)strokeColor
{
    NSSet *osiriXROIs = [self osiriXROIs];
    NSColor *color = nil;
    
    for (ROI *roi in osiriXROIs) {
        if (color == nil) {
            if ([roi type] != tPlain) {
                color = [[roi NSColor] colorWithAlphaComponent:[roi opacity]];
            }
        } else if ([color isEqual:[roi NSColor]] == NO) {
            return nil;
        }
    }
    
    return color;
}

- (void)setStrokeColor:(NSColor *)color
{
    NSSet *osiriXROIs = [self osiriXROIs];
    
    for (ROI *roi in osiriXROIs) {
        if ([roi type] != tPlain) {
            [roi setNSColor:[color colorWithAlphaComponent:1]];
            [roi setOpacity:[color alphaComponent]];
        }
    }
}

- (CGFloat)strokeThickness
{
    NSSet *osiriXROIs = [self osiriXROIs];
    CGFloat thickness = 0;
    
    for (ROI *roi in osiriXROIs) {
        if (thickness == 0) {
            thickness = [roi thickness];
        } else if (thickness != [roi thickness]) {
            return 0;
        }
    }
    
    return thickness;
}

- (void)setStrokeThickness:(CGFloat)strokeThickness
{
    NSSet *osiriXROIs = [self osiriXROIs];
    
    for (ROI *roi in osiriXROIs) {
        [roi setThickness:strokeThickness];
    }
}

- (NSArray *)convexHull
{
	assert(0);
	return nil;
}

- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume
{
	return nil;
}

- (NSSet *)osiriXROIs
{
	return [NSSet set];
}

- (NSArray *)metricNames
{
	return [NSArray arrayWithObjects:@"intensityMean", @"intensityMax", @"intensityMin", @"volume", nil];
}

- (NSString *)localizedNameForMetric:(NSString *)metric
{
	if ([metric isEqualToString:@"intensityMean"]) {
		return @"Mean Intensity"; // localize me!
	} else if ([metric isEqualToString:@"intensityMax"]) {
		return @"Maximum Intensity"; // localize me!
	} else if ([metric isEqualToString:@"intensityMin"]) {
		return @"Minimum Intensity"; // localize me!
	}
	return nil;
}

- (NSString *)unitForMetric:(NSString *)metric
{
	if ([metric isEqualToString:@"intensityMean"]) {
		return @"";
	} else if ([metric isEqualToString:@"intensityMax"]) {
		return @""; 
	} else if ([metric isEqualToString:@"intensityMin"]) {
		return @"";
    } else if ([metric isEqualToString:@"volume"]) {
        return @"mm3";
    }

	return nil;
}

- (id)valueForMetric:(NSString *)metric floatVolumeData:(OSIFloatVolumeData *)floatVolumeData;
{
	if ([metric isEqualToString:@"intensityMean"]) {
		return [NSNumber numberWithDouble:[self intensityMeanWithFloatVolumeData:floatVolumeData]];
	} else if ([metric isEqualToString:@"intensityMax"]) {
		return [NSNumber numberWithDouble:[self intensityMaxWithFloatVolumeData:floatVolumeData]];
	} else if ([metric isEqualToString:@"intensityMin"]) {
		return [NSNumber numberWithDouble:[self intensityMinWithFloatVolumeData:floatVolumeData]];
	} else if ([metric isEqualToString:@"volume"]) {
		return [NSNumber numberWithDouble:(double)[self volume]];
	}
	return nil;
}

- (CGFloat)intensityMeanWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData
{
    return [[self ROIFloatPixelDataForFloatVolumeData:floatVolumeData] intensityMean];
}

- (CGFloat)intensityMaxWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData
{
    return [[self ROIFloatPixelDataForFloatVolumeData:floatVolumeData] intensityMax];
}

- (CGFloat)intensityMinWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData
{
    return [[self ROIFloatPixelDataForFloatVolumeData:floatVolumeData] intensityMin];
}

- (CGFloat)volume
{
    return -1;
}

- (OSIROIFloatPixelData *)ROIFloatPixelDataForFloatVolumeData:(OSIFloatVolumeData *)floatVolume; // convenience method
{
    OSIROIMask *roiMask;
    roiMask = [self ROIMaskForFloatVolumeData:floatVolume];

    assert([floatVolume checkDebugROIMask:roiMask]);
    
	return [[[OSIROIFloatPixelData alloc] initWithROIMask:roiMask floatVolumeData:floatVolume] autorelease];
}

- (OSIMaskROI *)maskROIRepresention
{
    return nil;
}

- (N3BezierPath *)bezierPath
{
    return nil;
}

- (N3Vector)centerOfMass
{
    assert(0);
    return N3VectorZero;
}

- (void)drawRect:(NSRect)rect inSlab:(OSISlab)slab inCGLContext:(CGLContextObj)glContext pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform
{
    [self drawSlab:slab inCGLContext:glContext pixelFormat:pixelFormat dicomToPixTransform:dicomToPixTransform];
}

- (void)drawSlab:(OSISlab)slab inCGLContext:(CGLContextObj)glContext pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform;
{
    
}

- (id)initWithPasteboardPropertyList:(id)propertyList ofType:(NSString *)type
{
    [self autorelease];
    self = nil;
    if ([type isEqualToString:OSIPasteboardTypeMaskROI] || [type isEqualToString:OSIPasteboardTypeCodingROI]) {
        self = [[NSKeyedUnarchiver unarchiveObjectWithData:propertyList] retain];
    }
    return self;
}

+ (NSArray *)readableTypesForPasteboard:(NSPasteboard *)pasteboard
{
    return @[OSIPasteboardTypeCodingROI, OSIPasteboardTypeMaskROI];
}

+ (NSPasteboardReadingOptions)readingOptionsForType:(NSString *)type pasteboard:(NSPasteboard *)pasteboard
{
    return NSPasteboardReadingAsData;
}

- (NSArray *)writableTypesForPasteboard:(NSPasteboard *)pasteboard
{
    NSMutableArray *types = [NSMutableArray array];
    if ([self conformsToProtocol:@protocol(NSCoding)]) {
        [types addObject:OSIPasteboardTypeCodingROI];
    }
    if ([self maskROIRepresention] != nil) {
        [types addObject:OSIPasteboardTypeMaskROI];
    }
    return types;
}

- (id)pasteboardPropertyListForType:(NSString *)type
{
    if ([type isEqualToString:OSIPasteboardTypeMaskROI] && [self isKindOfClass:[OSIMaskROI class]] == NO) {
        OSIMaskROI *maskROI = [self maskROIRepresention];
        return [NSKeyedArchiver archivedDataWithRootObject:maskROI];
    } if ([self conformsToProtocol:@protocol(NSCoding)]) {
        return [NSKeyedArchiver archivedDataWithRootObject:(id<NSCoding>)self];
    } else {
        return nil;
    }
}


@end

@implementation OSIROI (Private)

+ (id)ROIWithOsiriXROI:(ROI *)roi pixToDICOMTransfrom:(N3AffineTransform)pixToDICOMTransfrom;
{
	switch ([roi type]) {
		case tMesure:
		case tOPolygon:
		case tCPolygon:
		case tOval:
		case tROI:
        case tPencil:
			return [[[OSIPlanarPathROI alloc] initWithOsiriXROI:roi pixToDICOMTransfrom:pixToDICOMTransfrom] autorelease];
			break;
        case tPlain:
            return [[[OSIPlanarBrushROI alloc] initWithOsiriXROI:roi pixToDICOMTransfrom:pixToDICOMTransfrom] autorelease];
		default:
			return nil;;
	}
}


+ (id)ROICoalescedWithSourceROIs:(NSArray *)rois
{
	return [[[OSICoalescedPlanarROI alloc] initWithSourceROIs:rois] autorelease];
}


@end
