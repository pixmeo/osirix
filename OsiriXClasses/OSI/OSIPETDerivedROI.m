//
//  OSIPETDerivedROI.m
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 12/1/14.
//  Copyright (c) 2014 OsiriX Team. All rights reserved.
//

#import "OSIPETDerivedROI.h"
#import "OSIROIFloatPixelData.h"
#import "OSIROIMask.h"
#import "OSIMaskROI.h"

@implementation OSIPETDerivedROI

- (instancetype)initWithBaseROI:(OSIROI *)roi floatVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name
{
    self = [super initWithBaseROI:roi floatVolumeData:floatVolumeData name:name derivingBlock:^OSIROI *(OSIROI *blockBaseROI, OSIFloatVolumeData *blockFloatVolumeData) {
        // the alorgithm is:
        // find the highest intesity pixel
        // make a 3 CM sphere around the pixel with that intensity

        OSIROIMask *roiMask = [blockBaseROI ROIMaskForFloatVolumeData:blockFloatVolumeData];
        if ([roiMask maskRunCount] == 0) { // if the mask is empty, return an empty mask
            return [[[OSIMaskROI alloc] initWithROIMask:[OSIROIMask ROIMask] volumeTransform:[floatVolumeData volumeTransform] name:@"derivedNULL"] autorelease];
        }

        OSIROIFloatPixelData *floatPixelData = [[[OSIROIFloatPixelData alloc] initWithROIMask:roiMask floatVolumeData:blockFloatVolumeData] autorelease];
        float intesityMax = [floatPixelData intensityMax];
        NSPredicate *predicate = [NSPredicate predicateWithFormat:@"self.intensity == %f", intesityMax];
        OSIROIMask *intesityMaxMask = [roiMask filteredROIMaskUsingPredicate:predicate floatVolumeData:blockFloatVolumeData];
        NSLog(@"%@", intesityMaxMask);
        OSIROIMaskIndex maskIndex = [[[intesityMaxMask maskIndexes] objectAtIndex:0] OSIROIMaskIndexValue];

        OSIMaskROI *derivedROI = [[[OSIMaskROI alloc] initWithROIMask:[OSIROIMask ROIMask] volumeTransform:[floatVolumeData volumeTransform] name:@"derived"] autorelease];
        OSIROIMask *sphereMask = [derivedROI ROIMaskWithSphereDiameter:30 atDicomVector:OSIROIMaskIndexApplyTransform(maskIndex, N3AffineTransformInvert([floatVolumeData volumeTransform]))];
        [derivedROI unionWithMask:sphereMask];

        return derivedROI;
    }];

    return self;
}

@end
