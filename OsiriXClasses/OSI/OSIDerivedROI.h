//
//  OSIDerivedROI.h
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 11/29/14.
//  Copyright (c) 2014 OsiriX Team. All rights reserved.
//

#import "OSIROI.h"


@interface OSIDerivedROI : OSIROI
{
    NSString *_name;
    OSIROI *_baseROI;
    OSIROI *_cachedROI;

    NSColor *_strokeColor;
    NSColor *_fillColor;
    CGFloat _strokeThickness;

    OSIFloatVolumeData *_floatVolumeData;
    OSIROI *(^_derivingBlock)(OSIROI * baseROI, OSIFloatVolumeData *floatVolumeData);
}

- (instancetype)initWithBaseROI:(OSIROI *)roi floatVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name
                  derivingBlock:(OSIROI *(^)(OSIROI * baseROI, OSIFloatVolumeData *floatVolumeData))derivingBlock;

- (void)derive;

@property (nonatomic, readonly, retain) OSIFloatVolumeData *floatVolumeData;
@property (nonatomic, readonly, retain) OSIROI *baseROI;

- (OSIROI *(^)(OSIROI * baseROI, OSIFloatVolumeData *floatVolumeData)) derivingBlock;

@property (nonatomic, readwrite, retain) NSColor *strokeColor;
@property (nonatomic, readwrite, assign) CGFloat strokeThickness;
@property (nonatomic, readwrite, retain) NSColor *fillColor;

@end
