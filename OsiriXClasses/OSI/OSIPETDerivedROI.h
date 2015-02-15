//
//  OSIPETDerivedROI.h
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 12/1/14.
//  Copyright (c) 2014 OsiriX Team. All rights reserved.
//

#import "OSIDerivedROI.h"

@interface OSIPETDerivedROI : OSIDerivedROI

- (instancetype)initWithBaseROI:(OSIROI *)roi floatVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name;

@end
