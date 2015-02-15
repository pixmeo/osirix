//
//  OSIDerivedROI.m
//  OsiriX_Lion
//
//  Created by Joël Spaltenstein on 11/29/14.
//  Copyright (c) 2014 OsiriX Team. All rights reserved.
//

#import "OSIDerivedROI.h"
#import "Notifications.h"

@interface OSIDerivedROI ()

- (void)_volumetricROIDidUpdateNotification:(NSNotification *)notification;
- (void)_updateVolumeDataNotification:(NSNotification *)notification;

- (OSIROI *)_cachedROI;

@end

@implementation OSIDerivedROI

@synthesize floatVolumeData = _floatVolumeData;
@synthesize baseROI = _baseROI;

@synthesize strokeColor = _strokeColor;
@synthesize strokeThickness = _strokeThickness;
@synthesize fillColor = _fillColor;

- (instancetype)initWithBaseROI:(OSIROI *)baseROI floatVolumeData:(OSIFloatVolumeData *)floatVolumeData name:(NSString *)name
                  derivingBlock:(OSIROI *(^)(OSIROI * baseROI, OSIFloatVolumeData *floatVolumeData))derivingBlock
{
    if ( (self = [super init])) {
        _name = [name copy];
        _baseROI = [baseROI retain];
        _floatVolumeData = [floatVolumeData retain];
        _derivingBlock = [derivingBlock copy];
        [self derive];

        [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(_volumetricROIDidUpdateNotification:) name:OsirixVolumetricROIDidUpdateNotification object:_baseROI];
        [[NSNotificationCenter defaultCenter] addObserver:self selector:@selector(_updateVolumeDataNotification:) name:OsirixUpdateVolumeDataNotification object:nil];
    }
    return self;
}

- (void)dealloc
{
    [[NSNotificationCenter defaultCenter] removeObserver:self];

    [_name release];
    _name = nil;
    [_baseROI release];
    _baseROI = nil;
    [_cachedROI release];
    _cachedROI = nil;
    [_floatVolumeData release];
    _floatVolumeData = nil;
    [_derivingBlock release];
    _derivingBlock = nil;
    [_fillColor release];
    _fillColor = nil;
    [_strokeColor release];
    _strokeColor = nil;

    [super dealloc];
}

- (OSIROI *(^)(OSIROI * baseROI, OSIFloatVolumeData *floatVolumeData))derivingBlock
{
    return _derivingBlock;
}

- (void)derive
{
    [_cachedROI release];
    _cachedROI = nil;
    _cachedROI = _derivingBlock(self.baseROI, self.floatVolumeData);
    [_cachedROI retain];
    _cachedROI.fillColor = _fillColor;
    _cachedROI.strokeColor = _strokeColor;
    _cachedROI.strokeThickness = _strokeThickness;
}

- (void)setStrokeColor:(NSColor *)strokeColor
{
    if (_strokeColor != strokeColor) {
        [_strokeColor release];
        _strokeColor = [strokeColor retain];
        [[self _cachedROI] setStrokeColor:_strokeColor];
    }
}

- (void)setFillColor:(NSColor *)fillColor
{
    if (_fillColor != fillColor) {
        [_fillColor release];
        _fillColor = [fillColor retain];
        [[self _cachedROI] setFillColor:_fillColor];
    }
}

- (void)setStrokeThickness:(CGFloat)strokeThickness
{
    if (_strokeThickness != strokeThickness) {
        _strokeThickness = strokeThickness;
        [[self _cachedROI] setStrokeThickness:_strokeThickness];
    }
}

- (OSIROI *)_cachedROI
{
    if (_cachedROI == nil) {
        [self derive];
    }
    return _cachedROI;
}

- (void)_volumetricROIDidUpdateNotification:(NSNotification *)notification
{
    [_cachedROI release];
    _cachedROI = nil;
    [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self];
}

- (void)_updateVolumeDataNotification:(NSNotification *)notification
{
    [_cachedROI release];
    _cachedROI = nil;
    [[NSNotificationCenter defaultCenter] postNotificationName:OsirixVolumetricROIDidUpdateNotification object:self];
}

- (NSString *)name
{
    return _name;
}

- (NSArray *)metricNames;
{
    return [[self _cachedROI] metricNames];
}

- (NSString *)unitForMetric:(NSString *)metric
{
    return [[self _cachedROI] unitForMetric:metric];
}

- (id)valueForMetric:(NSString *)metric floatVolumeData:(OSIFloatVolumeData *)floatVolumeData;
{
    return [[self _cachedROI] valueForMetric:metric floatVolumeData:floatVolumeData];
}

- (CGFloat)intensityMeanWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;
{
    return [[self _cachedROI] intensityMeanWithFloatVolumeData:floatVolumeData];
}

- (CGFloat)intensityMaxWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;
{
    return [[self _cachedROI] intensityMaxWithFloatVolumeData:floatVolumeData];
}

- (CGFloat)intensityMinWithFloatVolumeData:(OSIFloatVolumeData *)floatVolumeData;
{
    return [[self _cachedROI] intensityMinWithFloatVolumeData:floatVolumeData];
}

- (CGFloat)volume;
{
    return [[self _cachedROI] volume];
}

- (OSIROIFloatPixelData *)ROIFloatPixelDataForFloatVolumeData:(OSIFloatVolumeData *)floatVolume
{
    return [[self _cachedROI] ROIFloatPixelDataForFloatVolumeData:floatVolume];
}

- (OSIROIMask *)ROIMaskForFloatVolumeData:(OSIFloatVolumeData *)floatVolume;
{
    return [[self _cachedROI] ROIMaskForFloatVolumeData:floatVolume];
}

- (OSIMaskROI *)maskROIRepresention
{
    return [[self _cachedROI] maskROIRepresention];
}

- (NSArray *)convexHull
{
    return [[self _cachedROI] convexHull];
}

- (N3BezierPath *)bezierPath
{
    return [[self _cachedROI] bezierPath];
}

- (N3Vector)centerOfMass
{
    return [[self _cachedROI] centerOfMass];
}

- (void)drawRect:(NSRect)rect inSlab:(OSISlab)slab inCGLContext:(CGLContextObj)glContext pixelFormat:(CGLPixelFormatObj)pixelFormat dicomToPixTransform:(N3AffineTransform)dicomToPixTransform
{
    [[self _cachedROI] drawRect:rect inSlab:slab inCGLContext:glContext pixelFormat:pixelFormat dicomToPixTransform:dicomToPixTransform];
}


@end
