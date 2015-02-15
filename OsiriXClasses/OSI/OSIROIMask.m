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

#import "OSIROIMask.h"
#import "OSIFloatVolumeData.h"
#import "OSIROIMaskRunStack.h"
#include <Accelerate/Accelerate.h>

OSIROIMaskRun OSIROIMaskRunMake(NSRange widthRange, NSUInteger heightIndex, NSUInteger depthIndex, float intensity)
{
    OSIROIMaskRun maskRun = {widthRange, heightIndex, depthIndex, intensity};
    return maskRun;
}

NSUInteger OSIROIMaskRunFirstWidthIndex(OSIROIMaskRun maskRun)
{
    return maskRun.widthRange.location;
}

NSUInteger OSIROIMaskRunLastWidthIndex(OSIROIMaskRun maskRun)
{
    return NSMaxRange(maskRun.widthRange) - 1;
}

const OSIROIMaskRun OSIROIMaskRunZero = {{0.0, 0.0}, 0, 0, 1.0};

@interface OSIMaskIndexPredicateStandIn : NSObject
{
    float intensity;
    float ROIMaskIntensity;
    NSUInteger ROIMaskIndexX;
    NSUInteger ROIMaskIndexY;
    NSUInteger ROIMaskIndexZ;
}
@property (nonatomic, readwrite, assign) float intensity;
@property (nonatomic, readwrite, assign) float ROIMaskIntensity;
@property (nonatomic, readwrite, assign) NSUInteger ROIMaskIndexX;
@property (nonatomic, readwrite, assign) NSUInteger ROIMaskIndexY;
@property (nonatomic, readwrite, assign) NSUInteger ROIMaskIndexZ;
@end
@implementation OSIMaskIndexPredicateStandIn
@synthesize intensity;
@synthesize ROIMaskIntensity;
@synthesize ROIMaskIndexX;
@synthesize ROIMaskIndexY;
@synthesize ROIMaskIndexZ;
@end


NSComparisonResult OSIROIMaskCompareRunValues(NSValue *maskRun1Value, NSValue *maskRun2Value, void *context)
{
    OSIROIMaskRun maskRun1 = [maskRun1Value OSIROIMaskRunValue];
    OSIROIMaskRun maskRun2 = [maskRun2Value OSIROIMaskRunValue];
    
    return OSIROIMaskCompareRun(maskRun1, maskRun2);
}


NSComparisonResult OSIROIMaskCompareRun(OSIROIMaskRun maskRun1, OSIROIMaskRun maskRun2)
{
    if (maskRun1.depthIndex < maskRun2.depthIndex) {
        return NSOrderedAscending;
    } else if (maskRun1.depthIndex > maskRun2.depthIndex) {
        return NSOrderedDescending;
    }
    
    if (maskRun1.heightIndex < maskRun2.heightIndex) {
        return NSOrderedAscending;
    } else if (maskRun1.heightIndex > maskRun2.heightIndex) {
        return NSOrderedDescending;
    }
    
    if (maskRun1.widthRange.location < maskRun2.widthRange.location) {
        return NSOrderedAscending;
    } else if (maskRun1.widthRange.location > maskRun2.widthRange.location) {
        return NSOrderedDescending;
    }
    
    return NSOrderedSame;
}

int OSIROIMaskQSortCompareRun(const void *voidMaskRun1, const void *voidMaskRun2)
{
    const OSIROIMaskRun* maskRun1 = voidMaskRun1;
    const OSIROIMaskRun* maskRun2 = voidMaskRun2;
    
    if (maskRun1->depthIndex < maskRun2->depthIndex) {
        return NSOrderedAscending;
    } else if (maskRun1->depthIndex > maskRun2->depthIndex) {
        return NSOrderedDescending;
    }
    
    if (maskRun1->heightIndex < maskRun2->heightIndex) {
        return NSOrderedAscending;
    } else if (maskRun1->heightIndex > maskRun2->heightIndex) {
        return NSOrderedDescending;
    }
    
    if (maskRun1->widthRange.location < maskRun2->widthRange.location) {
        return NSOrderedAscending;
    } else if (maskRun1->widthRange.location > maskRun2->widthRange.location) {
        return NSOrderedDescending;
    }
    
    return NSOrderedSame;

}

BOOL OSIROIMaskRunsOverlap(OSIROIMaskRun maskRun1, OSIROIMaskRun maskRun2)
{
    if (maskRun1.depthIndex == maskRun2.depthIndex && maskRun1.heightIndex == maskRun2.heightIndex) {
        return NSIntersectionRange(maskRun1.widthRange, maskRun2.widthRange).length != 0;
    }
    
    return NO;
}

BOOL OSIROIMaskRunsAbut(OSIROIMaskRun maskRun1, OSIROIMaskRun maskRun2)
{
    if (maskRun1.depthIndex == maskRun2.depthIndex && maskRun1.heightIndex == maskRun2.heightIndex) {
        if (NSMaxRange(maskRun1.widthRange) == maskRun2.widthRange.location ||
            NSMaxRange(maskRun2.widthRange) == maskRun1.widthRange.location) {
            return YES;
        }
    }
    return NO;
}

N3Vector OSIROIMaskIndexApplyTransform(OSIROIMaskIndex maskIndex, N3AffineTransform transform)
{
    return N3VectorApplyTransform(N3VectorMake(maskIndex.x, maskIndex.y, maskIndex.z), transform);
}

BOOL OSIROIMaskIndexInRun(OSIROIMaskIndex maskIndex, OSIROIMaskRun maskRun)
{
	if (maskIndex.y != maskRun.heightIndex || maskIndex.z != maskRun.depthIndex) {
		return NO;
	}
	if (NSLocationInRange(maskIndex.x, maskRun.widthRange)) {
		return YES;
	} else {
		return NO;
	}
}

NSArray *OSIROIMaskIndexesInRun(OSIROIMaskRun maskRun)
{
	NSMutableArray *indexes;
	NSUInteger i;
	OSIROIMaskIndex index;
	
	indexes = [NSMutableArray array];
	index.y = maskRun.heightIndex;
	index.z = maskRun.depthIndex;
	
	for (i = maskRun.widthRange.location; i < NSMaxRange(maskRun.widthRange); i++) {
		index.x = i;
		[indexes addObject:[NSValue valueWithOSIROIMaskIndex:index]];
	}
	return indexes;
}

@interface OSIROIMask ()
- (void)checkdebug;
@end

@implementation OSIROIMask

+ (instancetype)ROIMask
{
    return [[[[self class] alloc] init] autorelease];
}


+ (instancetype)ROIMaskWithSphereDiameter:(NSUInteger)diameter
{
    return [self ROIMaskWithElipsoidWidth:diameter height:diameter depth:diameter];
}

+ (instancetype)ROIMaskWithCubeSize:(NSUInteger)size
{
    return [[self class] ROIMaskWithBoxWidth:size height:size depth:size];
}

+ (instancetype)ROIMaskWithBoxWidth:(NSUInteger)width height:(NSUInteger)height depth:(NSUInteger)depth;
{
    NSUInteger i = 0;
    NSUInteger j = 0;
    
    OSIROIMaskRun *maskRuns = malloc(width * height * sizeof(OSIROIMaskRun));
    memset(maskRuns, 0, width * height * sizeof(OSIROIMaskRun));
    
    for (j = 0; j < height; j++) {
        for (i = 0; i < depth; i++) {
            maskRuns[(i*depth)+j] = OSIROIMaskRunMake(NSMakeRange(0, width), i, j, 1);
        }
    }

    return [[[OSIROIMask alloc] initWithSortedMaskRunData:[NSData dataWithBytesNoCopy:maskRuns length:width * height * sizeof(OSIROIMaskRun) freeWhenDone:YES]] autorelease];
}

+ (instancetype)ROIMaskWithElipsoidWidth:(NSUInteger)width height:(NSUInteger)height depth:(NSUInteger)depth
{
    NSUInteger i = 0;
    NSUInteger j = 0;
    NSUInteger k = 0;

    OSIROIMaskRun *maskRuns = malloc(height * depth * sizeof(OSIROIMaskRun));
    memset(maskRuns, 0, height * depth * sizeof(OSIROIMaskRun));

    CGFloat widthRadius = 0.5*(CGFloat)width;
    CGFloat heightRadius = 0.5*(CGFloat)height;
    CGFloat depthRadius = 0.5*(CGFloat)depth;
    NSInteger whiteSpace = 0;

    for (j = 0; j < depth; j++) {
        for (i = 0; i < height; i++) {
#if CGFLOAT_IS_DOUBLE
            CGFloat x = abs(((CGFloat)i)+.5-heightRadius);
            CGFloat y = abs(((CGFloat)j)+.5-depthRadius);
            if (radius*radius < x*x + y*y) {
                whiteSpace = -1;
            } else {
                whiteSpace = round(widthRadius - sqrt(widthRadius*widthRadius - x*x - y*y));
            }
#else
            CGFloat x = fabs(((CGFloat)i)+.5f-heightRadius);
            CGFloat y = fabs(((CGFloat)j)+.5f-depthRadius);
            if (widthRadius*widthRadius < x*x + y*y) {
                whiteSpace = -1;
            } else {
                whiteSpace = round(widthRadius - sqrt(widthRadius*widthRadius - x*x - y*y));
            }
#endif
            if (whiteSpace >= 0) {
                maskRuns[k] = OSIROIMaskRunMake(NSMakeRange(whiteSpace, width - (2 * whiteSpace)), i, j, 1);
                k++;
            }
        }
    }

    return [[[OSIROIMask alloc] initWithSortedMaskRunData:[NSData dataWithBytesNoCopy:maskRuns length:k * sizeof(OSIROIMaskRun) freeWhenDone:YES]] autorelease];
}

+ (instancetype)ROIMaskFromVolumeData:(OSIFloatVolumeData *)floatVolumeData __deprecated
{
    return [self ROIMaskFromVolumeData:floatVolumeData volumeTransform:NULL];
}


+ (id)ROIMaskFromVolumeData:(OSIFloatVolumeData *)floatVolumeData volumeTransform:(N3AffineTransformPointer)volumeTransformPtr
{
    NSInteger i;
    NSInteger j;
    NSInteger k;
    float intensity;
    NSMutableArray *maskRuns;
    OSIROIMaskRun maskRun;
    CPRVolumeDataInlineBuffer inlineBuffer;
        
    maskRuns = [NSMutableArray array];
    maskRun = OSIROIMaskRunZero;
    maskRun.intensity = 0.0;
    
    [floatVolumeData aquireInlineBuffer:&inlineBuffer];
    for (k = 0; k < inlineBuffer.pixelsDeep; k++) {
        for (j = 0; j < inlineBuffer.pixelsHigh; j++) {
            for (i = 0; i < inlineBuffer.pixelsWide; i++) {
                intensity = CPRVolumeDataGetFloatAtPixelCoordinate(&inlineBuffer, i, j, k);
                intensity = roundf(intensity*255.0f)/255.0f;
                
                if (intensity != maskRun.intensity) { // maybe start a run, maybe close a run
                    if (maskRun.intensity != 0) { // we need to end the previous run
                        [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
                        maskRun = OSIROIMaskRunZero;
                        maskRun.intensity = 0.0;
                    }
                    
                    if (intensity != 0) { // we need to start a new mask run
                        maskRun.depthIndex = k;
                        maskRun.heightIndex = j;
                        maskRun.widthRange = NSMakeRange(i, 1);
                        maskRun.intensity = intensity;
                    }
                } else  { // maybe extend a run // maybe do nothing
                    if (intensity != 0) { // we need to extend the run
                        maskRun.widthRange.length += 1;
                    }
                }
            }
            // after each run scan line we need to close out any open mask run
            if (maskRun.intensity != 0) {
                [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
                maskRun = OSIROIMaskRunZero;
                maskRun.intensity = 0.0;
            }
        }
    }

    if (volumeTransformPtr) {
        *volumeTransformPtr = floatVolumeData.volumeTransform;
    }

    return [[[[self class] alloc] initWithMaskRuns:maskRuns] autorelease];    
}

- (instancetype)init
{
	if ( (self = [super init]) ) {
		_maskRuns = [[NSArray alloc] init];
	}
	return self;
}

- (instancetype)initWithMaskRuns:(NSArray *)maskRuns
{
	if ( (self = [super init]) ) {
		_maskRuns = [[maskRuns sortedArrayUsingFunction:OSIROIMaskCompareRunValues context:NULL] retain];
        [self checkdebug];
	}
	return self;
}

- (instancetype)initWithMaskRunData:(NSData *)maskRunData
{
    NSMutableData *mutableMaskRunData = [maskRunData mutableCopy];
    
    qsort([mutableMaskRunData mutableBytes], [mutableMaskRunData length]/sizeof(OSIROIMaskRun), sizeof(OSIROIMaskRun), OSIROIMaskQSortCompareRun);
	
    id maskRun = [self initWithSortedMaskRunData:mutableMaskRunData];
    [mutableMaskRunData release];
    return maskRun;
}

- (instancetype)initWithSortedMaskRunData:(NSData *)maskRunData
{
	if ( (self = [super init]) ) {
		_maskRunsData = [maskRunData retain];
        [self checkdebug];
	}
	return self;
}

- (instancetype)initWithSortedMaskRuns:(NSArray *)maskRuns
{
	if ( (self = [super init]) ) {
		_maskRuns = [maskRuns retain];
        [self checkdebug];
	}
	return self;
}

- (instancetype)initWithIndexes:(NSArray *)maskIndexes
{
    NSMutableData *maskData = [NSMutableData dataWithLength:[maskIndexes count] * sizeof(OSIROIMaskIndex)];
    OSIROIMaskIndex *maskIndexArray = [maskData mutableBytes];
    NSUInteger i;
    for (i = 0; i < [maskIndexes count]; i++) {
        maskIndexArray[i] = [[maskIndexes objectAtIndex:i] OSIROIMaskIndexValue];
    }
    
    return [self initWithIndexData:maskData];
}

- (instancetype)initWithIndexData:(NSData *)indexData
{
    if ( (self = [super init]) ) {
        OSIROIMaskIndex *indexes = (OSIROIMaskIndex *)[indexData bytes];
        NSUInteger indexCount = [indexData length] / sizeof(OSIROIMaskIndex);
        NSUInteger i;
        NSMutableArray *maskRuns = [NSMutableArray array];
        OSIROIMaskRun maskRun = OSIROIMaskRunZero;
        
        if (indexCount == 0) {
            _maskRuns = [[NSArray alloc] init];
            return self;
        }
        
        for (i = 0; i < indexCount; i++) {
            maskRun.widthRange.location = indexes[i].x;
            maskRun.widthRange.length = 1;
            maskRun.heightIndex = indexes[i].y;
            maskRun.depthIndex = indexes[i].z;
            [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
        }
        
        
        NSArray *sortedMaskRuns = [maskRuns sortedArrayUsingFunction:OSIROIMaskCompareRunValues context:NULL];
        NSMutableArray *newSortedRuns = [NSMutableArray array];
        
        maskRun = [[sortedMaskRuns objectAtIndex:0] OSIROIMaskRunValue];
        
        for (i = 1; i < indexCount; i++) {
            OSIROIMaskRun sortedRun = [[sortedMaskRuns objectAtIndex:i] OSIROIMaskRunValue];
            
            if (NSMaxRange(maskRun.widthRange) == sortedRun.widthRange.location &&
                maskRun.heightIndex == sortedRun.heightIndex &&
                maskRun.depthIndex == sortedRun.depthIndex) {
                maskRun.widthRange.length++;
            } else if (OSIROIMaskRunsOverlap(maskRun, sortedRun)) {
                NSLog(@"overlap?");
            } else {
                [newSortedRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
                maskRun = sortedRun;
            }
        }
        
        [newSortedRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
		_maskRuns = [newSortedRuns retain];
        [self checkdebug];
	}
	return self;
}

- (instancetype)initWithSortedIndexes:(NSArray *)maskIndexes
{
    NSMutableData *maskData = [NSMutableData dataWithLength:[maskIndexes count] * sizeof(OSIROIMaskIndex)];
    OSIROIMaskIndex *maskIndexArray = [maskData mutableBytes];
    NSUInteger i;
    for (i = 0; i < [maskIndexes count]; i++) {
        maskIndexArray[i] = [[maskIndexes objectAtIndex:i] OSIROIMaskIndexValue];
    }
    
    return [self initWithSortedIndexData:maskData];
}

- (instancetype)initWithSortedIndexData:(NSData *)indexData
{
    if ( (self = [super init]) ) {
        OSIROIMaskIndex *indexes = (OSIROIMaskIndex *)[indexData bytes];
        NSUInteger indexCount = [indexData length];
        NSUInteger i;
        NSMutableArray *maskRuns = [NSMutableArray array];
        
        if (indexCount == 0) {
            _maskRuns = [[NSArray alloc] init];
            return self;
        }
        
        OSIROIMaskRun maskRun = OSIROIMaskRunZero;
        maskRun.widthRange.location = indexes[0].x;
        maskRun.widthRange.length = 1;
        maskRun.heightIndex = indexes[0].y;
        maskRun.depthIndex = indexes[0].z;
        
        for (i = 1; i < indexCount; i++) {
            if (maskRun.widthRange.location + 1 == indexes[i].x &&
                maskRun.heightIndex == indexes[1].y &&
                maskRun.depthIndex == indexes[1].z) {
                maskRun.widthRange.length++;
            } else {
                [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
                maskRun.widthRange.location = indexes[i].x;
                maskRun.widthRange.length = 1;
                maskRun.heightIndex = indexes[i].y;
                maskRun.depthIndex = indexes[i].z;
            }
        }
        
        [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRun]];
		_maskRuns = [maskRuns retain];
        [self checkdebug];
	}
	return self;
}

+ (BOOL)supportsSecureCoding
{
    return YES;
}

- (instancetype)initWithCoder:(NSCoder *)aDecoder
{
    NSData *maskRunsData = [aDecoder decodeObjectOfClass:[NSData class] forKey:@"maskRunsData"];
    return [self initWithSortedMaskRunData:maskRunsData];
}

- (void)encodeWithCoder:(NSCoder *)aCoder
{
    [aCoder encodeObject:[self maskRunsData] forKey:@"maskRunsData"];
}

- (instancetype)copyWithZone:(NSZone *)zone
{
    return [[[self class] allocWithZone:zone] initWithSortedMaskRunData:[[[self maskRunsData] copy] autorelease]];
}

- (void)dealloc
{
    [_maskRunsData release];
    _maskRunsData = nil;
    [_maskRuns release];
    _maskRuns = nil;
    
    [super dealloc];
}

- (OSIROIMask *)ROIMaskByTranslatingByX:(NSInteger)x Y:(NSInteger)y Z:(NSInteger)z
{
    const OSIROIMaskRun *maskRuns = (const OSIROIMaskRun *)[[self maskRunsData] bytes];
    NSInteger maskRunCount = [self maskRunCount];

    OSIROIMaskRun *newMaskRuns = malloc(maskRunCount * sizeof(OSIROIMaskRun));
    memset(newMaskRuns, 0, maskRunCount * sizeof(OSIROIMaskRun));
    NSUInteger newMaskRunsIndex = 0;
    NSUInteger i;

    for (i = 0; i < maskRunCount; i++) {
        if ((NSInteger)OSIROIMaskRunLastWidthIndex(maskRuns[i]) >= -x &&
            (NSInteger)maskRuns[i].heightIndex >= -y &&
            (NSInteger)maskRuns[i].depthIndex >= -z) {

            newMaskRuns[newMaskRunsIndex] = maskRuns[i];

            if ((NSInteger)OSIROIMaskRunFirstWidthIndex(newMaskRuns[newMaskRunsIndex]) < -x) {
                newMaskRuns[newMaskRunsIndex].widthRange.length += x + (NSInteger)OSIROIMaskRunFirstWidthIndex(newMaskRuns[newMaskRunsIndex]);
                newMaskRuns[newMaskRunsIndex].widthRange.location = 0;
            } else {
                newMaskRuns[newMaskRunsIndex].widthRange.location += x;
            }

            newMaskRuns[newMaskRunsIndex].heightIndex += y;
            newMaskRuns[newMaskRunsIndex].depthIndex += z;

            newMaskRunsIndex++;
        }
    }

    return [[[OSIROIMask alloc] initWithSortedMaskRunData:[NSData dataWithBytesNoCopy:newMaskRuns length:newMaskRunsIndex * sizeof(OSIROIMaskRun) freeWhenDone:YES]] autorelease];
}

- (OSIROIMask *)ROIMaskByIntersectingWithMask:(OSIROIMask *)otherMask
{
    return [self ROIMaskBySubtractingMask:[self ROIMaskBySubtractingMask:otherMask]];
}

- (OSIROIMask *)ROIMaskByUnioningWithMask:(OSIROIMask *)otherMask
{
    NSUInteger index1 = 0;
    NSUInteger index2 = 0;
    
//    OSIROIMaskRun run1;
//    OSIROIMaskRun run2;
    
    OSIROIMaskRun runToAdd = OSIROIMaskRunZero;
    OSIROIMaskRun accumulatedRun = OSIROIMaskRunZero;
    accumulatedRun.widthRange.length = 0;
    
    NSData *maskRun1Data = [self maskRunsData];
    NSData *maskRun2Data = [otherMask maskRunsData];
    const OSIROIMaskRun *maskRunArray1 = [maskRun1Data bytes];
    const OSIROIMaskRun *maskRunArray2 = [maskRun2Data bytes];
    
    NSMutableData *resultMaskRuns = [NSMutableData data];
    
    
    while (index1 < [self maskRunCount] || index2 < [otherMask maskRunCount]) {
        if (index1 < [self maskRunCount] && index2 < [otherMask maskRunCount]) {
            if (OSIROIMaskCompareRun(maskRunArray1[index1], maskRunArray2[index2]) == NSOrderedAscending) {
                runToAdd = maskRunArray1[index1];
                index1++;
            } else {
                runToAdd = maskRunArray2[index2];
                index2++;
            }
        } else if (index1 < [self maskRunCount]) {
            runToAdd = maskRunArray1[index1];
            index1++;
        } else {
            runToAdd = maskRunArray2[index2];
            index2++;
        }
        
        if (accumulatedRun.widthRange.length == 0) {
            accumulatedRun = runToAdd;
        } else if (OSIROIMaskRunsOverlap(runToAdd, accumulatedRun) || OSIROIMaskRunsAbut(runToAdd, accumulatedRun)) {
            if (NSMaxRange(runToAdd.widthRange) > NSMaxRange(accumulatedRun.widthRange)) {
                accumulatedRun.widthRange.length = NSMaxRange(runToAdd.widthRange) - accumulatedRun.widthRange.location;
            }
        } else {
            [resultMaskRuns appendBytes:&accumulatedRun length:sizeof(OSIROIMaskRun)];
            accumulatedRun = runToAdd;
        }
    }
    
    if (accumulatedRun.widthRange.length != 0) {
        [resultMaskRuns appendBytes:&accumulatedRun length:sizeof(OSIROIMaskRun)];
    }
    
    return [[[OSIROIMask alloc] initWithSortedMaskRunData:resultMaskRuns] autorelease];
}

- (OSIFloatVolumeData *)floatVolumeDataRepresentationWithVolumeTransform:(N3AffineTransform)volumeTransform;
{
    NSUInteger maxHeight = NSIntegerMin;
    NSUInteger minHeight = NSIntegerMax;
    NSUInteger maxDepth = NSIntegerMin;
    NSUInteger minDepth = NSIntegerMax;
    NSUInteger maxWidth = NSIntegerMin;
    NSUInteger minWidth = NSIntegerMax;

    [self extentMinWidth:&minWidth maxWidth:&maxWidth minHeight:&minHeight maxHeight:&maxHeight minDepth:&minDepth maxDepth:&maxDepth];

    NSUInteger width = (maxWidth - minWidth) + 1;
    NSUInteger height = (maxHeight - minHeight) + 1;
    NSUInteger depth = (maxDepth - minDepth) + 1;

    float *floatBytes = calloc(width * height * depth, sizeof(float));
    if (floatBytes == 0) {
        NSLog(@"%s wasn't able to allocate a buffer of size %ld", __PRETTY_FUNCTION__, width * height * depth * sizeof(float));
        return nil;
    }

    OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[self maskRunsData] bytes];
    NSInteger maskRunCount = [self maskRunCount];
    NSInteger i;

    // draw in the runs
    for (i = 0; i < maskRunCount; i++) {
        NSInteger x = maskRuns[i].widthRange.location - minWidth;
        NSInteger y = maskRuns[i].heightIndex - minHeight;
        NSInteger z = maskRuns[i].depthIndex - minDepth;

        vDSP_vfill(&(maskRuns[i].intensity), &(floatBytes[x + y*width + z*width*height]), 1, maskRuns[i].widthRange.length);
    }
    NSData *floatData = [NSData dataWithBytesNoCopy:floatBytes length:width * height * depth * sizeof(float)];

    // since we shifted the data, we need to shift the volumeTransform as well.
    N3AffineTransform shiftedVolumeTransform = N3AffineTransformConcat(volumeTransform, N3AffineTransformMakeTranslation(-1.0*(CGFloat)minWidth, -1.0*(CGFloat)minHeight, -1.0*(CGFloat)minDepth));

    return [[[OSIFloatVolumeData alloc] initWithData:floatData pixelsWide:width pixelsHigh:height pixelsDeep:depth volumeTransform:shiftedVolumeTransform outOfBoundsValue:0] autorelease];
}

- (OSIROIMask *)ROIMaskBySubtractingMask:(OSIROIMask *)subtractMask
{
    OSIROIMaskRunStack *templateRunStack = [[OSIROIMaskRunStack alloc] initWithMaskRunData:[self maskRunsData]];
    OSIROIMaskRun newMaskRun;
    NSUInteger length;

    NSUInteger subtractIndex = 0;
    NSData *subtractData = [subtractMask maskRunsData];
    NSInteger subtractDataCount = [subtractData length]/sizeof(OSIROIMaskRun);
    const OSIROIMaskRun *subtractRunArray = [subtractData bytes];
   
    NSMutableData *resultMaskRuns = [NSMutableData data];
    OSIROIMaskRun tempMaskRun;

    while (subtractIndex < subtractDataCount && [templateRunStack count]) {
        if (OSIROIMaskRunsOverlap([templateRunStack currentMaskRun], subtractRunArray[subtractIndex]) == NO) {
            if (OSIROIMaskCompareRun([templateRunStack currentMaskRun], subtractRunArray[subtractIndex]) == NSOrderedAscending) {
                tempMaskRun = [templateRunStack currentMaskRun];
                [resultMaskRuns appendBytes:&tempMaskRun length:sizeof(OSIROIMaskRun)];
                [templateRunStack popMaskRun];
            } else {
                subtractIndex++;
            }
        } else {
            // run the 4 cases
            if (NSLocationInRange([templateRunStack currentMaskRun].widthRange.location, subtractRunArray[subtractIndex].widthRange)) {
                if (NSLocationInRange(NSMaxRange([templateRunStack currentMaskRun].widthRange) - 1, subtractRunArray[subtractIndex].widthRange)) {
                    // 1.
                    [templateRunStack popMaskRun];
                } else {
                    // 2.
                    newMaskRun = [templateRunStack currentMaskRun];
                    length = NSIntersectionRange([templateRunStack currentMaskRun].widthRange, subtractRunArray[subtractIndex].widthRange).length;
                    newMaskRun.widthRange.location += length;
                    newMaskRun.widthRange.length -= length;
                    [templateRunStack popMaskRun];
                    [templateRunStack pushMaskRun:newMaskRun];
                    assert(newMaskRun.widthRange.length > 0);
                }
            } else {
                if (NSLocationInRange(NSMaxRange([templateRunStack currentMaskRun].widthRange) - 1, subtractRunArray[subtractIndex].widthRange)) {
                    // 4.
                    newMaskRun = [templateRunStack currentMaskRun];
                    length = NSIntersectionRange([templateRunStack currentMaskRun].widthRange, subtractRunArray[subtractIndex].widthRange).length;
                    newMaskRun.widthRange.length -= length;
                    [templateRunStack popMaskRun];
                    [templateRunStack pushMaskRun:newMaskRun];
                    assert(newMaskRun.widthRange.length > 0);
                } else {
                    // 3.
                    OSIROIMaskRun originalMaskRun = [templateRunStack currentMaskRun];
                    [templateRunStack popMaskRun];
                    
                    newMaskRun = originalMaskRun;
                    length = NSMaxRange(subtractRunArray[subtractIndex].widthRange) - originalMaskRun.widthRange.location;
                    newMaskRun.widthRange.location += length;
                    newMaskRun.widthRange.length -= length;
                    [templateRunStack pushMaskRun:newMaskRun];
                    assert(newMaskRun.widthRange.length > 0);

                    
                    newMaskRun = originalMaskRun;
                    length = NSMaxRange(originalMaskRun.widthRange) - subtractRunArray[subtractIndex].widthRange.location;
                    newMaskRun.widthRange.length -= length;
                    [templateRunStack pushMaskRun:newMaskRun];
                    assert(newMaskRun.widthRange.length > 0);
                }
            }
        }
    }
    
    while ([templateRunStack count]) {
        tempMaskRun = [templateRunStack currentMaskRun];
        [resultMaskRuns appendBytes:&tempMaskRun length:sizeof(OSIROIMaskRun)];
        [templateRunStack popMaskRun];
    }

    [templateRunStack release];
    return [[[OSIROIMask alloc] initWithSortedMaskRunData:resultMaskRuns] autorelease];
}

- (OSIROIMask *)ROIMaskCroppedToWidth:(NSUInteger)width height:(NSUInteger)height depth:(NSUInteger)depth
{
    const OSIROIMaskRun *maskRuns = (const OSIROIMaskRun *)[[self maskRunsData] bytes];
    NSInteger maskRunCount = [self maskRunCount];
    NSInteger i;
    NSUInteger badRuns = 0; // runs that are totally outside the bounds
    NSUInteger clippedRuns = 0; // runs that are partially outside the bounds and will need to be clipped

    for (i = 0; i < maskRunCount; i++) {
        if (OSIROIMaskRunFirstWidthIndex(maskRuns[i]) >= width || maskRuns[i].heightIndex >= height || maskRuns[i].depthIndex >= depth) {
            badRuns++;
        } else if (OSIROIMaskRunLastWidthIndex(maskRuns[i]) >= width) {
            clippedRuns++;
        }
    }

    if (badRuns + clippedRuns == 0) {
        return self;
    }

    NSUInteger newMaskRunsCount = maskRunCount - badRuns;

    if (newMaskRunsCount == 0) {
        return [OSIROIMask ROIMask];
    }

    OSIROIMaskRun *newMaskRuns = malloc(newMaskRunsCount * sizeof(OSIROIMaskRun));
    memset(newMaskRuns, 0, newMaskRunsCount * sizeof(OSIROIMaskRun));
    NSUInteger newMaskRunsIndex = 0;

    for (i = 0; i < maskRunCount; i++) {
        if (OSIROIMaskRunFirstWidthIndex(maskRuns[i]) < width &&
            maskRuns[i].heightIndex < height &&
            maskRuns[i].depthIndex < depth) {

            newMaskRuns[newMaskRunsIndex] = maskRuns[i];

            if (OSIROIMaskRunLastWidthIndex(maskRuns[i]) >= width) {
                newMaskRuns[newMaskRunsIndex].widthRange.length = (width - newMaskRuns[newMaskRunsIndex].widthRange.location);
            }
            newMaskRunsIndex++;
        }
    }

    return [[[OSIROIMask alloc] initWithSortedMaskRunData:[NSData dataWithBytesNoCopy:newMaskRuns length:newMaskRunsCount * sizeof(OSIROIMaskRun) freeWhenDone:YES]] autorelease];
}

- (BOOL)intersectsMask:(OSIROIMask *)otherMask // probably could use a faster implementation...
{
    OSIROIMask *intersection = [self ROIMaskByIntersectingWithMask:otherMask];
    return [intersection maskRunCount] > 0;
}

- (BOOL)isEqualToMask:(OSIROIMask *)otherMask // super lazy implementation FIXME!
{
    OSIROIMask *intersection = [self ROIMaskByIntersectingWithMask:otherMask];
    OSIROIMask *subMask1 = [self ROIMaskBySubtractingMask:intersection];
    OSIROIMask *subMask2 = [otherMask ROIMaskBySubtractingMask:intersection];
    
    return [subMask1 maskRunCount] == 0 && [subMask2 maskRunCount] == 0;
}


- (OSIROIMask *)filteredROIMaskUsingPredicate:(NSPredicate *)predicate floatVolumeData:(OSIFloatVolumeData *)floatVolumeData
{
    NSMutableArray *newMaskArray = [NSMutableArray array];
    OSIROIMaskRun activeMaskRun;
    BOOL isMaskRunActive = NO;
    float intensity;
    OSIMaskIndexPredicateStandIn *standIn = [[[OSIMaskIndexPredicateStandIn alloc] init] autorelease];

    for (NSValue *maskRunValue in [self maskRuns]) {
        OSIROIMaskRun maskRun = [maskRunValue OSIROIMaskRunValue];
        
        OSIROIMaskIndex maskIndex;
        maskIndex.y = maskRun.heightIndex;
        maskIndex.z = maskRun.depthIndex;

        standIn.ROIMaskIntensity = maskRun.intensity;
        standIn.ROIMaskIndexY = maskIndex.y;
        standIn.ROIMaskIndexZ = maskIndex.z;
        
        for (maskIndex.x = maskRun.widthRange.location; maskIndex.x < NSMaxRange(maskRun.widthRange); maskIndex.x++) {
            [floatVolumeData getFloat:&intensity atPixelCoordinateX:maskIndex.x y:maskIndex.y z:maskIndex.z];
            standIn.ROIMaskIndexX = maskIndex.x;
            standIn.intensity = intensity;

            if ([predicate evaluateWithObject:standIn]) {
                if (isMaskRunActive) {
                    activeMaskRun.widthRange.length++;
                } else {
                    activeMaskRun.widthRange.location = maskIndex.x;
                    activeMaskRun.widthRange.length = 1;
                    activeMaskRun.heightIndex = maskIndex.y;
                    activeMaskRun.depthIndex = maskIndex.z;
                    activeMaskRun.intensity = maskRun.intensity;
                    isMaskRunActive = YES;
                }
            } else {
                if (isMaskRunActive) {
                    [newMaskArray addObject:[NSValue valueWithOSIROIMaskRun:activeMaskRun]];
                    isMaskRunActive = NO;
                }
            }
        }
        if (isMaskRunActive) {
            [newMaskArray addObject:[NSValue valueWithOSIROIMaskRun:activeMaskRun]];
            isMaskRunActive = NO;
        }
    }
    
    OSIROIMask *filteredMask = [[[OSIROIMask alloc] initWithSortedMaskRuns:newMaskArray] autorelease];
    [filteredMask checkdebug];
    return filteredMask;
}

- (NSArray *)maskRuns 
{
    if (_maskRuns == nil) {
        NSUInteger maskRunCount = [_maskRunsData length]/sizeof(OSIROIMaskRun);
        const OSIROIMaskRun *maskRunArray = [_maskRunsData bytes];
        NSMutableArray *maskRuns = [[NSMutableArray alloc] initWithCapacity:maskRunCount];
        NSUInteger i;
        for (i = 0; i < maskRunCount; i++) {
            [maskRuns addObject:[NSValue valueWithOSIROIMaskRun:maskRunArray[i]]];
        }
        _maskRuns = maskRuns;
    }
    
	return _maskRuns;
}

- (NSData *)maskRunsData
{
    OSIROIMaskRun *maskRunArray;
    NSInteger i;
    
    if (_maskRunsData == nil) {
        maskRunArray = malloc([_maskRuns count] * sizeof(OSIROIMaskRun));
        
        for (i = 0; i < [_maskRuns count]; i++) {
            maskRunArray[i] = [[_maskRuns objectAtIndex:i] OSIROIMaskRunValue];
        }
        
        _maskRunsData = [[NSData alloc] initWithBytesNoCopy:maskRunArray length:[_maskRuns count] * sizeof(OSIROIMaskRun) freeWhenDone:YES];
    }
    
    return _maskRunsData;
}

- (NSUInteger)maskRunCount
{
    if (_maskRuns) {
        return [_maskRuns count];
    } else {
        return [_maskRunsData length] / sizeof(OSIROIMaskRun);
    }
}

- (NSUInteger)maskIndexCount
{
    NSData *maskRunData = [self maskRunsData];
    const OSIROIMaskRun *maskRunArray = [maskRunData bytes];
    NSUInteger maskRunCount = [self maskRunCount];
    NSUInteger maskIndexCount = 0;
    NSUInteger i = 0;
    
    for (i = 0; i < maskRunCount; i++) {
        maskIndexCount += maskRunArray[i].widthRange.length;
    }
    
    return maskIndexCount;
}

- (NSArray *)maskIndexes
{
	NSValue *maskRunValue;
	NSMutableArray *indexes;
    OSIROIMaskRun maskRun;
	
	indexes = [NSMutableArray array];
			   
	for (maskRunValue in [self maskRuns]) {
        maskRun = [maskRunValue OSIROIMaskRunValue];
        if (maskRun.intensity) {
            [indexes addObjectsFromArray:OSIROIMaskIndexesInRun(maskRun)];
        }
	}
			   
	return indexes;
}

- (BOOL)indexInMask:(OSIROIMaskIndex)index
{
    return [self containsIndex:index];
}

- (BOOL)containsIndex:(OSIROIMaskIndex)index;
{
    // since the runs are sorted, we can binary search
    NSUInteger runIndex = 0;
    NSUInteger runCount = 0;
    
    NSData *maskRunsData = [self maskRunsData];
    const OSIROIMaskRun *maskRuns = [maskRunsData bytes];
    runCount = [self maskRunCount];

    while (runCount) {
        NSUInteger middleIndex = runIndex + (runCount / 2);
        if (OSIROIMaskIndexInRun(index, maskRuns[middleIndex])) {
            return YES;
        }
        
        BOOL before = NO;
        if (index.z < maskRuns[middleIndex].depthIndex) {
            before = YES;
        } else if (index.z == maskRuns[middleIndex].depthIndex && index.y < maskRuns[middleIndex].heightIndex) {
            before = YES;
        } else if (index.z == maskRuns[middleIndex].depthIndex && index.y == maskRuns[middleIndex].heightIndex && index.x < maskRuns[middleIndex].widthRange.location) {
            before = YES;
        }
        
        if (before) {
            runCount /= 2; 
        } else {
            runIndex = middleIndex + 1;
            runCount = (runCount - 1) / 2;
        }
    }
    
    return NO;
}

- (instancetype)ROIMaskByResamplingFromVolumeTransform:(N3AffineTransform)fromTransform toVolumeTransform:(N3AffineTransform)toTransform interpolationMode:(CPRInterpolationMode)interpolationsMode
{
    if (N3AffineTransformEqualToTransform(fromTransform, toTransform)) {
        return self;
    }

    if ([self maskRunCount] == 0) {
        return self;
    }

    // The implementation of this function can be made a lot less memory demanding my only sampling one slice at a time instead of the whole volume
    OSIROIMask *resampledMask = nil;
    N3AffineTransform toVolumeTransform = N3AffineTransformIdentity;
    N3Vector shift = N3VectorZero;
    @autoreleasepool {
        OSIFloatVolumeData *fromVolumeData = [self floatVolumeDataRepresentationWithVolumeTransform:fromTransform];
        OSIFloatVolumeData *toVolumeData = [fromVolumeData volumeDataResampledWithVolumeTransform:toTransform interpolationMode:interpolationsMode];
        resampledMask = [OSIROIMask ROIMaskFromVolumeData:toVolumeData volumeTransform:&toVolumeTransform];

        // volumeDataResampledWithVolumeTransform can shift the data so that it doesn't store more than it needs to, so figure out how much the shift was, and translate the mask so that is is at the right place
        shift = N3VectorApplyTransform(N3VectorZero, N3AffineTransformConcat(N3AffineTransformInvert(toTransform), toVolumeTransform));

#if CGFLOAT_IS_DOUBLE
        shift.x = round(shift.x);
        shift.y = round(shift.y);
        shift.z = round(shift.z);
#else
        shift.x = roundf(shift.x);
        shift.y = roundf(shift.y);
        shift.z = roundf(shift.z);
#endif

        resampledMask = [[resampledMask ROIMaskByTranslatingByX:(NSInteger)-shift.x Y:(NSInteger)-shift.y Z:(NSInteger)-shift.z] retain];
    }

    return [resampledMask autorelease];
}


- (void)extentMinWidth:(NSUInteger*)minWidthPtr maxWidth:(NSUInteger*)maxWidthPtr minHeight:(NSUInteger*)minHeightPtr maxHeight:(NSUInteger*)maxHeightPtr minDepth:(NSUInteger*)minDepthPtr maxDepth:(NSUInteger*)maxDepthPtr;
{
    NSUInteger maxWidth = 0;
    NSUInteger minWidth = NSUIntegerMax;
    NSUInteger maxHeight = 0;
    NSUInteger minHeight = NSUIntegerMax;
    NSUInteger maxDepth = 0;
    NSUInteger minDepth = NSUIntegerMax;

    OSIROIMaskRun *maskRuns = (OSIROIMaskRun *)[[self maskRunsData] bytes];
    NSInteger maskRunCount = [self maskRunCount];
    NSInteger i;

    if (maskRunCount == 0) {
        if (minWidthPtr) {
            *minWidthPtr = 0;
        }
        if (maxWidthPtr) {
            *maxWidthPtr = 0;
        }
        if (minHeightPtr) {
            *minHeightPtr = 0;
        }
        if (maxHeightPtr) {
            *maxHeightPtr = 0;
        }
        if (minDepthPtr) {
            *minDepthPtr = 0;
        }
        if (maxDepthPtr) {
            *maxDepthPtr = 0;
        }

        return;
    }

    for (i = 0; i < maskRunCount; i++) {
        maxWidth = MAX(maxWidth, (NSInteger)OSIROIMaskRunLastWidthIndex(maskRuns[i]));
        minWidth = MIN(minWidth, (NSInteger)OSIROIMaskRunFirstWidthIndex(maskRuns[i]));

        maxHeight = MAX(maxHeight, (NSInteger)maskRuns[i].heightIndex);
        minHeight = MIN(minHeight, (NSInteger)maskRuns[i].heightIndex);

        maxDepth = MAX(maxDepth, (NSInteger)maskRuns[i].depthIndex);
        minDepth = MIN(minDepth, (NSInteger)maskRuns[i].depthIndex);
    }

    if (minWidthPtr) {
        *minWidthPtr = minWidth;
    }
    if (maxWidthPtr) {
        *maxWidthPtr = maxWidth;
    }
    if (minHeightPtr) {
        *minHeightPtr = minHeight;
    }
    if (maxHeightPtr) {
        *maxHeightPtr = maxHeight;
    }
    if (minDepthPtr) {
        *minDepthPtr = minDepth;
    }
    if (maxDepthPtr) {
        *maxDepthPtr = maxDepth;
    }
}

- (NSArray *)convexHull
{
    NSUInteger maxHeight = NSIntegerMin;
    NSUInteger minHeight = NSIntegerMax;
    NSUInteger maxDepth = NSIntegerMin;
    NSUInteger minDepth = NSIntegerMax;
    NSUInteger maxWidth = NSIntegerMin;
    NSUInteger minWidth = NSIntegerMax;

    [self extentMinWidth:&minWidth maxWidth:&maxWidth minHeight:&minHeight maxHeight:&maxHeight minDepth:&minDepth maxDepth:&maxDepth];

    NSMutableArray *hull = [NSMutableArray arrayWithCapacity:8];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(minWidth, minDepth, minHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(minWidth, maxDepth, minHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(maxWidth, maxDepth, minHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(maxWidth, minDepth, minHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(minWidth, minDepth, maxHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(minWidth, maxDepth, maxHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(maxWidth, maxDepth, maxHeight)]];
    [hull addObject:[NSValue valueWithN3Vector:N3VectorMake(maxWidth, minDepth, maxHeight)]];
    
    return hull;
}

- (N3Vector)centerOfMass
{
    NSData *maskData = [self maskRunsData];
    NSInteger runCount = [maskData length]/sizeof(OSIROIMaskRun);
    const OSIROIMaskRun *runArray = [maskData bytes];
    NSUInteger i;
    CGFloat floatCount = 0;
    N3Vector centerOfMass = N3VectorZero;
    
    for (i = 0; i < runCount; i++) {
        centerOfMass.x += ((CGFloat)runArray[i].widthRange.location+((CGFloat)runArray[i].widthRange.length/2.0)) * (CGFloat)runArray[i].widthRange.length;
        centerOfMass.y += (CGFloat)runArray[i].heightIndex*(CGFloat)runArray[i].widthRange.length;
        centerOfMass.z += (CGFloat)runArray[i].depthIndex*(CGFloat)runArray[i].widthRange.length;
        floatCount += runArray[i].widthRange.length;
    }
    
    centerOfMass.x /= floatCount;
    centerOfMass.y /= floatCount;
    centerOfMass.z /= floatCount;
    
    return centerOfMass;
}

- (NSString *)description
{
    NSMutableString *desc = [NSMutableString string];
	[desc appendString:NSStringFromClass([self class])];
    [desc appendFormat:@"\nMask Run Count: %lld\n", (long long)[self maskRunCount]];
    [desc appendFormat:@"Index Count: %lld\n", (long long)[self maskIndexCount]];

    if ([self maskRunCount] > 0) {
        NSUInteger minWidth = 0;
        NSUInteger minHeight = 0;
        NSUInteger minDepth = 0;
        NSUInteger maxWidth = 0;
        NSUInteger maxHeight = 0;
        NSUInteger maxDepth = 0;
        [self extentMinWidth:&minWidth maxWidth:&maxWidth minHeight:&minHeight maxHeight:&maxHeight minDepth:&minDepth maxDepth:&maxDepth];
        [desc appendFormat:@"Width  Range: %4ld...%-4ld\n", (long )minWidth, (long)maxWidth];
        [desc appendFormat:@"Height Range: %4ld...%-4ld\n", (long )minHeight, (long)maxHeight];
        [desc appendFormat:@"Depth  Range: %4ld...%-4ld\n", (long )minDepth, (long)maxDepth];
    }

    [desc appendString:@"{\n"];

    NSUInteger maskRunsCount = [self maskRunCount];
    const OSIROIMaskRun *maskRuns = [[self maskRunsData] bytes];
    NSUInteger i;

    for (i = 0; i < maskRunsCount; i++) {
        [desc appendFormat:@"X:%4ld...%-4ld Y:%-4ld Z:%-4ld\n", (long)OSIROIMaskRunFirstWidthIndex(maskRuns[i]), (long)OSIROIMaskRunLastWidthIndex(maskRuns[i]), (long)maskRuns[i].heightIndex, (long)maskRuns[i].depthIndex];
    }

    [desc appendString:@"}"];
    return desc;
}

- (void)checkdebug
{
#ifndef NDEBUG
    // make sure that all the runs are in order.
    assert(_maskRuns || _maskRunsData);
    NSInteger i;
    if (_maskRunsData) {
        NSInteger maskRunsDataCount = [_maskRunsData length]/sizeof(OSIROIMaskRun);
        const OSIROIMaskRun *maskRunArray = [_maskRunsData bytes];
        for (i = 0; i < (maskRunsDataCount - 1); i++) {
            assert(OSIROIMaskCompareRun(maskRunArray[i], maskRunArray[i+1]) == NSOrderedAscending);
            assert(OSIROIMaskRunsOverlap(maskRunArray[i], maskRunArray[i+1]) == NO);
        }
        for (i = 0; i < maskRunsDataCount; i++) {
            assert(maskRunArray[i].widthRange.length > 0);
        }
    }
    
    if (_maskRuns) {
        for (i = 0; i < ((NSInteger)[_maskRuns count]) - 1; i++) {
            assert(OSIROIMaskCompareRunValues([_maskRuns objectAtIndex:i], [_maskRuns objectAtIndex:i+1], NULL) == NSOrderedAscending);
            assert(OSIROIMaskRunsOverlap([[_maskRuns objectAtIndex:i] OSIROIMaskRunValue], [[_maskRuns objectAtIndex:i+1] OSIROIMaskRunValue]) == NO);
        }
        for (i = 0; i < [_maskRuns count]; i++) {
            assert([[_maskRuns objectAtIndex:i] OSIROIMaskRunValue].widthRange.length > 0);
        }
    }
#endif
}





@end

@implementation NSValue (OSIMaskRun)

+ (NSValue *)valueWithOSIROIMaskRun:(OSIROIMaskRun)volumeRun
{
	return [NSValue valueWithBytes:&volumeRun objCType:@encode(OSIROIMaskRun)];
}

- (OSIROIMaskRun)OSIROIMaskRunValue
{
	OSIROIMaskRun run;
    assert(strcmp([self objCType], @encode(OSIROIMaskRun)) == 0);
    [self getValue:&run];
    return run;
}	

+ (NSValue *)valueWithOSIROIMaskIndex:(OSIROIMaskIndex)maskIndex
{
	return [NSValue valueWithBytes:&maskIndex objCType:@encode(OSIROIMaskIndex)];
}

- (OSIROIMaskIndex)OSIROIMaskIndexValue
{
	OSIROIMaskIndex index;
    assert(strcmp([self objCType], @encode(OSIROIMaskIndex)) == 0);
    [self getValue:&index];
    return index;
}

@end










