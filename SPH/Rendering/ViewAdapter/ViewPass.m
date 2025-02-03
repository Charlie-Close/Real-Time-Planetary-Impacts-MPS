//
//  ViewPass.m
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//

/**
 * 1) ViewAdapter -> ViewExtender (C++ -> Obj-C)
 * 2) ViewAdapter -> returns extension of MTK::View. (Obj-C -> C++)
 **/

#import <MetalKit/MetalKit.h>

//Class Camera; // Forward declare the C++ Camera class

@interface ViewExtender : MTKView {
    float _camera;
}
+ (void)load:(CGRect)frame;// withCamera:(void *)camera;
+ (ViewExtender *)get;
- (void)printDebug;
@end
