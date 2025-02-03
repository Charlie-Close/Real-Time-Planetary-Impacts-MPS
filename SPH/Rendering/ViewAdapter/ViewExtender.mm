//
//  ViewExtender.m
//  Lattice Boltzman
//
//  Created by Charlie Close on 08/12/2024.
//
//  I don't like this, but for some reason metalcpp does not have an interface for getting
//  user inputs. Therefore we have to code one ourselves using objective c++.
//

#import "ViewPass.m"
#import "ViewAdapter.hpp"

Camera* Explorer::ViewAdapter::_camera = nullptr;

Explorer::ViewAdapter::ViewAdapter(Camera* camera) {
    _camera = camera;
}

ViewExtender *adapter;

MTK::View *Explorer::ViewAdapter::get(CGRect frame) {
  [ViewExtender load: frame];
  return (__bridge MTK::View *)[ViewExtender get];
}

void Explorer::ViewAdapter::printDebug() {
  ViewExtender *ref = [ViewExtender get];
  [ref printDebug];
}



@implementation ViewExtender

+ (void)load:(CGRect)frame {
  NSLog(@"Loading Objective-c ViewAdapter ...");
  NSAutoreleasePool *pool = [[NSAutoreleasePool alloc] init];
  adapter = [[self alloc] initWithFrame:frame];
  [adapter init];
  [pool release];
}

+ (ViewExtender *)get {
  return adapter;
}

- (id)init {
  BOOL isFirstResponder = [self becomeFirstResponder];
  NSLog(@"Is first responder: %@", isFirstResponder ? @"Yes" : @"No");
  return self;
}

// Needs to be overwritten to accept keyboard events.
// Learn more here:
// https://developer.apple.com/library/archive/documentation/Cocoa/Conceptual/EventOverview/EventHandlingBasics/EventHandlingBasics.html#//apple_ref/doc/uid/10000060i-CH5
- (BOOL)acceptsFirstResponder {
  return YES;
}

- (void)keyDown:(NSEvent *)event {
  Explorer::ViewAdapter::_camera->handleKeyDown(event.keyCode);
  [event type];
}

- (void)keyUp:(NSEvent *)event {
  Explorer::ViewAdapter::_camera->handleKeyUp(event.keyCode);
  [event type];
}

- (void)mouseDown:(NSEvent *)event {
    Explorer::ViewAdapter::_camera->handleMouseDown(event.locationInWindow.x, event.locationInWindow.y);
}

- (void)mouseDragged:(NSEvent *)event {
    Explorer::ViewAdapter::_camera->handleMouseDrag(event.locationInWindow.x, event.locationInWindow.y);
}

- (void)printDebug {
  NSLog(@"ViewAdapter debug info ...");
}


@end

