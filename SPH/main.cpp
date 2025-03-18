//
//  main.cpp
//  SPH
//
//  Created by Charlie Close on 22/01/2025.
//

#define NS_PRIVATE_IMPLEMENTATION
#define MTL_PRIVATE_IMPLEMENTATION
#define MTK_PRIVATE_IMPLEMENTATION
#define CA_PRIVATE_IMPLEMENTATION

#include "AppDelegate.hpp"
#include <chrono>
#include <thread>
#include "Parameters.h"
#include "Headless.hpp"
#include <iostream>
#include <csignal>

int main( int argc, char* argv[] )
{
    if (!HEADLESS) {
        // Stops multiple windows (silly bug which apple need to fix).
        // https://developer.apple.com/forums/thread/765445
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        
        // The actual interesting code is in Render.cpp. This is just the code necessary
        // to run a metal application.
        NS::AutoreleasePool* pAutoreleasePool = NS::AutoreleasePool::alloc()->init();
        MyAppDelegate del;
        
        NS::Application* pSharedApplication = NS::Application::sharedApplication();
        pSharedApplication->setDelegate( &del );
        pSharedApplication->run();
        
        pAutoreleasePool->release();
    } else {
        Headless headless;
        headless.run();
        std::cout << "Destroying..." << std::endl;
        headless.~Headless();
    }

    return 0;
}
