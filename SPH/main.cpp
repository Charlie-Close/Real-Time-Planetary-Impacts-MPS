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
        int i = 0;
        auto start = std::chrono::high_resolution_clock::now();
        while (true) {
            headless.step();
            if (i == HEADLESS_ITTERATION_RATE_REFRESH) {
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end - start;
                float fps = (float)HEADLESS_ITTERATION_RATE_REFRESH * STEPS_PER_FRAME / elapsed.count();
                std::cout << "Steps Per Second: " << std::to_string(fps) << '\r' << std::flush;
                i = 0;
                start = std::chrono::high_resolution_clock::now();
            }
            i++;
        }
    }

    return 0;
}
