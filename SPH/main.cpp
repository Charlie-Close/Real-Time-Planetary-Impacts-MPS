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

bool run = true;

void handleSignal(int /*signal*/) {
    run = false;
}

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
        std::signal(SIGINT, handleSignal);
        int i = 0;
        auto start = std::chrono::high_resolution_clock::now();
        auto const constStart = start;
        float simTimStart = headless.getTime();
        const float constSimStart = simTimStart;
        int framesPerCout = HEADLESS_ITTERATION_RATE_REFRESH * fmax(1, STEPS_PER_FRAME);
        while (run) {
            headless.step();
            if (i == HEADLESS_ITTERATION_RATE_REFRESH) {
                auto end = std::chrono::high_resolution_clock::now();
                float simTimeEnd = headless.getTime();
                std::chrono::duration<double> elapsed = end - start;
                std::chrono::duration<double> totalElapsed = end - constStart;
                float simElapsed = simTimeEnd - simTimStart;
                float fps = (float)framesPerCout / elapsed.count();
                float secsPerSecs = simElapsed / elapsed.count();
                float avgSecsPerSecs = (simTimeEnd - constSimStart) / totalElapsed.count();
                float dtAvg = simElapsed / framesPerCout;
                std::cout << "Steps Per Second: " << std::to_string(fps);
                std::cout << "\nSimulations seconds per second: " << std::to_string(secsPerSecs);
                std::cout << "\nAverage simulations seconds per second: " << std::to_string(avgSecsPerSecs);
                std::cout << "\nAverage dt: " << std::to_string(dtAvg);
                std::cout << "\nSimulation time: " << std::to_string(simTimeEnd) << "\n" << std::endl;
                i = 0;
                start = std::chrono::high_resolution_clock::now();
                simTimStart = simTimeEnd;
            }
            i++;
        }
        std::cout << "Destroying..." << std::endl;
        headless.~Headless();
    }

    return 0;
}
