//
//  Parameters.h
//  SPH
//
//  Created by Charlie Close on 06/02/2025.
//

#ifndef Parameters_h
#define Parameters_h

// Hydro
#define VISCOSITY_ALPHA 1.5f
#define DENSITY_ETA .005f // half a percent error
#define MAX_SMOOTHING_LENGTH 1.f
#define MAX_DENSITY_NR_ITTERATIONS 50
#define N_NEIGHBOURS_ESTIM 128 //Estimate number of neighbours. Used for caching on density pass. If too small, can't use cache.
#define THREADGROUP_CACHE_SIZE 512
#define RESOLUTION_ETA 1.2348
#define PARTICLE_IN_RANGE_THRESHOLD 10 // For choosing which method to use with NR

// EOS
#define ANEOS_TEXTURE_RESOLUTION 2048
#define ANEOS_MIN_RHO 1e-6
#define ANEOS_MAX_RHO 1
#define ANEOS_MIN_U 1e-10
#define ANEOS_MAX_U 1

// Time stepping
#define CFL .2f
#define DT0 32.f
#define R_MAX 10
#define DT_MIN1 (DT0 / (1 << R_MAX))
#define DT_MIN 1 / DT_MIN1
#define STEPS_PER_FRAME 1 // Can set to 0 for headless mode to split frame up into multiple steps
//#define WENDLAND_C2_KERNEL
#define CUBIC_SPLINE_KERNEL

// Gravity
#define GRAVITY_ETA .005f // half a percent error
#define P 3 // Power of multipole expansion (max is 4). Higher is faster at lower etas.
#define N_EXPANSION_TERMS ((P+1)*(P+2)*(P+3)/6)
#define MAX_TREE_DEPTH 100
#define GRAVITY_SMOOTHING_LENGTH .48f
#define PLUMBER_EQUIVALENT 3
#define GRAVITY_MAX_RECURSION 3
#define MAX_UNCHECKED_POINTERS 512
#define MAX_CHILDREN_IN_LEAF 8
#define EXTRA_MEMORY_MULTIPLIER 1.1 // When claiming memory for an octree, claim some extra so if our tree grows, we don't need to reclaim
#define RED_MEMEORY_MULTIPLIER 0.8 // Scale memory down if we can

// Cells and sorting
#define CELL_WIDTH 0.3f
#define CELL_POWER 8 // total number of cells is 2^(3 * CELL_POWER)
#define SORTING_BLOCK_SIZE 256
#define SORTING_MASK_LENGTH 8
#define SORTING_BUCKET_NUMBER (1 << SORTING_MASK_LENGTH)
#define SORTING_BIT_MASK (SORTING_BUCKET_NUMBER - 1)
#define SORTING_ITTERATIONS (32 / SORTING_MASK_LENGTH)
#define SHUFFLE_FRAMES 1

// File
//#define FILEPATH "demo_impact_n60.hdf5"
//#define FILEPATH "demo_impact_n60_0023.hdf5"
//#define FILEPATH "demo_impact_n65.hdf5"
//#define FILEPATH "demo_impact_n70.hdf5"
#define FILEPATH "demo_impact_n50.hdf5"
//#define FILEPATH "demo_impact_n40.hdf5"
// #define FILEPATH "saves/state_79.hdf5"

// Rendering
#define PARTICLE_SIZE 0.13
#define PARTICLE_SUBDIVITIONS 0
#define STARTING_POSITION 315, 355, 125
#define STARTING_PITCH -10
#define STARTING_YAW 90
#define LIGHT_DIRECTION -1.0, -1.0, 0.5
#define SHADOW_MAP_RESOLUTION 1024
#define DENSITY_GRADIENT_SETTLING_ITTERATIONS 5

// Headless mode
#define HEADLESS false
#define SNAPSHOT_RESOLUTION 2048
#define HEADLESS_ITTERATION_RATE_REFRESH 16
#define START_SNAPSHOT 0 // if running from a saved state
#define SNAPSHOT_DIR "snapshots"
#define SAVES_DIR "saves"

#ifdef __METAL_VERSION__
#define M_1_PI 1 / M_PI_F
#endif
#ifdef WENDLAND_C2_KERNEL
#define GAMMA 1.936492
#define KERNEL_CONSTANT 10.5 * M_1_PI
#endif
#ifdef CUBIC_SPLINE_KERNEL
#define GAMMA 1.825742
#define KERNEL_CONSTANT 16 * M_1_PI
#endif

#define BOX_CENTER 318
#define BOX_SIZE 318

#endif /* Parameters_h */
