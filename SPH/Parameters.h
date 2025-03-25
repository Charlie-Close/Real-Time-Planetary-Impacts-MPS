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
#define VISCOSITY_BETA 3.f
#define ALPHA_MIN .1f
#define ALPHA_MAX 2.f
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
#define CFL .3f
#define MIN_DT 0.01f
#define MAX_DT 100.f
#define STEPS_PER_FRAME 2 // Can set to 0 for headless mode to split frame up into multiple steps
//#define WENDLAND_C2_KERNEL
#define CUBIC_SPLINE_KERNEL
//#define QUARTIC_SPLINE_KERNEL

// Gravity
#define GRAVITY_ETA .005f // half a percent error
#define P 2 // Power of multipole expansion (max is 4). Higher is faster at lower etas.
#define N_EXPANSION_TERMS ((P+1)*(P+2)*(P+3)/6)
#define MAX_TREE_DEPTH 100
#define GRAVITY_SMOOTHING_LENGTH .48f
#define PLUMBER_EQUIVALENT 3
#define GRAVITY_MAX_RECURSION 3
#define MAX_UNCHECKED_POINTERS 32 // For FMM, max pointers at lowest level
#define MAX_CHILDREN_IN_LEAF 8
#define EXTRA_MEMORY_MULTIPLIER 1.05 // When claiming memory for an octree, claim some extra so if our tree grows, we don't need to reclaim
#define RED_MEMEORY_MULTIPLIER 0.8 // Scale memory down if we can

// Cells and sorting
#define CELL_WIDTH .3f
#define CELL_POWER 8 // total number of cells is 2^(3 * CELL_POWER)
#define SORTING_BLOCK_SIZE 256
#define SORTING_MASK_LENGTH 8
#define SORTING_BUCKET_NUMBER (1 << SORTING_MASK_LENGTH)
#define SORTING_BIT_MASK (SORTING_BUCKET_NUMBER - 1)
#define SORTING_ITTERATIONS (32 / SORTING_MASK_LENGTH)
#define SHUFFLE_FRAMES 1

// File
//#define FILEPATH "demo_impact_n60.hdf5"
//#define FILEPATH "saves_n60/state_8.hdf5"
//#define FILEPATH "demo_impact_n60_0023.hdf5"
//#define FILEPATH "demo_impact_n65.hdf5"
//#define FILEPATH "saves_n65/state_12.hdf5"
//#define FILEPATH "demo_impact_n70.hdf5"
//#define FILEPATH "saves_n70_visc_v2/state_42.hdf5"
//#define FILEPATH "demo_impact_n75.hdf5"
//#define FILEPATH "saves_n75/state_17.hdf5"
//#define FILEPATH "saves_n70_v2/state_3.hdf5"
#define FILEPATH "demo_impact_n50.hdf5"
//#define FILEPATH "demo_impact_n40.hdf5"
// #define FILEPATH "saves_n50/state_3.hdf5"
//#define FILEPATH "mercury_n50.hdf5"
//#define FILEPATH "mercury_n50/state_35.hdf5"
//#define FILEPATH "venus_n50.hdf5"
//#define FILEPATH "earth_n50.hdf5"
//#define FILEPATH "mars_n50.hdf5"
//#define FILEPATH "jupiter_n50.hdf5"
//#define FILEPATH "saturn_n50.hdf5"
//#define FILEPATH "uranus_n50.hdf5"
//#define FILEPATH "neptune_n50.hdf5"
//#define FILEPATH "neptune_n50/state_150.hdf5"
//#define FILEPATH "demo_target_n50.hdf5"

// Rendering
#define PARTICLE_SIZE 150
//#define PARTICLE_SIZE 500
#define PARTICLE_SUBDIVITIONS 0
#define STARTING_POSITION 315, 355, 100
#define STARTING_PITCH -10
#define STARTING_YAW 90
#define LIGHT_DIRECTION -1.0, -1.0, 0.5
#define SHADOW_MAP_RESOLUTION 4096
#define DENSITY_GRADIENT_SETTLING_ITTERATIONS 3
#define ACTIVATE_ALL_STEPS 256

// Headless mode
#define HEADLESS false
#define SNAPSHOT_RESOLUTION 4096
#define HEADLESS_ITTERATION_RATE_REFRESH 16
#define START_SNAPSHOT 42 // if running from a saved state
#define SNAPSHOT_DIR "snapshots_n70_visc_v2"
#define SAVES_DIR "saves_n70_visc_v2"
#define N_SNAPSHOTTERS 4
#define POSITIONS { (simd::float3){315, 355, 100}, (simd::float3){350, 355, 260}, (simd::float3){262, 362, 280}, (simd::float3){248, 306, 238} }
#define PITCHES { -10, -30, -34, 3 }
#define YAWS { 90, 120, 40, 53 }
#define SNAPSHOT_PERIOD_SECONDS 25

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
#ifdef QUARTIC_SPLINE_KERNEL
#define GAMMA 2.018932
#define KERNEL_CONSTANT (15625 / 512) * M_1_PI
#endif


#define BOX_CENTER 318
#define BOX_SIZE 318

#endif /* Parameters_h */
