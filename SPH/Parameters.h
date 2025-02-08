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
#define DENSITY_ETA .001f
#define MAX_SMOOTHING_LENGTH .5f
#define MAX_DENSITY_NR_ITTERATIONS 50

// EOS
#define ANEOS_TEXTURE_RESOLUTION 2048
#define ANEOS_MIN_RHO 1e-6
#define ANEOS_MAX_RHO 2.98e-2
#define ANEOS_MIN_U 1e-10
#define ANEOS_MAX_U 8.8e-2

// Time stepping
#define CFL .5f
#define DT0 32.f
#define R_MAX 8
#define DT_MIN1 (DT0 / (1 << R_MAX))
#define DT_MIN 1 / DT_MIN1
#define STEPS_PER_FRAME 2

// Gravity
#define GRAVITY_ETA .005f
#define P 2
#define N_EXPANSION_TERMS ((P+1)*(P+2)*(P+3)/6)
#define MAX_TREE_DEPTH 6
#define GRAVITY_SMOOTHING_LENGTH 0.1
#define GRAVITY_MAX_RECURSION 1
#define MAX_UNCHECKED_POINTERS 1024
#define MAX_CHILDREN_IN_LEAF 4
#define EXTRA_MEMORY_MULTIPLIER 1.5 // When claiming memory for an octree, claim some extra so if our tree goes, we don't need to reclaim

// Cells and sorting
#define CELL_WIDTH 0.2f
#define CELL_POWER 8 // total number of cells is 2^(3 * CELL_POWER)
#define SORTING_BLOCK_SIZE 512
#define SORTING_MASK_LENGTH 8
#define SORTING_BUCKET_NUMBER (1 << SORTING_MASK_LENGTH)
#define SORTING_BIT_MASK (SORTING_BUCKET_NUMBER - 1)
#define SORTING_ITTERATIONS (32 / SORTING_MASK_LENGTH)

// File
//#define FILEPATH "demo_impact_n60.hdf5"
#define FILEPATH "demo_impact_n50.hdf5"
//#define FILEPATH "demo_impact_n40.hdf5"

// Rendering
#define PARTICLE_SIZE 0.12
#define PARTICLE_SUBDIVITIONS 1
#define STARTING_POSITION 305, 315, 175
#define STARTING_PITCH 0
#define STARTING_YAW 80
#define DENSITY_CUTOFF 0.005


#endif /* Parameters_h */
