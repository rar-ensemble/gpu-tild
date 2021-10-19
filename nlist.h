////////////////////////////////////
// Rob Riggleman        7/23/2021 //
// Calculates a neighbor list for //
// a specified group of particles.//
////////////////////////////////////

#include "group.h"

#ifndef _NLIST
#define _NLIST

// NOTE: the below rcut assumes whatever "skin" is needed is included
// in rcut; no additional skin is added, so rcut for NList should
// be larger than whatever functions use these lists.

class NList {
private:
    float rcut;         // Cut-off distance for neighbor list
    int GRID;           // Size of thread grid for cell operations
    int BLOCK;          // Num of blocks per thread
    string group_name;  // Name of the group in the n list
    int group_index;    // Index of the group

public:

    int total_ncells;   // total number of cells
    int CELL_MAX;       // max number of particles in each cell
    int* ncells;        // [Dim] number of cells in each direction
    int* d_ncells;      // [Dim] number of cells in each direction
    int* cell_count;    // [total_ncells] number of particles in each cell
    int* d_cell_count;  // [total_ncells] number of particles in each cell
    int* cell_inds;     // [total_ncells*CELL_MAX] 
    int* d_cell_inds;   //   id of particles in each cell
    int* particle_cell; // [ns] stores cell particle belongs to
    int* d_particle_cell; // [ns] stores cell particle belongs to
    int* neigh_count;   // [ns] number of neighbors per particles
    int* d_neigh_count; // [ns] number of neighbors per particles
    int* neigh_inds;    // [ns*CELL_MAX] indices of neighbors of a particle
    int* d_neigh_inds;  // [ns*CELL_MAX] indices of neighbors of a particle

    
    float* Lcell;   // [Dim] Cell size in each direction
    float* d_Lcell; // [Dim] Cell size in each direction

    void Initialize(float, string);
    void AssignCells(void);
    void BuildNList(void);
    NList();
    ~NList();
};

#endif