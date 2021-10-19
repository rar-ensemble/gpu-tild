#include "globals.h"
#include "nlist.h"

__global__ void d_zero_int_vector(int*, int);
__global__ void d_assign_cells(float*, int*, int*, float*, int*, int, 
    int*, int*, int, int );
__device__ int d_stack_input(const int*, const int*, const int);
__device__ void d_unstack(int, int*, const int*, const int);
__global__ void d_make_nlist(const float*, const float, const int*,
    const int*, const int*, const int, const int*, int*, int*,
    const int*, const float*, const float*, const int, const int);

void NList::Initialize(float rc, string grp_name) {

    this->group_index = -1;
    for (int i = 0; i < n_groups; i++) {
        if (Groups.at(i).name == grp_name)
            this->group_index = i;
    }

    if (this->group_index == -1) {
        char death[120];
        sprintf(death, "Group %s not found to make a neighbor list!", grp_name.c_str());
        die(death);
    }

    this->group_name = grp_name;

    cout << "Making neighbor list for group " << grp_name << ". ";
    this->rcut = rc;
    this->Lcell = (float*)calloc(Dim, sizeof(float));
    this->ncells = (int*)calloc(Dim, sizeof(int));

    this->total_ncells = 1;
    float cell_vol = 1.0f;
    for (int j = 0; j < Dim; j++) {
        this->ncells[j] = (int)(L[j] / rc);
        this->Lcell[j] = L[j] / (int)this->ncells[j];

        this->total_ncells *= this->ncells[j];
        cell_vol *= this->Lcell[j];
    }

    cout << "Has " << this->total_ncells << " cells." << endl;

    this->BLOCK = threads;
    this->GRID = (int)ceil((float)(this->total_ncells) / threads);

    // This factor of 7.5 is a dangerous hard-code
    // Need to figure out how to safely adjust this
    this->CELL_MAX = (int)((float)(ns * cell_vol / V) * 7.5);

    this->cell_count = (int*)calloc(this->total_ncells, sizeof(int));
    this->neigh_count = (int*)calloc(ns, sizeof(int));

    const int alloc_number = this->total_ncells * this->CELL_MAX;
    this->cell_inds = new int[alloc_number];

    const int neigh_alloc_num = ns * this->CELL_MAX;
    this->neigh_inds = new int[neigh_alloc_num];

    int grp_ns = Groups[this->group_index].nsites;
    this->particle_cell = new int[ns];

    cudaMalloc(&this->d_Lcell, Dim * sizeof(float));
    cudaMalloc(&this->d_ncells, Dim * sizeof(int));
    cudaMalloc(&this->d_cell_count, this->total_ncells * sizeof(int));
    cudaMalloc(&this->d_cell_inds, alloc_number * sizeof(int));
    cudaMalloc(&this->d_neigh_count, ns * sizeof(int));
    cudaMalloc(&this->d_neigh_inds, neigh_alloc_num * sizeof(int));
    cudaMalloc(&this->d_particle_cell, ns * sizeof(int));
    

    // Copy data to the device
    cudaMemcpy(this->d_Lcell, this->Lcell, Dim * sizeof(float),
        cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_ncells, this->ncells, Dim * sizeof(int),
        cudaMemcpyHostToDevice);

}

// Calls the device routines to assign particles to cells
void NList::AssignCells() {
    

    d_zero_int_vector<<<this->GRID, this->BLOCK>>>(this->d_cell_count, 
        this->total_ncells);


    int gid = this->group_index;

    d_assign_cells<<<Groups[gid].grid, Groups[gid].block>>>(
        d_x, this->d_cell_count, this->d_cell_inds,
        this->d_Lcell, this->d_ncells, this->CELL_MAX,
        this->d_particle_cell, Groups[gid].d_index, Groups[gid].nsites, Dim);
        
}

void NList::BuildNList() {
    this->NList::AssignCells();
    
    d_zero_int_vector<<<ns_Grid, ns_Block>>>(this->d_neigh_count, ns);

    int gid = this->group_index;

    d_make_nlist << <Groups[gid].grid, Groups[gid].block >> >(d_x, this->rcut, this->d_cell_count, this->d_cell_inds,
        this->d_ncells, this->CELL_MAX, this->d_particle_cell,
        this->d_neigh_count, this->d_neigh_inds, Groups[gid].d_index,
        d_L, d_Lh, Groups[gid].nsites, Dim);
}



// This code currently does 2x as many operations as needed:
// For particle i, when a particle j that is within the cutoff is found
// j is only added to i's neighbor list. Also adding i to j's list would
// require atomic operations, and unclear which is more efficient.
__global__ void d_make_nlist(const float *x, const float rc, const int* cell_ct,
    const int* cell_inds, const int* ncell, const int CMAX,
    const int* part_cell, int* neigh_ct, int* neigh_inds,
    const int* site_list, const float *L, const float *Lh, 
    const int ns, const int D) {
    
    const int list_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (list_id >= ns)
        return;

    const int id = site_list[list_id];
    const int cell_num = part_cell[id];

    // get cell x,y,z indices
    int nn[3];
    d_unstack(cell_num, nn, ncell, D);

    float rc2 = rc * rc;

    // Explicitly assign z-cell for 2D case (it'll be ignored)
    if (D == 2)
        nn[2] = 1;

    // Loop over neighboring cells in 2/3 dimensions
    for (int ix = nn[0] - 1; ix <= nn[0] + 1; ix++) {
        int nix = ix;
        if (nix < 0) nix = ncell[0] - 1;
        else if (nix == ncell[0]) nix = 0;

        for (int iy = nn[1] - 1; iy <= nn[1] + 1; iy++) {
            int niy = iy;
            if (niy < 0) niy = ncell[1] - 1;
            else if (niy == ncell[1]) niy = 0;

            for (int iz = nn[2] - 1; iz <= nn[2] + 1; iz++) {
                int niz = iz;

                if (D == 2)
                    niz = 0;
                else if (D == 3) {
                    if (niz < 0) niz = ncell[2] - 1;
                    else if (niz == ncell[2]) niz = 0;
                }

                int neigh_nn[3] = { nix, niy, niz };

                // Index of the current neighboring cell
                // ncn = neighbor_cell_number
                const int ncn = d_stack_input(neigh_nn, ncell, D); 

                // Loop over the particles in cell 'neigh_cell_num'
                for (int i = 0; i < cell_ct[ncn]; i++) {

                    // nid: neighbord particle id
                    // NOTE: cells only contain particles in this group, so
                    // additional check to make sure nid is in site_list not needed.
                    const int nid = cell_inds[ncn * CMAX + i];

                    // Magnitude of separation between two particles
                    float mdr2 = 0.0f;
                    float dr[3];
                    for (int j = 0; j < D; j++) {
                        dr[j] = x[id * D + j] - x[nid * D + j];
                        if (dr[j] < -Lh[j]) dr[j] += L[j];
                        else if (dr[j] > Lh[j]) dr[j] -= L[j];

                        mdr2 += dr[j] * dr[j];
                    }


                    // Check cut-off, add nid to id's list if close enough
                    if (mdr2 < rc2) {
                        neigh_inds[id * CMAX + neigh_ct[id]] = nid;
                        neigh_ct[id] += 1;
                    }
                }

                // Escape the z-loop if this is a 2D simulation
                if (D == 2)
                    break;
                    

            }// iz = nn[2]-1 : nn[2]+1
        }// iy=nn[1]-1 : nn[1]+1
    }//ix = nn[0]-1:nn[0]+1
}

__global__ void d_assign_cells(float* x, int* cell_ct, int* cell_inds,
    float* cellL, int* ncell, int CMAX, int* part_cell,
    int* site_list, int ns, int D) {

    const int list_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (list_id >= ns)
        return;

    int id = site_list[list_id];

    int cell_id[3] = {0,0,0};
    for (int j = 0; j < D; j++) {
        cell_id[j] = (int)(x[id * D + j] / cellL[j]);
    }

    int cell_num = d_stack_input(cell_id, ncell, D);

    // Assign particle 'id' to its cell
    int ct = atomicAdd(&cell_ct[cell_num], 1);
    cell_inds[cell_num * CMAX + ct] = id;
    
    // Let particle 'id' know its cell
    part_cell[id] = cell_num;
}





NList::NList() {

}
NList::~NList() {

}