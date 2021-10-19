///////////////////////////////////////////////
// Rob Riggleman                7/9/2021     //
// Defining group class to eventually enable //
// functions to operate on specific groups.  //
///////////////////////////////////////////////

#ifndef _GROUPS
#define _GROUPS


class Group {
private:
public:

    int nsites;             // Number of sites in this group
    int* index;             // Site indices of this group
    int* d_index;           // Site indices for this group on device
    int grid;               // block grid size for operations on this group on the GPU
    int block;              // size of the blocks for operations on this group on GPU
    string command_line;    // Stores the full command given in input file to create group
    string name;            // Text name of this group
    void Initialize(void);  // Allocates memory for lists
    int has_nlist;          // Flag to indicate whether this group has a neighbor list
    int nlist_num;          // Number associated with the neighbor list
    Group();
    ~Group();
};

#endif