#include "globals.h"
#include "Extraforce.h"
#include <sstream>

__global__ void d_ExtraForce_Langevin(float*, const float*,
    const float, const float, const int*, const int,
    const int, curandState*);

__global__ void d_ExtraForce_midpush(float*, const float*,
        const float, const int*, const float*, const int, const int);


void ExtraForce::AddExtraForce() {

    const int gid = this->group_index;
    const int nlid = Groups[gid].nlist_num;
    const int GRID = Groups[gid].grid;
    const int BLOCK = Groups[gid].block;

    if (this->style == "langevin") {

        // Noise in Langevin equation is sqrtf( 2.0 * gamma * delt);
        // We must divide by delt here b/c the VV algo multiplies f by delt
        const float lang_noise = sqrtf(this->params[0]) * noise_mag / delt;

        d_ExtraForce_Langevin<<<GRID, BLOCK>>>(d_f, d_v, lang_noise,
            this->params[0], Groups[gid].d_index, Groups[gid].nsites,
            Dim, d_states);

    }


    // Pushes molecules towards center of the box in the Dim-1 direction
    else if ( this->style == "midpush" ) {
        d_ExtraForce_midpush<<<GRID, BLOCK>>>(d_f, d_x, 
            this->params[0], Groups[gid].d_index, d_Lh, Groups[gid].nsites, Dim);
    }
}


// This initialize routine will be called from read_input
// It will parse the remainder of the input command to determine
// style of 'ExtraForces' and any relevant parameters.
void ExtraForce::Initialize(string input_cmd) {

    this->command_line = input_cmd;
    istringstream iss(input_cmd);
    string word;
    iss >> word; // skip extraforces keyword


    // Group on which this ExtraForce acts
    iss >> word;
    this->group_index = -1;
    for (int i = 0; i < n_groups; i++) {
        if (Groups.at(i).name == word) {
            this->group_index = i;
            this->group_name = word;
        }
    }

    if (this->group_index == -1) {
        char death[120];
        sprintf(death, "Group %s not found to apply ExtraForce!", word.c_str());
        die(death);
    }



    // ExtraForces style
    iss >> word;



    // For extraforces style "langevin":
    // params[0] is the drag coefficient
    if (word == "langevin") {
        this->style = "langevin";

        string tofloat;
        iss >> tofloat;
        this->params[0] = stof(tofloat);
    }

    // Extraforces midpush style
    // params[0] is force magnitude
    else if (word == "midpush") {
        this->style = "midpush";

        string tofloat;

        iss >> tofloat;
        this->params[0] = stof(tofloat);
    }// if extraforces == dpd

    // Invalid extraforces command
    else {
        char death[120];
        sprintf(death, "extraforces %s not a valid ExtraForce!", word.c_str());
        die(death);
    }
}


string ExtraForce::PrintCommand() {
    return this->command_line;
}

ExtraForce::ExtraForce() {
}

ExtraForce::~ExtraForce() {
}
