//////////////////////////////////////////////
// Rob Riggleman                7/21/2021   //
// Defining integrator class that will call //
// particle integration routines based on   //
// groups.                                  //
//////////////////////////////////////////////

#include "group.h"

#ifndef _INTEGRATOR
#define _INTEGRATOR

class Integrator {
private:
public:
    int group_index;  // index of group to be integrated
    string group_name;// name of the group to be integrated
    string name;      // name of the integrator to be used
    void Initialize(string, string);
    Integrator();     // constructor
    ~Integrator();    // destructor
    void Integrate_1(void); // Calls the pre-force integrator
    void Integrate_2(void); // Calls the pre-force integrator
};

#endif