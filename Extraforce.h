//////////////////////////////////////////////
// Rob Riggleman                7/26/2021   //
// Class for adding additional forces.      //
// Will initially be designed to add the DPD//
// dissipation and friction forces, but     //
// should be extensible to make groups      //
// interact with walls and the like.        //
//////////////////////////////////////////////

#ifndef _EXTRAFORCE
#define _EXTRAFORCE

class ExtraForce {
private:
    int group_index;    // index of group on which this acts
    string group_name;  // Name of the group on which this acts
    string style;    // Style of the extraforce
    float params[5];    // Parameters for any of the extra forces routine
                        //  Defined in more detail in ExtraForce.cu for each style
    string command_line;// Full line of the command from input file


public:
    void Initialize(string);
    void AddExtraForce(void);
    string PrintCommand(void);
    ExtraForce();
    ~ExtraForce();
};
#endif
