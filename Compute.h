//////////////////////////////////////////////
// Rob Riggleman                8/18/2021   //
// Defining compute group to calculate      //
// various quantities at regular intervals. //
// Examples: Average S(k) for a group, avg  //
// density field, etc.                      //
//////////////////////////////////////////////

#ifndef _COMPUTE
#define _COMPUTE

class Compute {
private:
    int GRID;               // Size of thread grid for cell operations
    int BLOCK;              // Num of blocks per thread
    string group_name;      // Name of the group in the n list
    int group_index;        // Index of the group
    string out_name;        // Name of the output file
    string style;           // Compute style (e.g., avg_sk, avg_rho)
    int iparams[5];         // Integer parameters
    int compute_id;
    float fparams[5];       // Float parameters
    float* fstore1;         // Variable to store data
    float* fstore2;         // Variable to store data
    vector<complex<float>> cpx;   // Variable to store data


public:


    int compute_freq;       // Frequency for calling doCompute (timesteps), default 100
    int compute_wait;       // Time to wait before calling doCompute (timesteps), default 0
    string input_command;   // Full text of input command.
    void Initialize(string);
    void allocStorage(void);
    void doCompute(void);
    void writeResults(int);
    Compute();
    ~Compute();

};

#endif
