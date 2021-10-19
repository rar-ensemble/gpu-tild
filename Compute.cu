#include "globals.h"
#include "timing.h"
#include "Compute.h"
#include <sstream>

void write_kspace_data(const char*, complex<float>*);
__global__ void d_multiplyComplex(cufftComplex*, cufftComplex*,
    cufftComplex*, int);
__global__ void d_prepareDensity(int, float*, cufftComplex*, int);

Compute::Compute() {};
Compute::~Compute() {};

void Compute::Initialize(string command) {
    this->input_command = command;
    this->compute_id = n_computes - 1;

    istringstream iss(command);
    string word;

    // Put the first word, which should be "compute"
    iss >> word;

    // Detect compute style
    iss >> word;
    if (word == "avg_sk") {
        
        // For avg_sk, will calculate the average S(k) for one type of particle.
        // iparams[0] = particle type
        // iparams[1] will track the number of calls to doCompute()

        this->style = "avg_sk";


        iss >> word;
        this->iparams[0] = stoi(word) - 1;    // Store relevant type, shifting by one to zero index.
        if (this->iparams[0] >= ntypes)
            die("Invalid type to calculate avg_sk");
        
        this->iparams[1] = 0;


        // Set defaults for optional arguments
        this->compute_wait = 0;
        this->compute_freq = 100;

        
        // Check for optional arguments
        while (!iss.eof()) {
            iss >> word;
            if (word == "freq") {
                iss >> word;
                this->compute_freq = stoi(word);
            }
            else if (word == "wait") {
                iss >> word;
                this->compute_wait = stoi(word);
            }
        }
        cout << "  Calculating <S(k)> for component " << iparams[0] + 1 << " every " << this->compute_freq <<
            " steps after " << this->compute_wait << " steps have passed." << endl;
        
    }// word == avg_sk

    else {
        die("Undefined compute style!");
    }
}


void Compute::allocStorage() {

    if ( this->style == "avg_sk" ) {

        this->cpx.resize(M);
        for (int i = 0; i < M; i++) 
            this->cpx.at(i) = 0.0f;
        
        cout << " this->cpx has initial size " << this->cpx.size() << " and capacity " << this->cpx.capacity() << endl;
    }

}
void Compute::doCompute() {
    compute_t_in = int(time(0));

    if (this->style == "avg_sk") {

        // Extract the density of the relevant type
        d_prepareDensity<<<M_Grid, M_Block>>>(this->iparams[0], d_all_rho, d_cpx1, M);
        check_cudaError("Compute->doCompute.style = avg_sk prepare density");

        // fourier from d_cpx1 to d_cpx2 forward
        cufftExecC2C(fftplan, d_cpx1, d_cpx2, CUFFT_FORWARD);
        check_cudaError("Compute->doCompute.style = avg_sk cufftExec");


        // Multiply by the complex conjugate and scale by 1/M
        // Store it in d_cpx1 as the values inside are not needed at this point
        d_multiplyComplex<<<M_Grid, M_Block>>> (d_cpx2, d_cpx2, d_cpx1, M);
        check_cudaError("Compute->doCompute.style = avg_sk multiplyComplex");


        // Copy data to host and store
        // NOTE: this should probably be stored on the device and only 
        // communicated when writing, but may be OK for now.
        cudaMemcpy(cpx1, d_cpx1, M * sizeof(cufftComplex), cudaMemcpyDeviceToHost);
    
        for (int i = 0; i < M; i++)
            this->cpx.at(i) += cpx1[i].x + I * cpx1[i].y;

        this->iparams[1] += 1;
    }

    compute_t_out = int(time(0));
    compute_tot_time += compute_t_out - compute_t_in;
    
}


// Results are going to be written at log_freq intervals
void Compute::writeResults(int compute_id) {


    if (this->style == "avg_sk") {
        for (int i = 0; i < M; i++) {
            if (this->iparams[1] > 0)
                k_tmp[i] = this->cpx[i] / float(this->iparams[1]);
            else
                k_tmp[i] = 0.0f;
        }
        char nm[50];

        // compute_id is used in the name instead of "type" in case multiple 
        // computes operate on the same type
        sprintf(nm, "avg_sk_%d.dat", compute_id);
        write_kspace_data(nm, k_tmp);
    }
}


