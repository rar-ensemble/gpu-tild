#include "globals.h"

void die(const char* msg) {
	cout << msg << endl;
	exit(1);
}

void check_cudaError(const char* tag) {
	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while %s ran\n", cudaGetErrorString(cudaReturn), tag);
		die(cherror);
	}
}