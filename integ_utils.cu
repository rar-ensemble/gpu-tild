#include "globals.h"

float integ_trapPBC(float* dat) {
	float sum = 0.f;

	for (int i = 0; i < M; i++)
		sum += dat[i];

	for (int i = 0; i < Dim; i++)
		sum *= dx[i];

	return sum;
}