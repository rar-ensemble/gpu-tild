#include "globals.h"

// dr = r1 - r2
float pbc_mdr2(float* r1, float* r2, float* dr) {
	float mdr2 = 0.0f;
	for (int i = 0; i < Dim; i++) {
		dr[i] = r1[i] - r2[i];
		if (dr[i] > Lh[i]) dr[i] -= L[i];
		else if (dr[i] < -Lh[i]) dr[i] += L[i];
		mdr2 += dr[i] * dr[i];
	}
	return mdr2;
}

