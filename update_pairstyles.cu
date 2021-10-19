#include "globals.h"

void update_pairstyles() {
	
	for (Gaussian& Iter : Gausses) {
		Iter.Update();
	}

	for (Erf& Iter : Erfs) {
		Iter.Update();
	}

	for (FieldPhase& Iter : Fields ) {
		Iter.Update();
	}

	for (GaussianErf& Iter : GaussErfs) {
		Iter.Update();
	}

}