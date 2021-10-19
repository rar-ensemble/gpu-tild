#include "pair_style.h"

#ifndef _PAIR_FIELDPHASE
#define _PAIR_FIELDPHASE


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class FieldPhase : public PairStyle {
public:
    FieldPhase();
    ~FieldPhase();
    void Initialize_FieldPhase(float, int, int, int, int, int, int);
    void Initialize();
	void CalcForces();

    int n_periods, dir, phase;
};


#endif