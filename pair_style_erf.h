#include "pair_style.h"

#ifndef _PAIR_ERF
#define _PAIR_ERF


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class Erf : public PairStyle {

public:
    Erf();
    ~Erf();
    void Initialize_Erf(float, float, float, int, int, int);
    void Initialize();
};


#endif