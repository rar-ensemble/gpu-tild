#include "pair_style.h"

#ifndef _PAIR_GAUSSIAN_ERF
#define _PAIR_GAUSSIAN_ERF


// This class is of type PairStyle and inherits all PairStyle
// routines. PairStyle is initialized first, then the Gaussian
// initialization is performed.
class GaussianErf : public PairStyle {
public:
    GaussianErf();
    ~GaussianErf();
    void Initialize_Gaussian_Erf(float, float, float, int, int, int);
    void Initialize();
};


#endif
