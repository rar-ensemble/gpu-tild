#include "pair_style.h"

#ifndef _PAIR_CHARGES
#define _PAIR_CHARGES

class Charges : public PairStyle {
public:
	void CalcCharges();
	void Initialize();
};

#endif