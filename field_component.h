#ifndef _FIELD_COMP
#define _FIELD_COMP

class FieldComponent {
public:
    float* rho, ** force, * d_rho, * d_force;

    //__global__ void ZeroDeviceGradient(int);
    void ZeroGradient();
    void Initialize(int);
    FieldComponent();
    ~FieldComponent();
};

#endif  