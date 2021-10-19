extern __device__ float d_get_k(int, float*, const float*,
    const int*, const int);

extern __device__ float d_pbc_mdr2(float*, float*, float*, const float*,
    const float*, const int);

extern __device__ void d_get_r(const int, float*, const int*,
    const float*, const int);

extern __device__ void d_unstack(int, int*, const int*, const int);