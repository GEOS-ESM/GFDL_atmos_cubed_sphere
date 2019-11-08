#include <iostream>
#include <Kokkos_Core.hpp>

extern "C" { 
  void kernel1_orig(int nrows, int ncols, const float* yfx, const float* fy2, float* fyy);
  void kernel1(int nrows, int ncols, const float* yfx, const float* fy2, float* fyy);
} 

void kernel1_orig(int nrows, int ncols, const float* yfx, const float* fy2, float* fyy) {
  for (auto index = 0; index < nrows*ncols; index++)
    fyy[index] = yfx[index] * fy2[index];
}

void kernel1(int nrows, int ncols, const float* yfx, const float* fy2, float* fyy) {
  typedef Kokkos::View<float*, Kokkos::MemoryUnmanaged> t_1d;
  typedef Kokkos::View<const float*, Kokkos::MemoryUnmanaged> t_1d_const;

  // Wrap incoming 2D arrays in Kokkos::View objects
  t_1d_const yfx_k(yfx, nrows*ncols);
  t_1d_const fy2_k(fy2, nrows*ncols);
  t_1d fyy_k(fyy, nrows*ncols);

  Kokkos::parallel_for(nrows*ncols, [=] (const size_t index) {
      fyy_k[index] = yfx_k[index] * fy2_k[index];
    });
}
