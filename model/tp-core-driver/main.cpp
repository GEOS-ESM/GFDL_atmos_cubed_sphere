#include <iostream>
#include <cstdlib>  // atoi
#include <cassert>

#include <Kokkos_Core.hpp>

// Implemented in Fortran
extern "C" {
  void run_driver(int*, int*);
}

int main(int argc, char** argv) {
  assert(argc == 3);
  int resolution = atoi(argv[1]);
  int n_iterations = atoi(argv[2]);

  run_driver(&resolution, &n_iterations);
  return 0;
}
