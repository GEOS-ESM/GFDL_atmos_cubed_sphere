#include <iostream>
#include <cstdlib>  // atoi
#include <cassert>

#include <Kokkos_Core.hpp>

// Implemented in Fortran
extern "C" {
  void run_driver(int*, int*);
}

void usage(char* program_name) {
  std::cout << "Usage: " << program_name << " <res> <niters>" << std::endl;
}

void error(char* errmsg) {
  std::cerr << "ERROR: " << errmsg << std::endl;
  std::exit(1);
}

int main(int argc, char** argv) {
  if (argc != 3) {
    usage(argv[0]);
    error("Incorrect number of arguments");
  }
  int resolution = atoi(argv[1]);
  int n_iterations = atoi(argv[2]);

  // Initialize Kokkos  
  Kokkos::InitArguments args;
  args.num_threads = 4;
  args.num_numa = 2;
  Kokkos::initialize(args);
  
  // Call driver
  run_driver(&resolution, &n_iterations);

  // Wrap up
  Kokkos::finalize();
  return 0;
}
