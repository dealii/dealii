# Summary of github-actions workflows

Please keep this feature matrix up-to-date when making changes to the yaml scripts in `./workflows`.

### linux

| Jobname               | Platform   | OS           | Compiler           | Flags                                            | External Dependencies                                    | Miscellaneous                                                                      |
|-----------------------|------------|--------------|--------------------|--------------------------------------------------|----------------------------------------------------------|------------------------------------------------------------------------------------|
| release-serial        | x64        | jammy, noble | g++ 11.4.0, 13.3.0 | -Werror -std=c++20                               |                                                          | container:dealii/dependencies                                                      |
| debug-parallel        | x64, arm64 | jammy, noble | g++ 11.4.0, 13.3.0 | -Werror -std=c++20 -mno-outline-atomics          | OpenMPI, CGAL, HDF5, NetCDF, Metis, PETSc, Trilinos, VTK | container:dealii/dependencies                                                      |
| debug-parallel-tpetra | x64        | jammy, noble | g++ 11.4.0, 13.3.0 | -std=c++20                                       | OpenMPI, p4est, Trilinos                                 | DEAL_II_WITH_64BIT_INDICES=ON, container:dealii/dependencies                       |
| debug-intel-oneapi    | x64        | jammy        | icpx 2026.0.0      | -Werror -Wno-error=tautological-constant-compare | IntelMPI, MKL, TBB                                       | uses:rscohn2/setup-oneapi                                                          |
| debug-cuda-12         | x64        | jammy        | g++ 11.4.0         | -Werror -Wno-non-template-friend                 | OpenMPI, p4est, Kokkos 4.0.01, CUDA 12.3                 | Uses nvcc_wrapper as compiler.                                                     |
| debug-cuda-12-clang   | x64        | jammy        | clang++ 19.1.7     | -std=c++17                                       | OpenMPI, p4est, Kokkos 4.0.01, CUDA 12.3                 |                                                                                    |
| clang-22-modules      | x64        | noble        | clang++ 22.1.8     |                                                  | Kokkos 5.0.1                                             | CMAKE_CXX_STANDARD=23, DEAL_II_WITH_CXX20_MODULE=ON, container:dealii/dependencies |

### macos

| Jobname    | Platform | OS           | Compiler               | Flags   | External Dependencies | Miscellaneous                 |
|------------|----------|--------------|------------------------|---------|-----------------------|-------------------------------|
| serial     | arm64    | macos-15, 26 | clang++ 17.0.0, 21.0.0 | -Werror |                       |                               |
| parallel64 | arm64    | macos-15, 26 | clang++ 17.0.0, 21.0.0 | -Werror | OpenMPI               | DEAL_II_WITH_64BIT_INDICES=ON |

### windows

| Jobname | Platform | OS                 | Compiler                          | Flags          | Miscellaneous                |
|---------|----------|--------------------|-----------------------------------|----------------|------------------------------|
| serial  | x64      | windows-2022, 2025 | MSVC 19.44.35228.0, 19.51.36248.0 | /WX /std:c++20 | FE_EVAL_FACTORY_DEGREE_MAX=2 |
