FROM dealii/dependencies:focal

LABEL maintainer="luca.heltai@gmail.com"

USER root
RUN cd /usr/src \
    && git clone https://github.com/dealii/dealii \
    && cd dealii && mkdir build && cd build \
    && cmake -GNinja \
    -DCMAKE_PREFIX_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi;/usr/include/hdf5/openmpi" \
    -DDEAL_II_ALLOW_AUTODETECTION=OFF \
    -DDEAL_II_COMPILE_EXAMPLES=OFF \
    -DDEAL_II_COMPONENT_PYTHON_BINDINGS=ON \
    -DCMAKE_CXX_FLAGS="-std=c++17" \
    -DDEAL_II_WITH_64BIT_INDICES=OFF \
    -DDEAL_II_WITH_ADOLC=ON \
    -DDEAL_II_WITH_ARBORX=OFF \
    -DDEAL_II_WITH_ARPACK=ON \
    -DDEAL_II_WITH_ASSIMP=ON \
    -DDEAL_II_WITH_BOOST=ON \
    -DDEAL_II_WITH_CGAL=ON \
    -DDEAL_II_WITH_COMPLEX_VALUES=ON \
    -DDEAL_II_WITH_CUDA=OFF \
    -DDEAL_II_WITH_GINKGO=OFF \
    -DDEAL_II_WITH_GMSH=ON \
    -DDEAL_II_WITH_GSL=ON \
    -DDEAL_II_WITH_HDF5=ON \
    -DDEAL_II_WITH_KOKKOS=ON \
    -DDEAL_II_WITH_LAPACK=ON \
    -DDEAL_II_WITH_METIS=ON \
    -DDEAL_II_WITH_MPI=ON \
    -DDEAL_II_WITH_MUPARSER=ON \
    -DDEAL_II_WITH_OPENCASCADE=ON \
    -DDEAL_II_WITH_P4EST=ON \
    -DDEAL_II_WITH_PETSC=ON \
    -DDEAL_II_WITH_SCALAPACK=ON \
    -DDEAL_II_WITH_SLEPC=ON \
    -DDEAL_II_WITH_SUNDIALS=ON \
    -DDEAL_II_WITH_SYMENGINE=ON \
    -DDEAL_II_WITH_TASKFLOW=OFF \
    -DDEAL_II_WITH_TBB=ON \
    -DDEAL_II_WITH_TRILINOS=ON \
    -DDEAL_II_WITH_UMFPACK=ON \
    -DDEAL_II_WITH_ZLIB=ON \
    .. \
    && ninja -j 2 install \
    && cd ../ && rm -rf .git build

USER $USER
WORKDIR $HOME
