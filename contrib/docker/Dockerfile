FROM dealii/dependencies:focal

LABEL maintainer="luca.heltai@gmail.com"

USER root
RUN cd /usr/src \
    && git clone https://github.com/dealii/dealii \
    && cd dealii && mkdir build && cd build \
    && cmake -GNinja \
    -DCMAKE_PREFIX_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi;/usr/include/hdf5/openmpi" \
    -DDEAL_II_COMPILE_EXAMPLES=OFF \
    -DDEAL_II_COMPONENT_PYTHON_BINDINGS=ON \
    -DDEAL_II_WITH_MPI=ON \
    -DCMAKE_CXX_FLAGS="-std=c++17" \
    .. \  
    && ninja install \
    && cd ../ && rm -rf .git build

USER $USER
WORKDIR $HOME
