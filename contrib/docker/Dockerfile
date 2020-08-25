FROM dealii/dependencies:focal

LABEL maintainer="luca.heltai@gmail.com"

USER root
RUN git clone https://github.com/dealii/dealii \
    && cd dealii && mkdir build && cd build \
    && cmake -GNinja \
    -DCMAKE_PREFIX_PATH="/usr/lib/x86_64-linux-gnu/hdf5/openmpi;/usr/include/hdf5/openmpi" \
    -DDEAL_II_WITH_MPI=ON \
    -DDEAL_II_WITH_SIMPLEX_SUPPORT=ON \
    -DCMAKE_INSTALL_PREFIX=/usr \
    -DDEAL_II_COMPILE_EXAMPLES=OFF \
    .. \  
    && ninja install \
    && cd ../../ && rm -rf dealii

USER $USER
WORKDIR $HOME

