FROM ubuntu

# Install necessary build tools
RUN add-apt-repository -y 'deb http://mirrors.kernel.org/ubuntu hirsute main universe'
RUN apt-get update && apt-get install -y gcc-11 make swig3.0 libfmt-dev libspdlog-dev
RUN snap install -y cmake --classic
RUN ln -s /snap/bin/cmake /usr/bin/

# Install necesary headers/libs
RUN apt-get install -y python-devel numpy

COPY . /opt/fasta

WORKDIR /opt/fasta

# Configure and make
RUN cmake -B -DFASTA_ENABLE_PYTHON=ON -DBUILD_TESTING=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release build .
RUN make -C build -j $(nproc) fasta
RUN make -C build -j swigfasta
RUN (cd build/faiss/python && python setup.py install)
RUN make -C build install
RUN make -C build demo_fasta_indexing && ./build/demos/demo_fasta_indexing
