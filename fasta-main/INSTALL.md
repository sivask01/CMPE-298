# Installing Fasta via conda

The recommended way to install Fasta is through [conda](https://docs.conda.io).
Stable releases are pushed regularly to the pytorch conda channel.

The `fasta` conda package is currently available on Linux, OSX, and
Windows. 

To install the latest stable release:

``` shell
$ conda install -c pytorch fasta
```

# Building from source

Fasta can be built from source using CMake. Fasta can be installed on x86_64 machines running Linux or OSX. It may also successfully run on other platforms or architectures, however it is not being actively tested for those.

The basic requirements are:
- a C++20 compiler (with support for OpenMP support version 2 or higher)
- the `fmt` library
- for the python bindings:
  - python 3,
  - numpy,
  - and swig.

For Debian, install required libraries with:

``` shell
$ sudo apt install libfmt-dev swig
$ python -m pip install numpy
```

## Step 1: Invoking CMake

``` shell
$ cmake -B build .
```

This generates the system-dependent configuration/build files in the `build/`
subdirectory.

Several options can be passed to CMake, among which:
- general options:
  - `-DFASTA_ENABLE_PYTHON=OFF` in order to disable building python bindings
  (possible values are `ON` and `OFF`),
  - `-DFASTA_ENABLE_WARNINGS=ON` in order to enable compilation warnings
  (possible values are `ON` and `OFF`),
  - `-DBUILD_TESTING=OFF` in order to disable building C++ tests,
  - `-DBUILD_SHARED_LIBS=ON` in order to build a shared library (possible values
  are `ON` and `OFF`),
- optimization-related options:
  - `-DCMAKE_BUILD_TYPE=Release` in order to enable generic compiler
  optimization options (enables `-O3` on gcc for instance; possible values are
  empty, `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`),
- python-related options:
  - `-DPython_EXECUTABLE=/path/to/python3.9` in order to build a python
  interface for a different python than the default one (see
  [CMake docs](https://cmake.org/cmake/help/latest/module/FindPython.html)).

## Step 2: Invoking Make

``` shell
$ make -C build -j fasta
```

This builds the C++ library (`libfasta.a` by default, and `libfasta.so` if
`-DBUILD_SHARED_LIBS=ON` was passed to CMake).

The `-j` option enables parallel compilation of multiple units, leading to a
faster build, but increasing the chances of running out of memory, in which case
it is recommended to set the `-j` option to a fixed value (such as `-j4`).

## Step 3: Building the python bindings (optional)

``` shell
$ make -C build -j swigfasta
$ (cd build/fasta/python && python setup.py install)
```

The first command builds the python bindings for Fasta, while the second one
generates and installs the python package.

## Step 4: Installing the C++ library and headers (optional)

``` shell
$ make -C build install
```

This will make the compiled library (either `libfasta.a` or `libfasta.so` on
Linux) available system-wide, as well as the C++ headers. This step is not
needed to install the python package only.

Optionally, you can specify a prefix to have the library installed in an alternate location than the system default.

``` shell
$ make -C build install --prefix="/path/to/install/location"
```

## Step 5: Testing (optional)

### Running the C++ test suite

To run the whole test suite, make sure that `cmake` was invoked with
`-DBUILD_TESTING=ON`, and run:

``` shell
$ cmake --build build --target fasta_test
$ make -C build test
```

### Running the python test suite

``` shell
$ (cd build/fasta/python && python setup.py build)
$ PYTHONPATH="$(ls -d ./build/fasta/python/build/lib*/)" pytest tests/python/test_*.py
```

### Basic example

A basic usage example is available in
[`demos/demo_fasta_indexing.cpp`](https://github.com/facebookresearch/fasta/blob/main/demos/demo_fasta_indexing.cpp).

It creates a small index performs some searches. A normal runtime
is around 10s. With a fast machine, it runs in <2s.

It can be built with
``` shell
$ make -C build demo_fasta_indexing
```
and subsequently ran with
``` shell
$ ./build/demos/demo_fasta_indexing
```

# Development local building

To compile the program with debuging symbols for development, you must enable the `Debug` CMake build type. Note that code converage symbols will also be added to the executable as they are required for the `csv-parser` library. To build a debug version of the library and executable in a directory called `debug`, execute the following.

``` shell
mkdir debug
cmake -DCMAKE_BUILD_TYPE=Debug -B debug
make -C debug exec
```

The `llvm` library is needed on MacOS and the `lcov` library on Linux to measure code coverage. On Debian, you can install `lcov` as follows.

``` shell
sudo apt install lcov
```

# Acknowledgment

The structure of the project, documentation, and many of the overall settings were heavily inspired from the [FAISS](https://github.com/facebookresearch/faiss) project from Facebook Research.
