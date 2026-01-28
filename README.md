# matrix

Small C++ matrix library with basic linear-algebra utilities (RREF used by tests and the example app).

## Features
- Matrix type with initializer-list construction.
- rref() to compute reduced-row echelon form.
- Header-only public include at `lib/include/matrix/matrix.hpp`.
- Built as `matrix` static library and used by `app` and `tests`.

## Project layout
- lib/ — library sources and public headers (`lib/include/matrix/matrix.hpp`)
- app/ — example application using the library
- tests/ — Google Test unit tests (tests/test_matrix.cpp)
- CMakeLists.txt — top-level build configuration

## Build (macOS / POSIX)
```bash
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
cmake --build . --config Debug -- -j
ctest --output-on-failure
# or run the test binary directly:
./tests/test_matrix
```

## Requirements
- CMake >= 3.20
- C++20 compiler
- Google Test (findable by CMake). On macOS you can install with Homebrew if needed:
```bash
brew install googletest
```

## Usage
Include the public header and link against the `matrix` target:
```cpp
#include "matrix/matrix.hpp"

Matrix A({{1,2,3},{4,5,6}});
A.rref();
double v = A(0,0);
```

## Notes
- Tests use Google Test and are discovered via `gtest_discover_tests`.
- Headers are provided under `lib/include/matrix/`.