cmake_minimum_required(VERSION 3.20)
project(matrix_lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_subdirectory(lib)
add_subdirectory(app)

enable_testing()
add_subdirectory(tests)