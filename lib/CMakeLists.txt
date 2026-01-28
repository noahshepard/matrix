add_library(matrix STATIC
    src/Matrix.cpp
)

target_include_directories(matrix
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)