add_library(linalg STATIC
    src/matrix.cpp
    src/vector.cpp
)

target_include_directories(linalg
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

