#include <iostream>
#include "matrix/matrix.hpp"
#include <string>

int main(int argc, char *argv[])
{
    if (argc > 3)
    {
        size_t rows = std::stoul(argv[1]);
        size_t cols = std::stoul(argv[2]);
        Matrix A(rows, cols);

        size_t expected_args = 3 + rows * cols;
        if (argc != expected_args)
        {
            std::cerr << "Expected " << expected_args - 1 << " values for a "
                      << rows << "x" << cols << " matrix, but got "
                      << argc - 3 << ".\n";
            return 1;
        }

        size_t arg_index = 3;
        for (size_t r = 0; r < rows; ++r)
        {
            for (size_t c = 0; c < cols; ++c)
            {
                A(r, c) = std::stod(argv[arg_index++]);
            }
        }
        std::cout << "Input matrix:\n"
                  << A << "\n";

        A.rref();
        std::cout << "RREF of the matrix:\n"
                  << A << "\n";
    }
    else
    {
        std::cerr << "Usage: " << argv[0]
                  << " <rows> <cols> <matrix values...>\n";
        return 1;
    }
}