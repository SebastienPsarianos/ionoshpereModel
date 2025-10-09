#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <functional>

#define MAX_ITTERATION_NUM 1000
#define RES_THRESHOLD  4
#define A 10

enum SourceTerm {
    EXPONENTIAL,
    CONSTANT,
    LAPLACE
};
enum Algorithm {
    GAUSS_SEIDEL,
    JACOBI
};

class Solver {
    public:
    Solver (int gridSize, 
              float x0, 
              float xMax,
              float u_x0,
              float u_xMax,
              SourceTerm sourceTerm ) {
        assert(gridSize > 2);
        assert(xMax > x0);

        _dx = (xMax - x0) / gridSize;
        _dx2 = std::pow(_dx, 2);

        // Select source term
        switch(sourceTerm) {
            case (EXPONENTIAL):
                source = [](float x) {return std::exp(x);};
                break;
            case (CONSTANT):
                source = [](float _) {return A;};
                break;
            case (LAPLACE):
                source = [](float _) {return 0;};
        }

        // Initialize current and previous guess u vectors
        _u[0] = std::vector<float>(gridSize, 0.0f);
        _u[1] = std::vector<float>(gridSize, 0.0f);

        // Add in boundary conditions
        _u[0][0] = u_x0;
        _u[0][gridSize - 1] = u_xMax;
        _u[1][0] = u_x0;
        _u[1][gridSize - 1] = u_xMax;

        // Set up grid x values and source term values
        for (int i = 0; i < gridSize; i++) {
            _x.push_back(i * _dx);
            _rho.push_back(source(_x[i]));
        }
    }

private: 
    std::vector<float> _x;
    std::array<std::vector<float>, 2> _u;
    std::vector<float> _rho;
    float _dx;
    float _dx2;

    std::function<float(float)> source;

    void solve(Algorithm algo) {
        // Use previous iteration for jacobi and current for gauss seidell
        int itteration_offset;
        switch(algo) {
            case GAUSS_SEIDEL:
                itteration_offset = 1;
                break;
            case JACOBI:
                itteration_offset = 0;
                break;
        }

        for(int itt = 0; itt < MAX_ITTERATION_NUM; itt++) {
            for (int x = 1; x < _u[0].size() - 1; x++) {
                _u[1][x] = (_u[0][x+1] + _u[itteration_offset][x-1] - _dx2  * _rho[x]) * 0.5;
            }

            // After each iteration, set u[0] as the previous iteration
            //     and calculate the residual
            float res = 0;
            for (int x = 1; x < _u.size() - 1; x++) {
                res += std::abs(_u[1][x+1] + _u[1][x-1] - 2 * _dx2 * _rho[x]);
                _u[0][x] = _u[1][x];
            }

            if (res < RES_THRESHOLD) {
                std::cout << "Simulation converged with " << itt << " iterations and a residual of " << res << std::endl;
                break;
            }
        }
    };
};

int main(int argc, char* argv[]) {

}
