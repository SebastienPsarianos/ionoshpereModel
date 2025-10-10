#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <functional>
#include <fstream>

#define MAX_ITERATION_NUM 10000000
#define RES_THRESHOLD  0.01
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
              double x0, 
              double xMax,
              double u_x0,
              double u_xMax,
              SourceTerm sourceTerm ) : _gridSize(gridSize) {
        assert(gridSize > 2);
        assert(xMax > x0);

        _dx = (xMax - x0) / gridSize;
        _dx2 = std::pow(_dx, 2);

        // Select source term
        switch(sourceTerm) {
            case (EXPONENTIAL):
                source = [](double x) {return std::exp(x);};
                break;
            case (CONSTANT):
                source = [](double _) {return A;};
                break;
            case (LAPLACE):
                source = [](double _) {return 0;};
        }

        // Initialize current and previous guess u vectors
        _u[0] = std::vector<double>(gridSize); // now has `gridSize` elements
        _u[1] = std::vector<double>(gridSize); // now has `gridSize` elements
        for (int i = 0; i < gridSize; i++) {
            double t = static_cast<double>(i) / (gridSize - 1); // normalized [0,1]
            double u_guess = (1 - t) * u_x0 + t * u_xMax;        // linear interpolation
            _u[0][i] = u_guess;
            _u[1][i] = u_guess;  // also initialize second buffer
        }

        // Add in boundary conditions
        _u[0][0] = u_x0;
        _u[0][gridSize - 1] = u_xMax;
        _u[1][0] = u_x0;
        _u[1][gridSize - 1] = u_xMax;

        // Set up grid x values and source term values
        for (int i = 0; i < gridSize; i++) {
            _x.push_back(x0 + i * _dx);
            _rho.push_back(source(_x[i]));
        }
    }

    void solve(Algorithm algo) {
        // Use previous iteration for jacobi and current for gauss seidell
        //     we don't use the second array for the gauss seidel to run faster
        for(int itt = 0; itt < MAX_ITERATION_NUM; itt++) {
            if (algo == GAUSS_SEIDEL) {
                for (int x = 1; x < _gridSize - 1; x++) {
                    _u[0][x] = (_u[0][x+1] + _u[0][x-1] - _dx2  * _rho[x]) * 0.5;
                } 
            }
            else {
                 for (int x = 1; x < _gridSize - 1; x++) {
                    _u[1][x] = (_u[0][x+1] + _u[0][x-1] - _dx2  * _rho[x]) * 0.5;
                }
                // Move all of the data from the previous itteration to the next
                std::swap(_u[0], _u[1]);
            }

            // Calculate the residual after each iteration. After the swap, _u[0] will contain the current iteration result regardless of algorithm
            double res = 0;
            for (int x = 1; x < _gridSize - 1; x++) {
                res += std::abs((_u[0][x+1] + _u[0][x-1] - 2 * _u[0][x]) / _dx2 - _rho[x]);
            }
            // Calculate average residual
            res /= (_gridSize - 2);
            _residuals.push_back(res);

            if (res < RES_THRESHOLD) {
                std::cout << "Simulation converged with " << itt << " iterations and a residual of " << res << std::endl;
                return;
            }
        }
        std::cout<<"No Convergence" <<std::endl;
    };

    void printResult(Algorithm algo, SourceTerm source) {
        std::string fname = (algo == GAUSS_SEIDEL ? "Gauss" : "Jacobi");
        fname += "_";
        switch(source) {
            case EXPONENTIAL: fname += "Exp"; break;
            case CONSTANT: fname += "Const"; break;
            case LAPLACE: fname += "Lap"; break;
        }
        fname += "_" + std::to_string    (_gridSize) + ".dat";

        std::ofstream dataFile(fname);
        if (!dataFile) {
            std::cerr << "Error opening file for output!\n";
            return;
        }

        for (int i = 0; i < _gridSize; i++) {
            dataFile << _x[i] << " " << _u[0][i] << "\n";
        }

        dataFile << "RESIDUALS" << "\n";

        for (long unsigned int i = 0; i < _residuals.size(); i++) {
            dataFile << i << " " << _residuals[i] << "\n"; 
        }

        dataFile.close();
        std::cout << "Results saved to " << fname << std::endl;
    }

private: 
    std::vector<double> _x;
    std::array<std::vector<double>, 2> _u;
    std::vector<double> _rho;
    std::vector<double> _residuals;
    double _dx;
    double _dx2;
    int _gridSize;


    std::function<double(double)> source;
};

int main(int argc, char* argv[]) {
    constexpr int GRID_SIZE = 400;
    constexpr double X0 = 0.0;
    constexpr double XMAX = 1.0;
    constexpr double U_X0 = 0.0;
    constexpr double U_XMAX = 1.0;

    Solver s1(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, CONSTANT);
    s1.solve(GAUSS_SEIDEL);
    s1.printResult(GAUSS_SEIDEL, CONSTANT);

    Solver s2(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, CONSTANT);
    s2.solve(JACOBI);
    s2.printResult(JACOBI,CONSTANT);

    Solver s3(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, EXPONENTIAL);
    s3.solve(GAUSS_SEIDEL);
    s3.printResult(GAUSS_SEIDEL,EXPONENTIAL);

    Solver s4(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, EXPONENTIAL);
    s4.solve(JACOBI);
    s4.printResult(JACOBI,EXPONENTIAL);

    Solver s5(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, LAPLACE);
    s5.solve(GAUSS_SEIDEL);
    s5.printResult(GAUSS_SEIDEL,LAPLACE);

    Solver s6(GRID_SIZE, X0, XMAX, U_X0, U_XMAX, LAPLACE);
    s6.solve(JACOBI);
    s6.printResult(JACOBI,LAPLACE);

    return 0;
}
