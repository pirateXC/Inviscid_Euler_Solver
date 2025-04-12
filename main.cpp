#include "gnuplot-iostream.h"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

// Structure definition using Eigen matrices for grid data.
struct Grid {
    int nx, ny;
    Eigen::MatrixXd x;
    Eigen::MatrixXd y;
};

bool readGridFile(const std::string &filename, Grid &grid);
void plotGrid(const Grid &grid);

int main() {
    Grid grid;
    if (!readGridFile("g641x065uf.dat", grid)) {
        std::cerr << "Error reading grid file.\n";
        return 1;
    }
    plotGrid(grid);
    return 0;
}

bool readGridFile(const std::string &filename, Grid &grid) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    // Read header and parse grid dimensions.
    std::string header;
    std::getline(file, header);
    
    size_t pos = header.find("i=");
    if (pos != std::string::npos) {
        pos += 2;
        size_t end = header.find(",", pos);
        std::string iVal = header.substr(pos, end - pos);
        grid.nx = std::stoi(iVal);
    } else {
        std::cerr << "Unable to find 'i=' in header." << std::endl;
        return false;
    }
    
    pos = header.find("j=");
    if (pos != std::string::npos) {
        pos += 2;
        size_t end = header.find_first_of(" \t", pos);
        std::string jVal = header.substr(pos, end - pos);
        grid.ny = std::stoi(jVal);
    } else {
        std::cerr << "Unable to find 'j=' in header." << std::endl;
        return false;
    }
    
    std::cout << "Parsed grid dimensions: nx = " << grid.nx << ", ny = " << grid.ny << "\n";

    // Allocate Eigen matrices to hold the grid points.
    grid.x.resize(grid.nx, grid.ny);
    grid.y.resize(grid.nx, grid.ny);

    // Loop through the file and populate the matrices.
    for (int i = 0; i < grid.nx; i++) {
        for (int j = 0; j < grid.ny; j++) {
            double xVal, yVal;
            char comma;
            if (!(file >> xVal >> comma >> yVal)) {
                std::cerr << "Error reading data at position [i=" << i << ", j=" << j << "].\n";
                return false;
            }
            grid.x(i, j) = xVal;
            grid.y(i, j) = yVal;
        }
    }
    return true;
}

void plotGrid(const Grid &grid) {
    // Create a Gnuplot object with the persistent flag.
    Gnuplot gp;
    
    // --- Plot horizontal lines ---
    // Each horizontal line corresponds to a fixed row in the grid matrices.
    for (int i = 0; i < grid.nx; i++) {
        std::vector<std::pair<double, double>> line;
        line.reserve(grid.ny);
        for (int j = 0; j < grid.ny; j++) {
            line.emplace_back(grid.x(i, j), grid.y(i, j));
        }
        // For the first line use "plot", then "replot" for subsequent ones.
        if (i == 0)
            gp << "plot '-' with lines lt rgb 'black' notitle\n";
        else
            gp << "replot '-' with lines lt rgb 'black' notitle\n";
        gp.send1d(line);
    }
    
    // --- Plot vertical lines ---
    // Each vertical line corresponds to a fixed column in the grid matrices.
    for (int j = 0; j < grid.ny; j++) {
        std::vector<std::pair<double, double>> line;
        line.reserve(grid.nx);
        for (int i = 0; i < grid.nx; i++) {
            line.emplace_back(grid.x(i, j), grid.y(i, j));
        }
        gp << "replot '-' with lines lt rgb 'black' notitle\n";
        gp.send1d(line);
    }
}
