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

    std::cout << "x(0,0) = " << grid.x(0,0) << ", y(0,0) = " << grid.y(0,0) << "\n";
    std::cout << "x(640,0) = " << grid.x(640,0) << ", y(640,0) = " << grid.y(640,0) << "\n";
    std::cout << "x(0,64) = " << grid.x(0,64) << ", y(0,64) = " << grid.y(0,64) << "\n";



    //plotGrid(grid);
    return 0;
}

bool readGridFile(const std::string &filename, Grid &grid) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return false;
    }

    // Read header line.
    std::string header;
    std::getline(file, header);
    
    // Parse grid dimensions from the header (assumes format "ZONE i=641, j=65").
    size_t pos = header.find("i=");
    if (pos != std::string::npos) {
        pos += 2;
        size_t end = header.find(",", pos);
        grid.nx = std::stoi(header.substr(pos, end - pos));
    } else {
        std::cerr << "Unable to find 'i=' in header.\n";
        return false;
    }
    
    pos = header.find("j=");
    if (pos != std::string::npos) {
        pos += 2;
        size_t end = header.find_first_of(" \t", pos);
        grid.ny = std::stoi(header.substr(pos, end - pos));
    } else {
        std::cerr << "Unable to find 'j=' in header.\n";
        return false;
    }
    
    std::cout << "Parsed grid dimensions: nx = " << grid.nx << ", ny = " << grid.ny << "\n";

    // Prepare vectors to hold the data.
    std::vector<double> xvals;
    std::vector<double> yvals;
    xvals.reserve(grid.nx * grid.ny);
    yvals.reserve(grid.nx * grid.ny);
    
    double xVal, yVal;
    char comma;
    while (file >> xVal >> comma >> yVal) {
        xvals.push_back(xVal);
        yvals.push_back(yVal);
    }
    
    if (static_cast<int>(xvals.size()) != grid.nx * grid.ny) {
        std::cerr << "Data size (" << xvals.size() 
                  << ") does not match expected grid dimensions (" 
                  << grid.nx << " * " << grid.ny << ").\n";
        return false;
    }
    
    // Eigen::Map converts a raw array into an Eigen matrix.
    // Eigen stores matrices in column-major order by default.
    grid.x = Eigen::Map<const Eigen::MatrixXd>(xvals.data(), grid.nx, grid.ny);
    grid.y = Eigen::Map<const Eigen::MatrixXd>(yvals.data(), grid.nx, grid.ny);
    
    return true;
}

void plotGrid(const Grid &grid) {
    // C:\Users\jmmcl\OneDrive\Desktop\Project\C_Plus_Plus_Libraries\gnuplot\bin
    Gnuplot gp("C:\\Users\\jmmcl\\OneDrive\\Desktop\\Project\\C_Plus_Plus_Libraries\\gnuplot\\bin\\gnuplot.exe");

    // Set labels and title similar to MATLAB.
    gp << "set xlabel 'x-axis'\n";
    gp << "set ylabel 'y-axis'\n";
    gp << "set title '2D FVM Grid'\n";
    gp << "set grid\n";

    // --- Plot horizontal lines ---
    // For each row (fixed i) extract the corresponding data from x and y matrices.
    for (int i = 0; i < grid.nx; i++) {
        // Create a vector of pairs for the horizontal line.
        std::vector<std::pair<double, double>> line;
        line.reserve(grid.ny);
        for (int j = 0; j < grid.ny; j++) {
            line.emplace_back(grid.x(i, j), grid.y(i, j));
        }
        // Use "plot" for the first line then "replot" for subsequent ones.
        if (i == 0)
            gp << "plot '-' with lines lt rgb 'black' notitle\n";
        else
            gp << "replot '-' with lines lt rgb 'black' notitle\n";
        gp.send1d(line);
    }
    
    // --- Plot vertical lines ---
    // For each column (fixed j) extract the corresponding data from x and y matrices.
    for (int j = 0; j < grid.ny; j++) {
        std::vector<std::pair<double, double>> line;
        line.reserve(grid.nx);
        for (int i = 0; i < grid.nx; i++) {
            line.emplace_back(grid.x(i, j), grid.y(i, j));
        }
        gp << "replot '-' with lines lt rgb 'black' notitle\n";
        gp.send1d(line);
    }

    // Optionally pause the window until a key is pressed.
    std::cout << "Press Enter to exit..." << std::endl;
    std::cin.get();
}
