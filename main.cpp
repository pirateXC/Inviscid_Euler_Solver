#include <matplot/matplot.h>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

using namespace matplot;

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
        std::cerr << "Failed to open file: " << filename << "\n";
        return false;
    }

    // Read header line (e.g., "ZONE i=641, j=65").
    std::string header;
    std::getline(file, header);

    // Parse grid dimensions from the header.
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

    // Reserve space in vectors.
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

    // Map the data into Eigen matrices (stored in column-major order by default).
    grid.x = Eigen::Map<const Eigen::MatrixXd>(xvals.data(), grid.nx, grid.ny);
    grid.y = Eigen::Map<const Eigen::MatrixXd>(yvals.data(), grid.nx, grid.ny);

    return true;
}

void plotGrid(const Grid &grid) {
    // Prepare a single vector for all row lines with NaN breaks.
    std::vector<double> X_rows, Y_rows;
    // Reserve estimated size: (grid.nx + 1) per row.
    X_rows.reserve((grid.nx + 1) * grid.ny);
    Y_rows.reserve((grid.nx + 1) * grid.ny);
    
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            X_rows.push_back(grid.x(i, j));
            Y_rows.push_back(grid.y(i, j));
        }
        // Insert NaN to break the line (requires <limits>).
        X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    // Plot all row lines at once in black.
    auto h_rows = plot(X_rows, Y_rows, "k");
    hold(on);  // Keep plotting subsequent lines on the same figure

    // Prepare a single vector for all column lines with NaN breaks.
    std::vector<double> X_cols, Y_cols;
    // Reserve estimated size: (grid.ny + 1) per column.
    X_cols.reserve((grid.ny + 1) * grid.nx);
    Y_cols.reserve((grid.ny + 1) * grid.nx);
    
    for (int i = 0; i < grid.nx; ++i) {
        for (int j = 0; j < grid.ny; ++j) {
            X_cols.push_back(grid.x(i, j));
            Y_cols.push_back(grid.y(i, j));
        }
        // Insert NaN to break the line.
        X_cols.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_cols.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    // Plot all column lines.
    auto h_cols = plot(X_cols, Y_cols, "k");

    // Set labels, title, and axis limits (preserving the title and black lines).
    title("Grid Mesh Plot");
    xlabel("x/L");
    ylabel("y/L");

    show();
}
