#include <matplot/matplot.h>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <limits>

using namespace matplot;

// Structure definition for grid data.
struct Grid {
    int nx, ny; // total number x and y components
    Eigen::MatrixXd x;   // x-coordinate grid
    Eigen::MatrixXd y;   // y-coordinate grid
};

bool readGridFile(const std::string &filename, Grid &grid);
void haloCell(Grid& grid);
void plotGrid(const Grid &grid, const std::string &plotTitle);

int main() {
    Grid grid;
    if (!readGridFile("g641x065uf.dat", grid)) {
        std::cerr << "Error reading grid file.\n";
        return 1;
    }

    // Plot original grid
    plotGrid(grid, "2D-Mesh");

    // Augment the grid with halo cells
    haloCell(grid);

    // Plot the augmented grid (with halo cells)
    plotGrid(grid, "2D-Mesh with Halo Cells");

    return 0;
}

bool readGridFile(const std::string &filename, Grid &grid) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return false;
    }

    // Read header line (example: "ZONE i=641, j=65").
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

void haloCell(Grid& grid) {
    // Create augmented matrices with 2 extra rows and 2 extra columns.
    Eigen::MatrixXd xAugGrid = Eigen::MatrixXd::Zero(grid.nx + 2, grid.ny + 2);
    Eigen::MatrixXd yAugGrid = Eigen::MatrixXd::Zero(grid.nx + 2, grid.ny + 2);

    // Copy the original grid into the center of the augmented grid.
    xAugGrid.block(1, 1, grid.nx, grid.ny) = grid.x;
    yAugGrid.block(1, 1, grid.nx, grid.ny) = grid.y;

    // Horizontal boundaries: left and right sides.
    for (int i = 0; i < grid.nx; ++i) {
        // Left Boundary: reflect the first column.
        xAugGrid(i + 1, 0) = 2.0 * grid.x(i, 0) - grid.x(i, 1);
        yAugGrid(i + 1, 0) = 2.0 * grid.y(i, 0) - grid.y(i, 1);

        // Right Boundary: reflect the last column.
        // Note: grid.x has grid.ny columns (last index is grid.ny-1).
        xAugGrid(i + 1, grid.ny + 1) = 2.0 * grid.x(i, grid.ny - 1) - grid.x(i, grid.ny - 2);
        yAugGrid(i + 1, grid.ny + 1) = 2.0 * grid.y(i, grid.ny - 1) - grid.y(i, grid.ny - 2);
    }

    // Vertical boundaries: top and bottom sides.
    // The augmented grid has grid.nx+2 rows (indices 0 to grid.nx+1).
    for (int j = 0; j < grid.ny + 2; ++j) {
        // Top Boundary: reflect the first interior rows.
        xAugGrid(0, j) = 2.0 * xAugGrid(1, j) - xAugGrid(2, j);
        yAugGrid(0, j) = 2.0 * yAugGrid(1, j) - yAugGrid(2, j);

        // Bottom Boundary: reflect the last interior rows.
        xAugGrid(grid.nx + 1, j) = 2.0 * xAugGrid(grid.nx, j) - xAugGrid(grid.nx - 1, j);
        yAugGrid(grid.nx + 1, j) = 2.0 * yAugGrid(grid.nx, j) - yAugGrid(grid.nx - 1, j);
    }

    // Update the grid structure with the augmented matrices.
    grid.x = xAugGrid;
    grid.y = yAugGrid;
    grid.nx += 2;
    grid.ny += 2;
}

void plotGrid(const Grid &grid, const std::string &plotTitle) {
    figure();

    // Prepare vectors for row lines with NaN separators.
    std::vector<double> X_rows, Y_rows;
    X_rows.reserve((grid.nx + 1) * grid.ny);
    Y_rows.reserve((grid.nx + 1) * grid.ny);
    
    for (int j = 0; j < grid.ny; ++j) {
        for (int i = 0; i < grid.nx; ++i) {
            X_rows.push_back(grid.x(i, j));
            Y_rows.push_back(grid.y(i, j));
        }
        // Insert NaN to break the line.
        X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    plot(X_rows, Y_rows, "k");
    hold(on);

    // Prepare vectors for column lines with NaN separators.
    std::vector<double> X_cols, Y_cols;
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
    plot(X_cols, Y_cols, "k");

    title(plotTitle);
    xlabel("x/L");
    ylabel("y/L");

    show();
}