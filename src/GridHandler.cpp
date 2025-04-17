#include "GridHandler.h"
#include <matplot/matplot.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <limits>

using namespace matplot;

GridHandler::GridHandler() : nx(0), ny(0) { }

bool GridHandler::readGridFile(const std::string &filename) {
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
        nx = std::stoi(header.substr(pos, end - pos));
    } else {
        std::cerr << "Unable to find 'i=' in header.\n";
        return false;
    }

    pos = header.find("j=");
    if (pos != std::string::npos) {
        pos += 2;
        size_t end = header.find_first_of(" \t", pos);
        ny = std::stoi(header.substr(pos, end - pos));
    } else {
        std::cerr << "Unable to find 'j=' in header.\n";
        return false;
    }

    // Reserve space in vectors.
    std::vector<double> xvals;
    std::vector<double> yvals;
    xvals.reserve(nx * ny);
    yvals.reserve(nx * ny);

    double xVal, yVal;
    char comma;
    while (file >> xVal >> comma >> yVal) {
        xvals.push_back(xVal);
        yvals.push_back(yVal);
    }

    if (static_cast<int>(xvals.size()) != nx * ny) {
        std::cerr << "Data size (" << xvals.size() 
                  << ") does not match expected grid dimensions (" 
                  << nx << " * " << ny << ").\n";
        return false;
    }

    // Map the data into Eigen matrices (data stored in column-major order).
    x = Eigen::Map<const Eigen::MatrixXd>(xvals.data(), nx, ny);
    y = Eigen::Map<const Eigen::MatrixXd>(yvals.data(), nx, ny);

    return true;
}

void GridHandler::haloCell() {
    // Create augmented matrices with 2 extra rows and 2 extra columns.
    Eigen::MatrixXd xAugGrid = Eigen::MatrixXd::Zero(nx + 2, ny + 2);
    Eigen::MatrixXd yAugGrid = Eigen::MatrixXd::Zero(nx + 2, ny + 2);

    // Copy the original grid into the center of the augmented grid.
    xAugGrid.block(1, 1, nx, ny) = x;
    yAugGrid.block(1, 1, nx, ny) = y;

    // Horizontal boundaries: left and right sides.
    for (int i = 0; i < nx; ++i) {
        // Left Boundary: reflect the first column.
        xAugGrid(i + 1, 0) = 2.0 * x(i, 0) - x(i, 1);
        yAugGrid(i + 1, 0) = 2.0 * y(i, 0) - y(i, 1);

        // Right Boundary: reflect the last column.
        // Note: valid original column indices are 0 to (ny-1).
        xAugGrid(i + 1, ny + 1) = 2.0 * x(i, ny - 1) - x(i, ny - 2);
        yAugGrid(i + 1, ny + 1) = 2.0 * y(i, ny - 1) - y(i, ny - 2);
    }

    // Vertical boundaries: top and bottom sides.
    for (int j = 0; j < ny + 2; ++j) {
        // Top Boundary: reflect the first interior rows.
        xAugGrid(0, j) = 2.0 * xAugGrid(1, j) - xAugGrid(2, j);
        yAugGrid(0, j) = 2.0 * yAugGrid(1, j) - yAugGrid(2, j);

        // Bottom Boundary: reflect the last interior rows.
        xAugGrid(nx + 1, j) = 2.0 * xAugGrid(nx, j) - xAugGrid(nx - 1, j);
        yAugGrid(nx + 1, j) = 2.0 * yAugGrid(nx, j) - yAugGrid(nx - 1, j);
    }

    // Update member variables with the augmented grid.
    x = xAugGrid;
    y = yAugGrid;
    nx += 2;
    ny += 2;
}

void GridHandler::computeCellMetrics() {
    // Extend the grid by adding ghost cells.
    haloCell();

    // Compute cell-centered coordinates for the grid
    xCenter = x.block(0, 0, nx - 1, ny - 1) +
              0.5 * (x.block(1, 1, nx - 1, ny - 1) - x.block(0, 0, nx - 1, ny - 1));
    yCenter = y.block(0, 0, nx - 1, ny - 1) +
              0.5 * (y.block(1, 1, nx - 1, ny - 1) - y.block(0, 0, nx - 1, ny - 1));

    // Compute cell volumes using the determinant method
    cellVolume = Eigen::MatrixXd::Zero(nx - 1, ny - 1);

    for (int i = 0; i < nx - 1; ++i) {
        for (int j = 0; j < ny - 1; ++j) {
            cellVolume(i, j) = 0.5 * ((x(i + 1, j + 1) - x(i, j)) * (y(i, j + 1) - y(i + 1, j)) -
                                      (y(i + 1, j + 1) - y(i, j)) * (x(i, j + 1) - x(i + 1, j)));
        }
    }

    // Compute face areas in the ξ (xi) direction
    xArea_Xi = Eigen::MatrixXd::Zero(nx - 2, ny - 3);
    yArea_Xi = Eigen::MatrixXd::Zero(nx - 2, ny - 3);

    for (int i = 1; i < nx - 1; ++i) {
        for (int j = 1; j < ny - 2; ++j) {
            xArea_Xi(i - 1, j - 1) = y(i + 1, j + 1) - y(i + 1, j);
            yArea_Xi(i - 1, j - 1) = x(i + 1, j + 1) - x(i + 1, j);
        }
    }

    // Compute face areas in the η (eta) direction
    xArea_Eta = Eigen::MatrixXd::Zero(nx - 3, ny - 1);
    yArea_Eta = Eigen::MatrixXd::Zero(nx - 3, ny - 1);

    for (int i = 1; i < nx - 2; ++i) {
        for (int j = 1; j < ny - 1; ++j) {
            xArea_Eta(i - 1, j - 1) = y(i + 1, j + 1) - y(i, j + 1);
            yArea_Eta(i - 1, j - 1) = x(i + 1, j + 1) - x(i, j + 1);
        }
    }
}


void GridHandler::plotGrid(const std::string &plotTitle) {
    figure();

    // Prepare vectors for row lines with NaN separators.
    std::vector<double> X_rows, Y_rows;
    X_rows.reserve((nx + 1) * ny);
    Y_rows.reserve((nx + 1) * ny);
    
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            X_rows.push_back(x(i, j));
            Y_rows.push_back(y(i, j));
        }
        // Insert NaN to break the line.
        X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    plot(X_rows, Y_rows, "k");
    hold(on);

    // Prepare vectors for column lines with NaN separators.
    std::vector<double> X_cols, Y_cols;
    X_cols.reserve((ny + 1) * nx);
    Y_cols.reserve((ny + 1) * nx);
    
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X_cols.push_back(x(i, j));
            Y_cols.push_back(y(i, j));
        }
        // Insert NaN to break the line.
        X_cols.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_cols.push_back(std::numeric_limits<double>::quiet_NaN());
    }
    plot(X_cols, Y_cols, "k");

    title(plotTitle);
    xlabel("x/L");
    ylabel("y/L");
}
