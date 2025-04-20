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
    // build augmented grids
    Eigen::MatrixXd xAug = Eigen::MatrixXd::Zero(nx+2, ny+2);
    Eigen::MatrixXd yAug = Eigen::MatrixXd::Zero(nx+2, ny+2);
    xAug.block(1,1,nx,ny) = x;
    yAug.block(1,1,nx,ny) = y;

    // left ghost‐column: rows 1..nx, col 0
    xAug.block(1, 0, nx, 1) = 2*x.col(0)   - x.col(1);
    yAug.block(1, 0, nx, 1) = 2*y.col(0)   - y.col(1);

    // right ghost‐column: rows 1..nx, col ny+1
    xAug.block(1, ny+1, nx, 1) = 2*x.col(ny-1) - x.col(ny-2);
    yAug.block(1, ny+1, nx, 1) = 2*y.col(ny-1) - y.col(ny-2);

    // top ghost‐row: row 0, cols 0..ny+1
    xAug.row(0)    = 2*xAug.row(1)   - xAug.row(2);
    yAug.row(0)    = 2*yAug.row(1)   - yAug.row(2);

    // bottom ghost‐row: row nx+1, cols 0..ny+1
    xAug.row(nx+1) = 2*xAug.row(nx)  - xAug.row(nx-1);
    yAug.row(nx+1) = 2*yAug.row(nx)  - yAug.row(nx-1);

    x = std::move(xAug);
    y = std::move(yAug);
    nx += 2;  ny += 2;
}

void GridHandler::computeCellMetrics() {
    // extend the grid by adding ghost cells
    haloCell();

    // compute cell-centered coordinates for the grid
    xCenter = x.block(0, 0, nx - 1, ny - 1) + 0.5 * (x.block(1, 1, nx - 1, ny - 1) - x.block(0, 0, nx - 1, ny - 1));
    yCenter = y.block(0, 0, nx - 1, ny - 1) + 0.5 * (y.block(1, 1, nx - 1, ny - 1) - y.block(0, 0, nx - 1, ny - 1));

    // compute cell volumes using the determinant method
    cellVolume = Eigen::MatrixXd::Zero(nx - 1, ny - 1);

    auto A = (x.block(1, 1, nx - 1, ny - 1) - x.block(0, 0, nx - 1, ny - 1)).array();
    auto B = (y.block(0, 1, nx - 1, ny - 1) - y.block(1, 0, nx - 1, ny - 1)).array();
    auto C = (y.block(1, 1, nx - 1, ny - 1) - y.block(0, 0, nx - 1, ny - 1)).array();
    auto D = (x.block(0, 1, nx - 1, ny - 1) - x.block(1, 0, nx - 1, ny - 1)).array();
    cellVolume = (0.5 * (A * B - C * D)).matrix();

    // compute face areas in the xi-direction
    xArea_Xi = Eigen::MatrixXd::Zero(nx - 2, ny - 3);
    yArea_Xi = Eigen::MatrixXd::Zero(nx - 2, ny - 3);

    xArea_Xi =  y.block(2, 2, nx - 2, ny - 3)
              - y.block(2, 1, nx - 2, ny - 3);
    yArea_Xi =  x.block(2, 2, nx - 2, ny - 3)
              - x.block(2, 1, nx - 2, ny - 3);

    // compute face areas in the eta-direction
    xArea_Eta = Eigen::MatrixXd::Zero(nx - 3, ny - 1);
    yArea_Eta = Eigen::MatrixXd::Zero(nx - 3, ny - 1);

    xArea_Eta = y.block(2, 2, nx - 3, ny - 1) - y.block(1, 2, nx - 3, ny - 1);
    yArea_Eta = x.block(2, 2, nx - 3, ny - 1) - x.block(1, 2, nx - 3, ny - 1);

    // compute unit normals in xi-direction
    xUnitNorm_Xi =  xArea_Xi.array() / ((xArea_Xi.array().square() + yArea_Xi.array().square()).sqrt());
    yUnitNorm_Xi = -yArea_Xi.array() / ((xArea_Xi.array().square() + yArea_Xi.array().square()).sqrt());

      // compute unit normals in eta-direction
    xUnitNorm_Eta =  xArea_Eta.array() / ((xArea_Eta.array().square() + yArea_Eta.array().square() ).sqrt());
    yUnitNorm_Eta = -yArea_Eta.array() / ((xArea_Eta.array().square() + yArea_Eta.array().square() ).sqrt());
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
