// GridHandler.cpp
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

    // Parse dimensions
    size_t pos = header.find("i=");
    if (pos == std::string::npos) {
        std::cerr << "Unable to find 'i=' in header.\n";
        return false;
    }
    pos += 2;
    size_t end = header.find(",", pos);
    nx = std::stoi(header.substr(pos, end - pos));

    pos = header.find("j=");
    if (pos == std::string::npos) {
        std::cerr << "Unable to find 'j=' in header.\n";
        return false;
    }
    pos += 2;
    end = header.find_first_of(" \t", pos);
    ny = std::stoi(header.substr(pos, end - pos));

    // Read x,y pairs
    std::vector<double> xvals, yvals;
    xvals.reserve(nx * ny);
    yvals.reserve(nx * ny);

    double xv, yv;
    char comma;
    while (file >> xv >> comma >> yv) {
        xvals.push_back(xv);
        yvals.push_back(yv);
    }
    if (static_cast<int>(xvals.size()) != nx * ny) {
        std::cerr << "Data size mismatch: got " << xvals.size()
                  << " but expected " << (nx * ny) << "\n";
        return false;
    }

    // Map into Eigen (column‐major)
    x = Eigen::Map<const Eigen::MatrixXd>(xvals.data(), nx, ny);
    y = Eigen::Map<const Eigen::MatrixXd>(yvals.data(), nx, ny);
    return true;
}

void GridHandler::haloCell() {
    Eigen::MatrixXd xAug = Eigen::MatrixXd::Zero(nx+2, ny+2);
    Eigen::MatrixXd yAug = Eigen::MatrixXd::Zero(nx+2, ny+2);

    xAug.block(1,1,nx,ny) = x;
    yAug.block(1,1,nx,ny) = y;

    // left/right ghosts (reflect)
    xAug.block(1, 0,    nx, 1) =  2*x.col(0)   - x.col(1);
    yAug.block(1, 0,    nx, 1) =  2*y.col(0)   - y.col(1);
    xAug.block(1, ny+1, nx, 1) =  2*x.col(ny-1)- x.col(ny-2);
    yAug.block(1, ny+1, nx, 1) =  2*y.col(ny-1)- y.col(ny-2);

    // top/bottom ghosts
    xAug.row(0)    = 2*xAug.row(1)   - xAug.row(2);
    yAug.row(0)    = 2*yAug.row(1)   - yAug.row(2);
    xAug.row(nx+1) = 2*xAug.row(nx)  - xAug.row(nx-1);
    yAug.row(nx+1) = 2*yAug.row(nx)  - yAug.row(nx-1);

    x = std::move(xAug);
    y = std::move(yAug);
    nx += 2;
    ny += 2;
}

void GridHandler::computeCellMetrics() {
    // 1) build ghosts
    haloCell();

    // 2) cell centers
    xCenter = x.block(0,0,nx-1,ny-1)
            + 0.5*(x.block(1,1,nx-1,ny-1) - x.block(0,0,nx-1,ny-1));
    yCenter = y.block(0,0,nx-1,ny-1)
            + 0.5*(y.block(1,1,nx-1,ny-1) - y.block(0,0,nx-1,ny-1));

    // 3) cell volumes (areas)
    {
        auto A = (x.block(1,1,nx-1,ny-1) - x.block(0,0,nx-1,ny-1)).array();
        auto B = (y.block(0,1,nx-1,ny-1) - y.block(1,0,nx-1,ny-1)).array();
        auto C = (y.block(1,1,nx-1,ny-1) - y.block(0,0,nx-1,ny-1)).array();
        auto D = (x.block(0,1,nx-1,ny-1) - x.block(1,0,nx-1,ny-1)).array();
        cellVolume = (0.5*(A*B - C*D)).matrix();
    }

    // 4) ξ‑face areas: interior faces only
    {
        // compute all xi-face segments (difference in j)
        int R = nx - 1;
        int C = ny - 1;
        auto dyXi = y.block(1,1,R,C) - y.block(1,0,R,C);  // size (R x C)
        auto dxXi = x.block(1,1,R,C) - x.block(1,0,R,C);
        // extract interior faces (exclude boundary ghosts)
        xArea_Xi = dyXi.block(1, 0, R-1, C-1);
        yArea_Xi = dxXi.block(1, 0, R-1, C-1);
    }

    
    // 5) η‑face areas: interior faces only
    {
        // compute all eta-face segments (difference in i)
        int R = nx - 1;
        int C = ny - 1;
        auto dyEt = y.block(1,1,R,C) - y.block(0,1,R,C);
        auto dxEt = x.block(1,1,R,C) - x.block(0,1,R,C);
        // extract interior faces
        xArea_Eta = dyEt.block(0, 1, R-1, C-1);
        yArea_Eta = dxEt.block(0, 1, R-1, C-1);
    }

    // 6) unit normals
    {
        auto xi_mag = (xArea_Xi.array().square() + yArea_Xi.array().square()).sqrt();
        xUnitNorm_Xi =  xArea_Xi.array() / xi_mag;
        yUnitNorm_Xi = -yArea_Xi.array() / xi_mag;

        // ── DEBUG CHECK: any NaNs in xi_mag? ──
        if ((xi_mag.array().isNaN()).any()) {
            std::cerr << "*** computeCellMetrics: NaN detected in xi_mag ***\n";
        } else {
            std::cerr << "computeCellMetrics: xi_mag OK (no NaNs)\n";
        }

        auto eta_mag = (xArea_Eta.array().square() + yArea_Eta.array().square()).sqrt();
        xUnitNorm_Eta =  xArea_Eta.array() / eta_mag;
        yUnitNorm_Eta = -yArea_Eta.array() / eta_mag;

        // ── DEBUG CHECK: any NaNs in eta_mag? ──
        if ((eta_mag.array().isNaN()).any()) {
            std::cerr << "*** computeCellMetrics: NaN detected in eta_mag ***\n";
        } else {
            std::cerr << "computeCellMetrics: eta_mag OK (no NaNs)\n";
        }
    }
}

void GridHandler::buildFaceMasks() {
    xiPlusMask  = Eigen::ArrayXXi::Ones(nx-1, ny);
    xiMinusMask = Eigen::ArrayXXi::Ones(nx-1, ny);
    etaPlusMask = Eigen::ArrayXXi::Ones(nx, ny-1);
    etaMinusMask= Eigen::ArrayXXi::Ones(nx, ny-1);

    for (int i = 0; i < nx-1; ++i) {
        for (int j = 0; j < ny; ++j) {
            if (i == 0 || i+1 == nx-1)    xiPlusMask(i,j)  = 0;
            if (i == 1 || i   == nx-2)    xiMinusMask(i,j) = 0;
        }
    }
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny-1; ++j) {
            if (j == 0 || j+1 == ny-1)    etaPlusMask(i,j)  = 0;
            if (j == 1 || j   == ny-2)    etaMinusMask(i,j) = 0;
        }
    }
}

void GridHandler::plotGrid(const std::string &plotTitle) {
    figure();
    std::vector<double> Xr, Yr, Xc, Yc;
    Xr.reserve((nx+1)*ny); Yr.reserve((nx+1)*ny);
    Xc.reserve((ny+1)*nx); Yc.reserve((ny+1)*nx);

    // rows
    for (int j=0; j<ny; ++j) {
        for (int i=0; i<nx; ++i) {
            Xr.push_back(x(i,j));
            Yr.push_back(y(i,j));
        }
        Xr.push_back(NAN); Yr.push_back(NAN);
    }
    plot(Xr,Yr,"k"); hold(on);

    // columns
    for (int i=0; i<nx; ++i) {
        for (int j=0; j<ny; ++j) {
            Xc.push_back(x(i,j));
            Yc.push_back(y(i,j));
        }
        Xc.push_back(NAN); Yc.push_back(NAN);
    }
    plot(Xc,Yc,"k");

    title(plotTitle);
    xlabel("x/L");
    ylabel("y/L");
}
