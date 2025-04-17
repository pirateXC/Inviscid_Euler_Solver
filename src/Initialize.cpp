#include "Initialize.h"
#include <cmath>
#include <iostream>

Initialize::Initialize(GridHandler &grid_,
                       double R_,
                       double gamma_,
                       double Cp_)
  : grid(grid_), R(R_), gamma(gamma_), Cp(Cp_)
{
    // after grid.computeCellMetrics(), grid has been halo‑extended,
    const int &ni = grid.getNX();
    const int &nj = grid.getNY();

    rho  = Eigen::MatrixXd::Zero(ni, nj);
    rhoU = Eigen::MatrixXd::Zero(ni, nj);
    rhoV = Eigen::MatrixXd::Zero(ni, nj);
    E    = Eigen::MatrixXd::Zero(ni, nj);
}

void Initialize::setInitialConditions(double P0, double T0, double M0) {
    // compute freestream from P0, T0, M0
    double a0 = std::sqrt(gamma * R * T0);
    double u0 = M0 * a0;
    double v0 = 0.0;
    double rho0 = P0 / (R * T0);
    double e0 = P0 / (gamma - 1.0) + 0.5 * rho0 * (u0*u0 + v0*v0);

    // fill interior (skip ghost‑layer indices!)
    int ng = 1;                    // if you used 1‑cell halo
    int imax = rho.rows()  - ng;
    int jmax = rho.cols()  - ng;
    for (int i = ng; i < imax; ++i) {
        for (int j = ng; j < jmax; ++j) {
            rho(i,j)  = rho0;
            rhoU(i,j) = rho0 * u0;
            rhoV(i,j) = rho0 * v0;
            E(i,j)    = e0;
        }
    }
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
}

void Initialize::setInletConditions() {
    // supersonic inflow: prescribe all freestream -> copy into i=0 ghosts
    int jmax = rho.cols();
    for (int j = 0; j < jmax; ++j) {
        rho(0,j)  = rho(1,j);
        rhoU(0,j) = rhoU(1,j);
        rhoV(0,j) = rhoV(1,j);
        E(0,j)    = E(1,j);
    }
}

void Initialize::setOutletConditions() {
    // supersonic outflow: zero‐gradient -> copy last interior into ghost
    int imax = rho.rows();
    int jmax = rho.cols();
    for (int j = 0; j < jmax; ++j) {
        rho(imax-1,j)  = rho(imax-2,j);
        rhoU(imax-1,j) = rhoU(imax-2,j);
        rhoV(imax-1,j) = rhoV(imax-2,j);
        E(imax-1,j)    = E(imax-2,j);
    }
}

void Initialize::setWallConditions() {
    // inviscid slip wall on top/bottom: reflect normal momentum
    int imax = rho.rows();
    int jmax = rho.cols();

    // bottom (j=0) and top (j=jmax-1)
    for (int i = 0; i < imax; ++i) {
        // bottom wall: j=0 -> j=1, flip normal velocity (rhoV)
        rho(i,0)   = rho(i,1);
        rhoU(i,0)  = rhoU(i,1);
        rhoV(i,0)  = -rhoV(i,1);
        E(i,0)     = E(i,1);

        // top wall
        rho(i,jmax-1)  = rho(i,jmax-2);
        rhoU(i,jmax-1) = rhoU(i,jmax-2);
        rhoV(i,jmax-1) = -rhoV(i,jmax-2);
        E(i,jmax-1)    = E(i,jmax-2);
    }
}
