#include "Initialize.h"
#include <Eigen/Dense>
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

    for (auto &mat : Q) {
        mat.setZero(ni, nj);
    }
}

void Initialize::setInitialConditions(double P0, double T0, double M0) {
    // compute freestream from P0, T0, M0
    double a0 = std::sqrt(gamma * R * T0);
    double u0 = M0 * a0;
    double v0 = 0.0;
    double rho0 = P0 / (R * T0);
    double e0 = P0 / (gamma - 1.0) + 0.5 * rho0 * (u0*u0 + v0*v0);

    int ni = grid.getNX();
    int nj = grid.getNY();

    // Allocate full-size fields, initialize to zero
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd v = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(ni, nj);

    // Fill only interior cells
    P.block(1, 1, ni-2, nj-2).setConstant(P0);
    u.block(1, 1, ni-2, nj-2).setConstant(u0);
    v.block(1, 1, ni-2, nj-2).setConstant(v0);
    T.block(1, 1, ni-2, nj-2).setConstant(T0);

    // populate state vector
    packToQ(P, u, v, T);
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
}

void Initialize::setInletConditions() {
    // supersonic inflow: prescribe all freestream -> copy into i=0 ghosts
    int nj = grid.getNY();

    for (int j = 0; j < nj; ++j) {
        Q[RHO](0,j) = Q[RHO](1,j); 
        Q[RHO_U](0,j) = Q[RHO_U](1,j);
        Q[RHO_V](0,j) = Q[RHO_V](1,j);
        Q[ENERGY](0,j) = Q[ENERGY](1,j);
    }
}

void Initialize::setOutletConditions() {
    // supersonic outflow: zero‐gradient -> copy last interior into ghost
    int ni = grid.getNX();
    int nj = grid.getNY();
    for (int j = 0; j < nj; ++j) {
        Q[RHO](ni-1,j) = Q[RHO](ni-2,j);
        Q[RHO_U](ni-1,j) = Q[RHO_U](ni-2,j);
        Q[RHO_V](ni-1,j) = Q[RHO_V](ni-2,j);
        Q[ENERGY](ni-1,j) = Q[ENERGY](ni-2,j);

    }
}

void Initialize::setWallConditions() {
    // inviscid slip wall on top/bottom: reflect normal momentum
    int ni = grid.getNX();
    int nj = grid.getNY();

    // bottom (j=0) and top (j=jmax-1)
    for (int i = 0; i < ni; ++i) {
        // bottom wall: j=0 -> j=1, flip normal velocity (rhoV)
        Q[RHO](i,0) = Q[RHO](i,1);
        Q[RHO_U](i,0) = Q[RHO_U](i,1);
        Q[RHO_V](i,0) = Q[RHO_V](i,1);
        Q[ENERGY](i,0) = Q[ENERGY](i,1);

        // top wall
        Q[RHO](i,nj-1) = Q[RHO](i,nj-2);
        Q[RHO_U](i,nj-1) = Q[RHO_U](i,nj-2);
        Q[RHO_V](i,nj-1) = -Q[RHO_V](i,nj-2);
        Q[ENERGY](i,nj-1) = Q[ENERGY](i,nj-2);
    }
    }

Eigen::MatrixXd Initialize::computePressure() const {
    const auto &rho = Q[RHO].array();
    const auto &rho_u = Q[RHO_U].array();
    const auto &rho_v = Q[RHO_V].array();
    const auto &energy = Q[ENERGY].array();

    return ((gamma - 1.0)*(rho * energy - 0.5*rho*( (rho_u/rho).square() + (rho_v/rho).square() ))).matrix();
}

Eigen::MatrixXd Initialize::computeTemp() const {

}

Eigen::MatrixXd Initialize::computeU_Velo() const {

}

Eigen::MatrixXd Initialize::computeV_Velo() const {
    
}
