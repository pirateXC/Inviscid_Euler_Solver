#include "Initialize.h"
#include <cmath>

Initialize::Initialize(GridHandler &grid_,
                       FluxState       &flux_,
                       double           R_,
                       double           gamma_,
                       double           Cp_)
  : grid(grid_)
  , flux(flux_)
  , R(R_)
  , gamma(gamma_)
  , Cp(Cp_)
{ }

void Initialize::setInitialConditions(double P0, double T0, double M0) {
    double a0 = std::sqrt(gamma * R * T0);
    double u0 = M0 * a0;
    double v0 = 0.0;

    int ni = grid.getNX();
    int nj = grid.getNY();

    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd v = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(ni, nj);

    P.block(1, 1, ni-2, nj-2).setConstant(P0);
    u.block(1, 1, ni-2, nj-2).setConstant(u0);
    v.block(1, 1, ni-2, nj-2).setConstant(v0);
    T.block(1, 1, ni-2, nj-2).setConstant(T0);

    flux.packToQ(P, u, v, T, R, gamma);
    applyBoundaryConditions();
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
}

void Initialize::setInletConditions() {
    int nj = grid.getNY();
    auto &Q = flux.getQ();

    for (int j = 0; j < nj; ++j) {
        Q[FluxState::RHO]   (0, j) = Q[FluxState::RHO]   (1, j);
        Q[FluxState::RHO_U] (0, j) = Q[FluxState::RHO_U] (1, j);
        Q[FluxState::RHO_V] (0, j) = Q[FluxState::RHO_V] (1, j);
        Q[FluxState::ENERGY](0, j) = Q[FluxState::ENERGY](1, j);
    }
}

void Initialize::setOutletConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto &Q = flux.getQ();

    for (int j = 0; j < nj; ++j) {
        Q[FluxState::RHO]   (ni-1, j) = Q[FluxState::RHO]   (ni-2, j);
        Q[FluxState::RHO_U] (ni-1, j) = Q[FluxState::RHO_U] (ni-2, j);
        Q[FluxState::RHO_V] (ni-1, j) = Q[FluxState::RHO_V] (ni-2, j);
        Q[FluxState::ENERGY](ni-1, j) = Q[FluxState::ENERGY](ni-2, j);
    }
}

void Initialize::setWallConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto &Q = flux.getQ();

    for (int i = 0; i < ni; ++i) {
        // top wall
        Q[FluxState::RHO]   (i, 0)    = Q[FluxState::RHO]   (i, 1);
        Q[FluxState::RHO_U] (i, 0)    = Q[FluxState::RHO_U] (i, 1);
        Q[FluxState::RHO_V] (i, 0)    = -Q[FluxState::RHO_V] (i, 1);
        Q[FluxState::ENERGY](i, 0)    = Q[FluxState::ENERGY](i, 1);

        // bottom wall
        Q[FluxState::RHO]   (i, nj-1) = Q[FluxState::RHO]   (i, nj-2);
        Q[FluxState::RHO_U] (i, nj-1) = Q[FluxState::RHO_U] (i, nj-2);
        Q[FluxState::RHO_V] (i, nj-1) = -Q[FluxState::RHO_V] (i, nj-2);
        Q[FluxState::ENERGY](i, nj-1)= Q[FluxState::ENERGY](i, nj-2);
    }
}
