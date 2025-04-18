// src/Initialize.cpp

#include "Initialize.h"
#include <Eigen/Dense>
#include <cmath>

Initialize::Initialize(GridHandler &grid_,
                       State       &state_,
                       double      R_,
                       double      gamma_,
                       double      Cp_)
  : grid(grid_)
  , state(state_)
  , R(R_)
  , gamma(gamma_)
  , Cp(Cp_)
{ }

void Initialize::setInitialConditions(double P0, double T0, double M0) {
    // compute freestream from P0, T0, M0
    double a0   = std::sqrt(gamma * R * T0);
    double u0   = M0 * a0;
    double v0   = 0.0;

    int ni = grid.getNX();
    int nj = grid.getNY();

    // full-size primitives
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd v = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(ni, nj);

    // interior only
    P.block(1,1,ni-2,nj-2).setConstant(P0);
    u.block(1,1,ni-2,nj-2).setConstant(u0);
    v.block(1,1,ni-2,nj-2).setConstant(v0);
    T.block(1,1,ni-2,nj-2).setConstant(T0);

    // pack primitives → Q in State
    state.packToQ(P, u, v, T, R, gamma);
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
}

void Initialize::setInletConditions() {
    int nj = grid.getNY();
    auto &Q = state.getQ();

    // supersonic inlet: copy row 1 → row 0 ghosts
    for (int j = 0; j < nj; ++j) {
        Q[State::RHO]   (0,j) = Q[State::RHO]   (1,j);
        Q[State::RHO_U] (0,j) = Q[State::RHO_U] (1,j);
        Q[State::RHO_V] (0,j) = Q[State::RHO_V] (1,j);
        Q[State::ENERGY](0,j) = Q[State::ENERGY](1,j);
    }
}

void Initialize::setOutletConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto &Q = state.getQ();

    // supersonic outlet: zero‐gradient i=ni-1
    for (int j = 0; j < nj; ++j) {
        Q[State::RHO]   (ni-1,j) = Q[State::RHO]   (ni-2,j);
        Q[State::RHO_U] (ni-1,j) = Q[State::RHO_U] (ni-2,j);
        Q[State::RHO_V] (ni-1,j) = Q[State::RHO_V] (ni-2,j);
        Q[State::ENERGY](ni-1,j) = Q[State::ENERGY](ni-2,j);
    }
}

void Initialize::setWallConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto &Q = state.getQ();

    // inviscid slip wall on bottom (j=0) and top (j=nj-1)
    for (int i = 0; i < ni; ++i) {
        // bottom
        Q[State::RHO]   (i,0)   = Q[State::RHO]   (i,1);
        Q[State::RHO_U] (i,0)   = Q[State::RHO_U] (i,1);
        Q[State::RHO_V] (i,0)   = -Q[State::RHO_V] (i,1);
        Q[State::ENERGY](i,0)   = Q[State::ENERGY](i,1);

        // top
        Q[State::RHO]   (i,nj-1)= Q[State::RHO]   (i,nj-2);
        Q[State::RHO_U] (i,nj-1)= Q[State::RHO_U] (i,nj-2);
        Q[State::RHO_V] (i,nj-1)= -Q[State::RHO_V] (i,nj-2);
        Q[State::ENERGY](i,nj-1)= Q[State::ENERGY](i,nj-2);
    }
}
