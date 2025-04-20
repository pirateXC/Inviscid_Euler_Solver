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

    flux.getPressure().block(1, 1, ni-2, nj-2).setConstant(P0);
    flux.getVelo_U().block(1, 1, ni-2, nj-2).setConstant(u0);
    flux.getVelo_V().block(1, 1, ni-2, nj-2).setConstant(v0);
    flux.getTemp().block(1, 1, ni-2, nj-2).setConstant(T0);

    applyBoundaryConditions();
    flux.packToQ(R, gamma);
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
}

void Initialize::setInletConditions() {
    int nj = grid.getNY();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& U = flux.getVelo_U();
    auto& V = flux.getVelo_V();

    for(int j=0; j< nj; ++j) {
        P(0,j) = P(1,j);
        T(0,j) = T(1,j);
        U(0,j) = U(1,j);
        V(0,j) = V(1,j);
    }
}

void Initialize::setOutletConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& U = flux.getVelo_U();
    auto& V = flux.getVelo_V();

    for(int j=0; j< nj; ++j) {
        P(ni-1,j) = P(ni-2,j);
        T(ni-1,j) = T(ni-2,j);
        U(ni-1,j) = U(ni-2,j);
        V(ni-1,j) = V(ni-2,j);
    }
}

void Initialize::setWallConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto& nx_e = grid.getXUnitNormEta();
    auto& ny_e = grid.getYUnitNormEta();
    auto& P    = flux.getPressure();
    auto& T    = flux.getTemp();
    auto& U    = flux.getVelo_U();
    auto& V    = flux.getVelo_V();

    // slip condition, inviscid flow
    for(int i=1; i<ni-1; ++i) {
        // bottom wall (j=0), interior at j=1
        {
            double nx = nx_e(i,0);
            double ny = ny_e(i,0);
            double u1 = U(i,1);
            double v1 = V(i,1);

            U(i,0) =  (1 - 2*nx*nx)*u1 - 2*nx*ny*v1;
            V(i,0) = -2*nx*ny*u1 + (1 - 2*ny*ny)*v1;
            P(i,0) = P(i,1);
            T(i,0) = T(i,1);
        }

        // top wall (j=nj-1), use normals at j=nj-2, interior at j=nj-2
        {
            double nx = nx_e(i,nj-2);
            double ny = ny_e(i,nj-2);
            double u1 = U(i,nj-2);
            double v1 = V(i, nj-2);

            U(i,nj-1) =  (1 - 2*nx*nx)*u1 - 2*nx*ny*v1;
            V(i,nj-1) = -2*nx*ny*u1 + (1 - 2*ny*ny)*v1;
            P(i,nj-1) = P(i,nj-2);
            T(i,nj-1) = T(i,nj-2);
        }
    }
}
