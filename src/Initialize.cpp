#include "Initialize.h"
#include <cmath>

Initialize::Initialize(GridHandler &grid_,
                       FluxState &flux_,
                       double R_,
                       double gamma_,
                       double Cp_)
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
}

void Initialize::applyBoundaryConditions() {
    setInletConditions();
    setOutletConditions();
    setWallConditions();
    flux.packToQ(R, gamma);
}

void Initialize::setInletConditions() {
    int nj = grid.getNY();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    for(int j=0; j< nj; ++j) {
        P(0,j) = P(1,j);
        T(0,j) = T(1,j);
        u(0,j) = u(1,j);
        v(0,j) = v(1,j);
    }
}

void Initialize::setOutletConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    for(int j=0; j< nj; ++j) {
        P(ni-1,j) = P(ni-2,j);
        T(ni-1,j) = T(ni-2,j);
        u(ni-1,j) = u(ni-2,j);
        v(ni-1,j) = v(ni-2,j);
    }
}

void Initialize::setWallConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto& nx_e = grid.getXUnitNormEta();
    auto& ny_e = grid.getYUnitNormEta();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    // slip condition, inviscid flow
    for(int i=1; i<ni-1; ++i) {
        // bottom wall (j=0), interior at j=1
        {
            double nx = nx_e(i,0);
            double ny = ny_e(i,0);
            double u1 = u(i,1);
            double v1 = v(i,1);

            u(i,0) = (1 - 2*nx*nx)*u1 - 2*nx*ny*v1;
            v(i,0) = -2*nx*ny*u1 + (1 - 2*ny*ny)*v1;
            P(i,0) = P(i,1);
            T(i,0) = T(i,1);
        }

        // top wall (j=nj-1), use normals at j=nj-2, interior at j=nj-2
        {
            double nx = nx_e(i,nj-2);
            double ny = ny_e(i,nj-2);
            double u1 = u(i,nj-2);
            double v1 = v(i, nj-2);

            u(i,nj-1) = (1 - 2*nx*nx)*u1 - 2*nx*ny*v1;
            v(i,nj-1) = -2*nx*ny*u1 + (1 - 2*ny*ny)*v1;
            P(i,nj-1) = P(i,nj-2);
            T(i,nj-1) = T(i,nj-2);
        }
    }
}


void Initialize::computeTimeStep(double CFL) {
    int ni = grid.getNX();
    int nj = grid.getNY();
    const auto& u = flux.getVelo_U();        // size: (ni × nj)
    const auto& v = flux.getVelo_V();
    const auto& nx_xi = grid.getXUnitNormXi();   // size: (ni-2 × nj-3)
    const auto& ny_xi = grid.getYUnitNormXi();
    const auto& nx_eta = grid.getXUnitNormEta();  // size: (ni-3 × nj-2)
    const auto& ny_eta = grid.getYUnitNormEta();
    const auto& T = flux.getTemp();

    // convective velocities
    Eigen::ArrayXXd U_xi = 
       u.block(1,1, ni-2, nj-3).array() * nx_xi.array()
     + v.block(1,1, ni-2, nj-3).array() * ny_xi.array();

    Eigen::ArrayXXd U_et = 
       u.block(1,1, ni-3, nj-2).array() * nx_eta.array()
     + v.block(1,1, ni-3, nj-2).array() * ny_eta.array();
     
     // acoustic speed    
     Eigen::ArrayXXd c_xi = (gamma * R * T.block(1,1,ni-2,nj-3).array()).sqrt();
     Eigen::ArrayXXd c_eta = (gamma * R * T.block(1,1,ni-3,nj-2).array()).sqrt();

    // spectral radii
    Eigen::ArrayXXd spectralRadii_Xi  = U_xi.array().abs() + c_xi;
    Eigen::ArrayXXd spectralRadii_Eta = U_et.array().abs() + c_eta;

    // directional local dt
    Eigen::ArrayXXd dt_xi = spectralRadii_Xi.inverse();
    Eigen::ArrayXXd dt_eta = spectralRadii_Eta.inverse();

    // directional local dt min
    double dt_xi_min  = dt_xi .minCoeff();
    double dt_eta_min = dt_eta.minCoeff();

    // global dt
    dt = CFL * std::min(dt_xi_min, dt_eta_min);
}