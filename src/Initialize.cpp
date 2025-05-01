#include "Initialize.h"
#include <cmath>
#include <iostream>

Initialize::Initialize(GridHandler &grid_,
                       FluxState &flux_,
                       double R_,
                       double gamma_,
                       double Cp_,
                    double P_inf_,
                    double T_inf_,
                    double M_inf_)
  : grid(grid_)
  , flux(flux_)
  , R(R_)
  , gamma(gamma_)
  , Cp(Cp_)
  , P_inf(P_inf_)
  , T_inf(T_inf_)
  , M_inf(M_inf_)
{}

void Initialize::setInitialConditions() {
    double a_inf = std::sqrt(gamma * R * T_inf);
    u_inf = M_inf * a_inf;
    v_inf = 0.0;

    int ni = grid.getNX();
    int nj = grid.getNY();

    flux.getPressure().block(1, 1, ni-2, nj-2).setConstant(P_inf);
    flux.getVelo_U().block(1, 1, ni-2, nj-2).setConstant(u_inf);
    flux.getVelo_V().block(1, 1, ni-2, nj-2).setConstant(v_inf);
    flux.getTemp().block(1, 1, ni-2, nj-2).setConstant(T_inf);
}

void Initialize::applyBoundaryConditions()
{
    // build Q everywhere from freestream primitives
    flux.packToQ(R, gamma);
    setInletConditions();
    setOutletConditions();
    copyGhostColsToQ();
    setWallConditions();
}

void Initialize::enforceBoundaryConditions()
{
    flux.unpackFromQ(R, gamma);   // primitives from new interior Q
    setInletConditions();         // update ghost col 0
    setOutletConditions();        // update ghost col ni_c-1
    copyGhostColsToQ();           // put those two cols back into Q
    setWallConditions();          // reflect wall, final unpackFromQ
}


void Initialize::setInletConditions() {
    const int nj_c = grid.getNY() - 1;
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    for(int j = 0; j < nj_c; ++j) {
        /*
        P(0,j) = P_inf;
        T(0,j) = T_inf;
        u(0,j) = u_inf;
        v(0,j) = v_inf;
        */
       
       P(0,j) = P(1,j);
       T(0,j) = T(1,j);
       u(0,j) = u(1,j);
       v(0,j) = v(1,j);
       
    }
}

void Initialize::setOutletConditions() {
    const int ni_c = grid.getNX() - 1;
    const int nj_c = grid.getNY() - 1;
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    for(int j = 0; j < nj_c; ++j) {
        P(ni_c-1,j) = P(ni_c-2,j);
        T(ni_c-1,j) = T(ni_c-2,j);
        u(ni_c-1,j) = u(ni_c-2,j);
        v(ni_c-1,j) = v(ni_c-2,j);
    }
}

void Initialize::setWallConditions() {
    auto& Q    = flux.getQ();
    const auto& nx_e = grid.getXUnitNormEta();
    const auto& ny_e = grid.getYUnitNormEta();

    int nEi = nx_e.rows();
    int nEj = nx_e.cols();
    int ni  = Q[0].rows();
    int nj  = Q[0].cols();

    for(int ie = 0; ie < nEi; ++ie) {
        int i   = ie + 1;
        int jG  = 0;
        int jI  = 1;
        double nx = nx_e(ie, 0), ny = ny_e(ie, 0);
        double rho  = Q[0](i,jI), rhou = Q[1](i,jI), rhov = Q[2](i,jI);
        double d    = rhou*nx + rhov*ny;
        Q[0](i,jG) = rho;
        Q[1](i,jG) = rhou - 2.0*d*nx;
        Q[2](i,jG) = rhov - 2.0*d*ny;
        Q[3](i,jG) = Q[3](i,jI);
    }

    int jG_top = nj - 1, jI_top = nj - 2, je = nEj - 1;
    for(int ie = 0; ie < nEi; ++ie) {
        int i = ie + 1;
        double nx = nx_e(ie, je), ny = ny_e(ie, je);
        double rho  = Q[0](i,jI_top), rhou = Q[1](i,jI_top), rhov = Q[2](i,jI_top);
        double d    = rhou*nx + rhov*ny;
        Q[0](i,jG_top) = rho;
        Q[1](i,jG_top) = rhou - 2.0*d*nx;
        Q[2](i,jG_top) = rhov - 2.0*d*ny;
        Q[3](i,jG_top) = Q[3](i,jI_top);
    }

    flux.unpackFromQ(R, gamma);
}


void Initialize::copyGhostColsToQ()
{
    const int ni_c = grid.getNX() - 1;   // last interior column + 1
    const int nj_c = grid.getNY() - 1;   // last interior row    + 1

    auto& Q = flux.getQ();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    auto copyCol = [&](int iCol)
    {
        /*  j runs 0 … nj_c  (≤ , not <) so the
            *bottom* and *top* ghost rows are covered         */
        for (int j = 0; j < nj_c; ++j) {
            const double rho  = P(iCol,j) / (R * T(iCol,j));
            const double rhou = rho * u(iCol,j);
            const double rhov = rho * v(iCol,j);
            const double rhoE = P(iCol,j)/(gamma-1.0)
                              + 0.5*rho*(u(iCol,j)*u(iCol,j) + v(iCol,j)*v(iCol,j));

            Q[0](iCol,j) = rho;
            Q[1](iCol,j) = rhou;
            Q[2](iCol,j) = rhov;
            Q[3](iCol,j) = rhoE;
        }
    };

    copyCol(0);          // inlet ghost column
    copyCol(ni_c-1);     // outlet ghost column
}

void Initialize::computeTimeStep(double CFL) {
    const auto& vol   = grid.getCellVolume();   // (ni-1)×(nj-1)
    const auto& Sx_i  = grid.getXAreaXi();      // (ni-1)×nj 
    const auto& Sy_i  = grid.getYAreaXi();
    const auto& Sx_e  = grid.getXAreaEta();     //  ni ×(nj-1)
    const auto& Sy_e  = grid.getYAreaEta();
    const auto& Umat  = flux.getVelo_U();       //  ni × nj
    const auto& Vmat  = flux.getVelo_V();
    const auto& Tmat  = flux.getTemp();

    int ni = vol.rows(), nj = vol.cols();
    double dt_min = std::numeric_limits<double>::infinity();

    for (int i = 0; i <= ni-3; ++i) {
        for (int j = 0; j <= nj-3; ++j) {
            int ii = i+1, jj = j+1;
            double Vi = vol(i, j);

            // xi‐face
            double sxW = Sx_i(i,j), syW = Sy_i(i,j);
            double AfW = std::hypot(sxW, syW);
            double unW = (sxW*Umat(ii-1,jj) + syW*Vmat(ii-1,jj)) / AfW;
            double cW  = std::sqrt(gamma * R * Tmat(ii,jj));
            double term_xi = (std::abs(unW) + cW) * AfW / Vi;

            // eta‐face
            double sxS = Sx_e(i,j), syS = Sy_e(i,j);
            double AfS = std::hypot(sxS, syS);
            double unS = (sxS*Umat(ii,jj-1) + syS*Vmat(ii,jj-1)) / AfS;
            double term_eta = (std::abs(unS) + cW) * AfS / Vi;

            // <<< change: sum, not max >>>
            double dt_local = CFL / (term_xi + term_eta);
            dt_min = std::min(dt_min, dt_local);
        }
    }

    dt = dt_min;
}

