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
        
        P(0,j) = P_inf;
        T(0,j) = T_inf;
        u(0,j) = u_inf;
        v(0,j) = v_inf;
        /*
       
       P(0,j) = P(1,j);
       T(0,j) = T(1,j);
       u(0,j) = u(1,j);
       v(0,j) = v(1,j);
      */ 
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
    auto& Q = flux.getQ();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    int nI = Q[0].rows();   // Q’s actual row count
    int nJ = Q[0].cols();   // Q’s actual col count

    auto copyCol = [&](int iCol)
    {
        for (int j = 0; j < nJ; ++j) {
            double rho  = P(iCol,j) / (R * T(iCol,j));
            double rhou = rho * u(iCol,j);
            double rhov = rho * v(iCol,j);
            double rhoE = P(iCol,j)/(gamma-1.0)
                        + 0.5*rho*(u(iCol,j)*u(iCol,j) + v(iCol,j)*v(iCol,j));

            Q[0](iCol,j) = rho;
            Q[1](iCol,j) = rhou;
            Q[2](iCol,j) = rhov;
            Q[3](iCol,j) = rhoE;
        }
    };

    copyCol(0);       // inlet ghost
    copyCol(nI-1);    // outlet ghost
}


void Initialize::computeTimeStep(double CFL) {
    // grab all the metrics  
    const auto& vol   = grid.getCellVolume();   // (ni-1)×(nj-1)
    const auto& Sx_i  = grid.getXAreaXi();      // (ni-1)×nj 
    const auto& Sy_i  = grid.getYAreaXi();      // (ni-1)×nj
    const auto& Sx_e  = grid.getXAreaEta();     // ni×(nj-1)
    const auto& Sy_e  = grid.getYAreaEta();     // ni×(nj-1)
    const auto& Umat  = flux.getVelo_U();       // ni×nj
    const auto& Vmat  = flux.getVelo_V();       // ni×nj
    const auto& Tmat  = flux.getTemp();         // ni×nj

    int ni = vol.rows();
    int nj = vol.cols();
    double dt_min = std::numeric_limits<double>::infinity();

    // Loop over cell-centers (skip halos)
    for (int i = 0; i < ni-1; ++i) {
        for (int j = 0; j < nj-1; ++j) {
            // map to flux indices
            int ic = i+1, jc = j+1;
            double Vcell = vol(i,j);

            // local speed of sound at cell center
            double c = std::sqrt(gamma * R * Tmat(ic,jc));

            // Xi-face contribution (west face of cell)
            double sx_i = Sx_i(i,j), sy_i = Sy_i(i,j);
            double A_i  = std::hypot(sx_i, sy_i);
            double U_n_i = (sx_i*Umat(ic,jc) + sy_i*Vmat(ic,jc)) / A_i;
            // exactly: term_xi = (|U_n|+c)*sqrt((Sx/V)^2+(Sy/V)^2)
            double term_xi = (std::abs(U_n_i) + c) * std::sqrt((sx_i/Vcell)*(sx_i/Vcell) + (sy_i/Vcell)*(sy_i/Vcell));

            // Eta-face contribution (south face of cell)
            double sx_e = Sx_e(ic,j-1), sy_e = Sy_e(ic,j-1);
            double A_e  = std::hypot(sx_e, sy_e);
            double U_n_e = (sx_e*Umat(ic,jc) + sy_e*Vmat(ic,jc)) / A_e;
            double term_eta = (std::abs(U_n_e) + c) * std::sqrt((sx_e/Vcell)*(sx_e/Vcell) + (sy_e/Vcell)*(sy_e/Vcell));

            double dt_local = CFL / std::max(term_xi, term_eta);

            dt_min = std::min(dt_min, dt_local);
        }
    }

    dt = dt_min;
}


