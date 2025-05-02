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
    flux.unpackFromQ(R, gamma); // primitives from new interior Q
    setInletConditions(); // update ghost col 0
    setOutletConditions(); // update ghost col ni_c-1
    copyGhostColsToQ(); // put those two cols back into Q
    setWallConditions(); // reflect wall, final unpackFromQ
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
    const auto& V = grid.getCellVolume(); // (ni−1)×(nj−1)
    const auto& nxXi = grid.getXUnitNormXi(); // (ni−2)×(nj−3)
    const auto& nyXi = grid.getYUnitNormXi(); // (ni−2)×(nj−3)
    const auto& xiLen = grid.getXiFaceLength(); // (ni−2)×(nj−3)
    const auto& nxEta = grid.getXUnitNormEta(); // (ni−3)×(nj−2)
    const auto& nyEta = grid.getYUnitNormEta(); // (ni−3)×(nj−2)
    const auto& etaLen = grid.getEtaFaceLength(); // (ni−3)×(nj−2)
    const auto& Umat = flux.getVelo_U(); // ni_full×nj
    const auto& Vmat = flux.getVelo_V(); // ni×nj
    const auto& Tmat = flux.getTemp(); // ni×nj

    int nc_i = grid.getNX() - 1; // = ni−1
    int nc_j = grid.getNY() - 1; // = nj−1

    double dt_min = std::numeric_limits<double>::infinity();

    for (int i = 0; i < nc_i - 2; ++i) {
        for (int j = 0; j < nc_j - 2; ++j) {
            int ic = i + 1, jc = j + 1;
            double Vcell = V(i, j);
            double c = std::sqrt(gamma * R * Tmat(ic, jc));

            // xi‐face
            double Un_xi = nxXi(i, j)*Umat(ic, jc) + nyXi(i, j)*Vmat(ic, jc);
            double term_xi = (std::abs(Un_xi) + c) * xiLen(i, j) / Vcell;

            // eta‐face
            double Un_eta = nxEta(i, j)*Umat(ic, jc) + nyEta(i, j)*Vmat(ic, jc);
            double term_eta = (std::abs(Un_eta) + c) * etaLen(i, j) / Vcell;

            double dt_local = CFL / std::max(term_xi, term_eta);
            dt_min = std::min(dt_min, dt_local);
        }
    }

    dt = dt_min;
}



