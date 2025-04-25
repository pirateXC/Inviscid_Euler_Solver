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

    for(int j = 0; j < nj; ++j) {
        P(0,j) = P_inf;
        T(0,j) = T_inf;
        u(0,j) = u_inf;
        v(0,j) = v_inf;
    }
}

void Initialize::setOutletConditions() {
    int ni = grid.getNX();
    int nj = grid.getNY();
    auto& P = flux.getPressure();
    auto& T = flux.getTemp();
    auto& u = flux.getVelo_U();
    auto& v = flux.getVelo_V();

    for(int j = 0; j < nj; ++j) {
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

    for(int i = 1; i < ni - 3; ++i) {
        // bottom wall
        {
            double nx = nx_e(i,0);
            double ny = ny_e(i,0);
            double u1 = u(i,1);
            double v1 = v(i,1);
        
            u(i,0) = (ny*ny - nx*nx)*u1 - 2*nx*ny*v1;
            v(i,0) = (nx*nx - ny*ny)*v1 - 2*nx*ny*u1;
            P(i,0) = P(i,1);
            T(i,0) = T(i,1);
        }

        // top wall
        {
            double nx = nx_e(i,nj-3);
            double ny = ny_e(i,nj-3);
            double u1 = u(i,nj-2);
            double v1 = v(i, nj-2);

            u(i,nj-1) = (ny*ny - nx*nx)*u1 - 2*nx*ny*v1;
            v(i,nj-1) = (nx*nx - ny*ny)*v1 - 2*nx*ny*u1;
            P(i,nj-1) = P(i,nj-2);
            T(i,nj-1) = T(i,nj-2);
        }
    }
}

void Initialize::computeTimeStep(double CFL) {
    std::cout << "Op_Init_computeTimeStep: start\n";

    // Fetch grid & field sizes
    int ni = grid.getNX();
    int nj = grid.getNY();
    std::cout << "  ni=" << ni << "  nj=" << nj << "\n";

    // Grab references
    const auto& u = flux.getVelo_U();        // (ni × nj)
    const auto& v = flux.getVelo_V();        // (ni × nj)
    const auto& T = flux.getTemp();          // (ni × nj)
    const auto& nx_xi  = grid.getXUnitNormXi();   // (ni-2 × nj-3)
    const auto& ny_xi  = grid.getYUnitNormXi();   // (ni-2 × nj-3)
    const auto& nx_eta = grid.getXUnitNormEta();  // (ni-3 × nj-2)
    const auto& ny_eta = grid.getYUnitNormEta();  // (ni-3 × nj-2)

    // Print out all the shapes
    std::cout << "  u: " << u.rows() << "×" << u.cols()
              << "  nx_xi: " << nx_xi.rows() << "×" << nx_xi.cols() << "\n";
    std::cout << "  v: " << v.rows() << "×" << v.cols()
              << "  ny_xi: " << ny_xi.rows() << "×" << ny_xi.cols() << "\n";
    std::cout << "  T: " << T.rows() << "×" << T.cols() << "\n";
    std::cout << "  nx_eta: " << nx_eta.rows() << "×" << nx_eta.cols()
              << "  ny_eta: " << ny_eta.rows() << "×" << ny_eta.cols() << "\n";

    // Now use those shapes for the blocks
    int NXI = nx_xi.rows(), NYI = nx_xi.cols();
    std::cout << "  using xi‑block size: " << NXI << "×" << NYI << "\n";
    auto uXi  = u.block(1, 1, NXI, NYI).array();
    auto vXi  = v.block(1, 1, NXI, NYI).array();

    Eigen::ArrayXXd U_xi = uXi * nx_xi.array()
                        + vXi * ny_xi.array();
    std::cout << "  U_xi computed (" 
              << U_xi.rows() << "×" << U_xi.cols() << ")\n";

    int NXE = nx_eta.rows(), NYE = nx_eta.cols();
    std::cout << "  using eta‑block size: " << NXE << "×" << NYE << "\n";
    auto uEt  = u.block(1, 1, NXE, NYE).array();
    auto vEt  = v.block(1, 1, NXE, NYE).array();

    Eigen::ArrayXXd U_eta = uEt * nx_eta.array()
                         + vEt * ny_eta.array();
    std::cout << "  U_eta computed (" 
              << U_eta.rows() << "×" << U_eta.cols() << ")\n";

    auto T_xi  = T.block(1, 1, NXI, NYI).array();
    auto T_eta = T.block(1, 1, NXE, NYE).array();

    Eigen::ArrayXXd c_xi  = (gamma * R * T_xi).sqrt();
    Eigen::ArrayXXd c_eta = (gamma * R * T_eta).sqrt();
    std::cout << "  c_xi, c_eta computed\n";

    Eigen::ArrayXXd rho_xi  = (U_xi.abs() + c_xi).inverse();
    Eigen::ArrayXXd rho_eta = (U_eta.abs() + c_eta).inverse();
    double dt_xi_min  = rho_xi.minCoeff();
    double dt_eta_min = rho_eta.minCoeff();
    std::cout << "  dt_xi_min=" << dt_xi_min
              << "  dt_eta_min=" << dt_eta_min << "\n";

    dt = CFL * std::min(dt_xi_min, dt_eta_min);
    if(!std::isfinite(dt) || dt<=0){
        std::cerr<<"Invalid dt in computeTimeStep: "<<dt<<"\n";
        dt=1e-8;
    }
    std::cout << "Op_Init_computeTimeStep: dt=" << dt << "\n";
    std::cout << "Op_Init_computeTimeStep: end\n";
}

