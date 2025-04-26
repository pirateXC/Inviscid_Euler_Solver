// main.cpp

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <cmath>
#include <limits>

#include "GridHandler.h"
#include "FluxState.h"
#include "Initialize.h"
#include "StegerWarming.h"
#include "PostProcess.h"

// forward declarations of norm functions
double computeL2Norm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
);
double computeLinfNorm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
);

int main(){
    // redirect all output to debug.txt
    std::ofstream log("debug.txt");
    std::cout.rdbuf(log.rdbuf());
    std::cerr.rdbuf(log.rdbuf());

    std::cout << "+++ START main +++\n";

    // physical constants
    const double R     = 287.0,
                 Cp    = 1005.0,
                 gamma = 1.4,
                 CFL   = 0.5;
    const double P_inf = 11664.0,
                 T_inf = 216.7,
                 M_inf = 3.0;

    // --- GRID SETUP ---
    std::cout << "-> Reading grid file...\n";
    GridHandler grid;
    if(!grid.readGridFile("data/g641x065uf.dat")) return 1;

    std::cout << "-> Computing cell metrics...\n";
    grid.computeCellMetrics();

    std::cout << "-> Building face masks...\n";
    grid.buildFaceMasks();

    int ni = grid.getNX();
    int nj = grid.getNY();
    std::cout << "--- NX,NY = " << ni << "," << nj << " ---\n";

    // --- INITIALIZATION ---
    FluxState flux(ni, nj);
    Initialize init(grid, flux, R, gamma, Cp, P_inf, T_inf, M_inf);

    std::cout << "-> Applying initial conditions...\n";
    init.setInitialConditions();

    std::cout << "-> Applying boundary conditions...\n";
    init.applyBoundaryConditions();

    std::cout << "-> Packing primitives into Q...\n";
    flux.packToQ(R, gamma);

    std::cout << "-> Getting masks, volumes, normals...\n";
    // --- grab masks, volumes, normals once
    const auto& xiMask = grid.getXiMask();
    const auto& etaMask = grid.getEtaMask();
    const auto& V = grid.getCellVolume();    // (ni × nj)
    const auto& nxXi = grid.getXUnitNormXi();   // (ni-1 × nj)
    const auto& nyXi = grid.getYUnitNormXi();
    const auto& nxEt = grid.getXUnitNormEta();  // (ni × nj-1)
    const auto& nyEt = grid.getYUnitNormEta();

    std::cout
    << "CellVolume V: "<< V.rows() <<"×"<< V.cols() << "\n"
    << "xiMask: "  << xiMask.rows()  <<"×"<< xiMask.cols()  << "\n"
    << "etaMask: " << etaMask.rows() <<"×"<< etaMask.cols() << "\n";

    std::cout << "-> Establishing residuals...\n";
    // residual storage
    Eigen::MatrixXd res_rho    = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd res_rhou   = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd res_rhov   = Eigen::MatrixXd::Zero(ni, nj);
    Eigen::MatrixXd res_energy = Eigen::MatrixXd::Zero(ni, nj);

    double l2   = 1;
    double linf = 1;
    int iter    = 0;
    const int max_iter = 5000;

    std::cout << "+++ ENTER TIME MARCH LOOP +++\n";
    while (iter < max_iter && l2 > 1e-3) {
        std::cout << "***** Iter " << iter << " *****\n";

        // zero residuals
        res_rho.setZero();
        res_rhou.setZero();
        res_rhov.setZero();
        res_energy.setZero();

        // conserved variables reference
        auto& Q = flux.getQ();
        
        std::cout << "Starting grid sweep..\n";
        // assemble residual on interior cells
        for (int i = 1; i <= ni - 3; ++i) {
            for (int j = 1; j <= nj - 3; ++j) {
                //std::cout << "iter="<< iter << "  cell(i,j)=("<<i<<","<<j<<")\n";
                // west/east faces
                int maskW = xiMask(i-1, j);
                int maskE = xiMask(i, j);
                Eigen::Vector2d nW(nxXi(i-1, j-1), nyXi(i-1,j-1));
                Eigen::Vector2d nE(nxXi(i, j-1), nyXi(i, j-1));

                //std::cout << "Initializing west face left state..\n";
                Eigen::Vector4d QW_L(Q[0](i-1,j), Q[1](i-1,j), Q[2](i-1,j), Q[3](i-1,j));

                //std::cout << "Initializing west face right state..\n";
                Eigen::Vector4d QW_R(Q[0](i  ,j), Q[1](i  ,j), Q[2](i  ,j), Q[3](i  ,j));
                Eigen::Vector4d QE_L = QW_R;

                //std::cout << "Initializing east face right state..\n";
                Eigen::Vector4d QE_R(Q[0](i+1,j), Q[1](i+1,j), Q[2](i+1,j), Q[3](i+1,j));

                //std::cout << "Computing west face flux\n";
                auto FW = StegerWarming::computeFaceFlux(QW_L, QW_R, gamma, nW, maskW);
                
                //std::cout << "Computing east face flux\n";
                auto FE = StegerWarming::computeFaceFlux(QE_L, QE_R, gamma, nE, maskE);

                //std::cout << "Grabbing noth south masks\n";
                // south/north faces
                int maskS = etaMask(i, j-1);
                int maskN = etaMask (i, j);
                Eigen::Vector2d nS(nxEt(i-1, j-1), nyEt(i-1, j-1));
                Eigen::Vector2d nN(nxEt(i-1, j), nyEt(i-1, j));

                //std::cout << "Initializing south face left state..\n";
                Eigen::Vector4d QS_L(Q[0](i, j-1), Q[1](i, j-1), Q[2](i, j-1), Q[3](i, j-1));
                Eigen::Vector4d QS_R(Q[0](i, j), Q[1](i, j), Q[2](i, j), Q[3](i, j));
                Eigen::Vector4d QN_L = QS_R;

                //std::cout << "Initializing north face left state..\n";
                Eigen::Vector4d QN_R(Q[0](i,j+1), Q[1](i,j+1), Q[2](i,j+1), Q[3](i,j+1));

                //std::cout << "Computing south face flux\n";
                auto FS = StegerWarming::computeFaceFlux(QS_L, QS_R, gamma, nS, maskS);

                //std::cout << "Computing north face flux\n";
                auto FN = StegerWarming::computeFaceFlux(QN_L, QN_R, gamma, nN, maskN);

                // volume check
                double vol = V(i, j);
                if (!std::isfinite(vol) || vol <= 0.0) {
                    std::cerr << "*** Invalid volume at (" << i << "," << j << "): V=" << vol << "\n";
                    std::exit(EXIT_FAILURE);
                }
                double invV = 1.0 / vol;

                //std::cout << "Residual update..\n";
                // residual update
                res_rho   (i,j) = invV * ((FW(0) - FE(0)) + (FS(0) - FN(0)));
                res_rhou  (i,j) = invV * ((FW(1) - FE(1)) + (FS(1) - FN(1)));
                res_rhov  (i,j) = invV * ((FW(2) - FE(2)) + (FS(2) - FN(2)));
                res_energy(i,j) = invV * ((FW(3) - FE(3)) + (FS(3) - FN(3)));
            }
        }

        std::cout << "Computing norms...";
        // compute norms
        l2   = computeL2Norm  (res_rho, res_rhou, res_rhov, res_energy);
        linf = computeLinfNorm(res_rho, res_rhou, res_rhov, res_energy);

        std::cout << "Computing time step...";
        // compute dt and update
        init.computeTimeStep(CFL);
        double dt = init.getTimeStep();
        if (std::isnan(dt)) {
            std::cerr << "*** NaN in dt! ***\n";
            return 1;
        }

        Q[0].block(1,1,ni-2,nj-2) -= dt * res_rho  .block(1,1,ni-2,nj-2);
        Q[1].block(1,1,ni-2,nj-2) -= dt * res_rhou .block(1,1,ni-2,nj-2);
        Q[2].block(1,1,ni-2,nj-2) -= dt * res_rhov .block(1,1,ni-2,nj-2);
        Q[3].block(1,1,ni-2,nj-2) -= dt * res_energy.block(1,1,ni-2,nj-2);

        flux.unpackFromQ(R, gamma);
        init.applyBoundaryConditions();
        flux.packToQ(R, gamma);

        if (++iter % 100 == 0) {
            std::cout << "iter=" << iter
                      << "  L2=" << l2
                      << "  Linf=" << linf << "\n";
        }
    }

    // final plot
    auto P = flux.getPressure().block(1,1,ni-2,nj-2);
    plotFieldContour(grid, P, "Pressure Contours");

    std::cout << "+++ END main +++\n";
    return 0;
}

// definitions of the norm functions

double computeL2Norm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
) {
    int ni = r_rho.rows(), nj = r_rho.cols();
    int interior = (ni - 2) * (nj - 2);
    auto mag = (r_rho.array().square()
              + r_rhou.array().square()
              + r_rhov.array().square()
              + r_energy.array().square()).sqrt();

    // debug NaNs
    for (int i = 1; i <= ni - 2; ++i) {
        for (int j = 1; j <= nj - 2; ++j) {
            if (!std::isfinite(mag(i,j))) {
                std::cerr << "*** NaN detected in residual magnitude at ("
                          << i << "," << j << ")\n";
            }
        }
    }

    return std::sqrt(mag.block(1,1,ni-2,nj-2).array().square().sum() / interior);
}

double computeLinfNorm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
) {
    int ni = r_rho.rows(), nj = r_rho.cols();
    auto mag = (r_rho.array().square()
               + r_rhou.array().square()
               + r_rhov.array().square()
               + r_energy.array().square()).sqrt();
    return mag.block(1,1,ni-2,nj-2).maxCoeff();
}
