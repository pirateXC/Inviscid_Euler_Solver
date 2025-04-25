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

// forward declarations
void computeResidual(
    const GridHandler &grid,
    const FluxState   &flux,
    Eigen::MatrixXd   &res_rho,
    Eigen::MatrixXd   &res_rhou,
    Eigen::MatrixXd   &res_rhov,
    Eigen::MatrixXd   &res_energy,
    double            gamma
);
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
    // redirect all output to log file
    std::ofstream log("debug.txt");
    std::cout.rdbuf(log.rdbuf());
    std::cerr.rdbuf(log.rdbuf());

    std::cout << "+++ START main +++\n";
    const double R     = 287.0,
                 Cp    = 1005.0,
                 gamma = 1.4,
                 CFL   = 0.5;
    const double P_inf = 11664.0,
                 T_inf = 216.7,
                 M_inf = 3.0;

    std::cout << "-> Reading grid file...\n";
    GridHandler grid;
    if(!grid.readGridFile("data/g641x065uf.dat")) return 1;
    std::cout << "-> Computing cell metrics...\n";
    grid.computeCellMetrics();
    std::cout << "-> Building face masks...\n";
    grid.buildFaceMasks();
    int ni = grid.getNX(), nj = grid.getNY();
    std::cout << "--- NX,NY = "<<ni<<","<<nj<<" ---\n";

    std::cout << "-> Setting up FluxState and Initialize...\n";
    FluxState flux(ni,nj);
    Initialize init(grid, flux, R, gamma, Cp, P_inf, T_inf, M_inf);

    std::cout << "-> Applying initial conditions...\n";
    init.setInitialConditions();
    std::cout << "-> Applying boundary conditions...\n";
    init.applyBoundaryConditions();
    std::cout << "--- initial & BC done ---\n";

    std::cout << "-> Packing primitives into Q...\n";
    flux.packToQ(R, gamma);
    std::cout << "--- Q packed ---\n";

    std::cout << "-> Allocating residual arrays...\n";
    Eigen::MatrixXd res_rho(ni, nj), res_rhou(ni, nj), res_rhov(ni, nj), res_energy(ni, nj);

    double l2 = 1e9, linf = 1e9;
    int iter = 0;
    const int max_iter = 5000;

    std::cout << "+++ ENTER TIME MARCH LOOP +++\n";
    while(iter < max_iter && l2 > 1e-3) {
        std::cout << "***** Iter " << iter << " *****\n";

        std::cout << "-> Computing residual...\n";
        computeResidual(grid, flux, res_rho, res_rhou, res_rhov, res_energy, gamma);
        std::cout << "--- Residual computed ---\n";

        std::cout << "-> Computing L2 norm...\n";
        l2 = computeL2Norm(res_rho, res_rhou, res_rhov, res_energy);
        std::cout << "--- L2 norm = " << l2 << " ---\n";

        std::cout << "-> Computing Linf norm...\n";
        linf = computeLinfNorm(res_rho, res_rhou, res_rhov, res_energy);
        std::cout << "--- Linf norm = " << linf << " ---\n";

        std::cout << "-> Computing time step...\n";
        init.computeTimeStep(CFL);
        double dt = init.getTimeStep();
        std::cout << "--- dt = " << dt << " ---\n";
        if(std::isnan(dt)) {
            std::cerr << "*** NaN in dt! ***\n";
            return 1;
        }

        std::cout << "-> Updating conserved variables...\n";
        auto &Q = flux.getQ();
        Q[0].block(1,1,ni-2,nj-2) -= dt * res_rho  .block(1,1,ni-2,nj-2);
        Q[1].block(1,1,ni-2,nj-2) -= dt * res_rhou .block(1,1,ni-2,nj-2);
        Q[2].block(1,1,ni-2,nj-2) -= dt * res_rhov .block(1,1,ni-2,nj-2);
        Q[3].block(1,1,ni-2,nj-2) -= dt * res_energy.block(1,1,ni-2,nj-2);
        std::cout << "--- Q updated ---\n";

        std::cout << "-> Unpacking Q to primitives...\n";
        flux.unpackFromQ(R, gamma);
        std::cout << "-> Reapplying boundary conditions...\n";
        init.applyBoundaryConditions();
        std::cout << "-> Repacking Q...\n";
        flux.packToQ(R, gamma);
        std::cout << "--- Unpack/BC/pack done ---\n";

        if(iter % 100 == 0) {
            std::cout << "iter="<<iter<<"  L2="<<l2<<"  Linf="<<linf<<"\n";
        }
        ++iter;
    }

    std::cout << "-> Plotting final pressure contour...\n";
    auto P = flux.getPressure().block(1,1,ni-2,nj-2);
    plotFieldContour(grid, P, "Pressure Contours");

    std::cout << "+++ END main +++\n";
    return 0;
}

void computeResidual(
    const GridHandler &grid,
    const FluxState   &flux,
    Eigen::MatrixXd   &res_rho,
    Eigen::MatrixXd   &res_rhou,
    Eigen::MatrixXd   &res_rhov,
    Eigen::MatrixXd   &res_energy,
    double            gamma
) {
    const auto& xiMinusM = grid.getXiMinusMask();   // shape (nx-1 × ny)
    const auto& xiPlusM  = grid.getXiPlusMask();    // shape (nx-1 × ny)
    const auto& etaMinusM= grid.getEtaMinusMask();  // shape (nx × ny-1)
    const auto& etaPlusM = grid.getEtaPlusMask();   // shape (nx × ny-1)
    const auto& V    = grid.getCellVolume();
    const auto& nxXi = grid.getXUnitNormXi();
    const auto& nyXi = grid.getYUnitNormXi();
    const auto& nxEt = grid.getXUnitNormEta();
    const auto& nyEt = grid.getYUnitNormEta();
    const auto& Q    = flux.getQ();
    int ni = Q[0].rows(), nj = Q[0].cols();

    std::cout << "computeResidual: ni, nj = "<<ni<<","<<nj<<"\n";
    std::cout << "computeResidual: V size = "<<V.rows()<<"x"<<V.cols()<<"\n";
    std::cout << "computeResidual: nxXi size = "<<nxXi.rows()<<"x"<<nxXi.cols()<<"\n";
    std::cout << "computeResidual: nxEt size = "<<nxEt.rows()<<"x"<<nxEt.cols()<<"\n";
    std::cout << "computeResidual: starting loops\n";

    res_rho.setZero(ni,nj);
    res_rhou.setZero(ni,nj);
    res_rhov.setZero(ni,nj);
    res_energy.setZero(ni,nj);

    int eta_cols = nxEt.cols();
    for(int i=1; i<ni-2; ++i) {
        for(int j=1; j<nj-2; ++j) {
            int maskW = xiMinusM(i-1, j);    // west face between i-1↔i
            int maskE = xiPlusM (i,   j);

            Eigen::Vector2d nW{ nxXi(i-1,j-1), nyXi(i-1,j-1) };
            int ie = std::min(i, (int)nxXi.rows()-1);
            Eigen::Vector2d nE{ nxXi(ie,j-1), nyXi(ie,j-1) };

            Eigen::Vector4d QW_L{ Q[0](i-1,j), Q[1](i-1,j), Q[2](i-1,j), Q[3](i-1,j) };
            Eigen::Vector4d QW_R{ Q[0](i  ,j), Q[1](i  ,j), Q[2](i  ,j), Q[3](i  ,j) };
            Eigen::Vector4d QE_L=QW_R;
            Eigen::Vector4d QE_R{ Q[0](i+1,j), Q[1](i+1,j), Q[2](i+1,j), Q[3](i+1,j) };

            auto FW = StegerWarming::computeFaceFlux(QW_L, QW_R, gamma, nW, maskW);
            auto FE = StegerWarming::computeFaceFlux(QE_L, QE_R, gamma, nE, maskE);

            int jm = j-1;
            int maskS = etaMinusM(i, jm);
            int maskN = etaPlusM (i, j);
            int jn = std::min(j, eta_cols-1);
            Eigen::Vector2d nS{ nxEt(i-1,jm), nyEt(i-1,jm) };
            Eigen::Vector2d nN{ nxEt(i-1,jn), nyEt(i-1,jn) };

            Eigen::Vector4d QS_L{ Q[0](i,j-1), Q[1](i,j-1), Q[2](i,j-1), Q[3](i,j-1) };
            Eigen::Vector4d QS_R=QW_R;
            Eigen::Vector4d QN_L=QW_R;
            Eigen::Vector4d QN_R{ Q[0](i,j+1), Q[1](i,j+1), Q[2](i,j+1), Q[3](i,j+1) };

            auto FS = StegerWarming::computeFaceFlux(QS_L, QS_R, gamma, nS, maskS);
            auto FN = StegerWarming::computeFaceFlux(QN_L, QN_R, gamma, nN, maskN);

            double vol = V(i-1,j-1);
            if(!std::isfinite(vol) || vol <=0.0) {
                std::cerr<<"*** Invalid volume at ("<<i<<","<<j<<"): V="<<vol<<"\n";
                std::exit(EXIT_FAILURE);
            }
            double invV = 1.0/vol;
            res_rho(i,j)    = invV*((FW(0)-FE(0)) + (FS(0)-FN(0)));
            res_rhou(i,j)   = invV*((FW(1)-FE(1)) + (FS(1)-FN(1)));
            res_rhov(i,j)   = invV*((FW(2)-FE(2)) + (FS(2)-FN(2)));
            res_energy(i,j) = invV*((FW(3)-FE(3)) + (FS(3)-FN(3)));
        }
    }
    std::cout<<"computeResidual: loops complete\n";
}

double computeL2Norm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
) {
    int ni = r_rho.rows(), nj = r_rho.cols();
    int interior = (ni-3)*(nj-3);
    auto mag = (r_rho.array().square()
              + r_rhou.array().square()
              + r_rhov.array().square()
              + r_energy.array().square()).sqrt();

    for(int i=1;i<ni-2;++i){
        for(int j=1;j<nj-2;++j){
            if(!std::isfinite(mag(i,j))){
                std::cerr<<"*** NaN detected in residual magnitude at ("<<i<<","<<j<<") ***\n";
                std::cerr<<"r_rho="<<r_rho(i,j)
                          <<", r_rhou="<<r_rhou(i,j)
                          <<", r_rhov="<<r_rhov(i,j)
                          <<", r_energy="<<r_energy(i,j)<<"\n";
            }
        }
    }

    return std::sqrt(mag.block(1,1,ni-3,nj-3).square().sum()/interior);
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
    return mag.block(1,1,ni-3,nj-3).maxCoeff();
}
