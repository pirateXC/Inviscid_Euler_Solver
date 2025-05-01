#define EIGEN_STACK_ALLOCATION_LIMIT 0

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

int main() {
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

    // INITIALIZATION
    FluxState flux(ni, nj);
    Initialize init(grid, flux, R, gamma, Cp, P_inf, T_inf, M_inf);

    std::cout << "-> Applying initial conditions...\n";
    init.setInitialConditions();

    std::cout << "-> Applying boundary conditions...\n";
    init.applyBoundaryConditions();

    std::cout << "-> Getting masks, volumes, normals...\n";
    // grab masks, volumes, normals, cell face lengths
    const auto& xiMask = grid.getXiMask();
    const auto& etaMask = grid.getEtaMask();
    const auto& V = grid.getCellVolume();    // (ni × nj)
    const auto& nxXi = grid.getXUnitNormXi();   // (ni-1 × nj)
    const auto& nyXi = grid.getYUnitNormXi();
    const auto& nxEt = grid.getXUnitNormEta();  // (ni × nj-1)
    const auto& nyEt = grid.getYUnitNormEta();
    const auto xiFaceLength = grid.getXiFaceLength();
    const auto etaFaceLength = grid.getEtaFaceLength();

    std::cout
    << "CellVolume V: "<< V.rows() <<"×"<< V.cols() << "\n"
    << "xiMask: "  << xiMask.rows()  <<"×"<< xiMask.cols()  << "\n"
    << "etaMask: " << etaMask.rows() <<"×"<< etaMask.cols() << "\n";

    std::cout << "-> Establishing residuals...\n";
    // residual storage
    Eigen::MatrixXd res_rho = Eigen::MatrixXd::Zero(ni-1, nj-1);
    Eigen::MatrixXd res_rhou = Eigen::MatrixXd::Zero(ni-1, nj-1);
    Eigen::MatrixXd res_rhov = Eigen::MatrixXd::Zero(ni-1, nj-1);
    Eigen::MatrixXd res_energy = Eigen::MatrixXd::Zero(ni-1, nj-1);

    double l2 = 1;
    double linf = 1;
    int iter = 0;
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
        // cell centered grid sweep
        for (int i = 1; i <= ni - 3; ++i) {
            for (int j = 1; j <= nj - 3; ++j) {

               //std::cout << "iter="<< iter << "  cell(i,j)=("<<i<<","<<j<<")\n";
                // reference indexes
                int im_xi  = i-1;
                int jm_xi  = j-1;
                int im_eta = i-1;
                int jm_eta = j-1;
                

                // xi-integer masks
                int maskW = xiMask(im_xi, jm_xi);
                int maskE = xiMask(im_xi+1, jm_xi);
            
                // eta-integer masks
                int maskS = etaMask(im_eta, jm_eta);
                int maskN = etaMask(im_eta, jm_eta+1);
                /*
                std::cout<<"maskW="<<maskW<<" maskE="<<maskE
                <<" maskS="<<maskS<<" maskN="<<maskN<<"\n";
                */
                // west/east faces
                Eigen::Vector2d nW(nxXi(im_xi, jm_xi), nyXi(im_xi, jm_xi));
                Eigen::Vector2d nE(nxXi(im_xi + 1, jm_xi), nyXi(im_xi + 1, jm_xi));

                // south/north faces
                Eigen::Vector2d nS(nxEt(im_eta, jm_eta), nyEt(im_eta, jm_eta));
                Eigen::Vector2d nN(nxEt(im_eta, jm_eta + 1), nyEt(im_eta, jm_eta + 1));

               /* std::cout<<" (" << i << ", " << j << ") normals: nW="<<nW.transpose()
                <<"  nE="<<nE.transpose()
                <<"  nS="<<nS.transpose()
                <<"  nN="<<nN.transpose()<<"\n";
                */
                // determine cell face length
                double Lw = xiFaceLength(im_xi, jm_xi); 
                double Le = xiFaceLength(im_xi+1, jm_xi);
                double Ls = etaFaceLength(im_eta, jm_eta); 
                double Ln =  etaFaceLength(im_eta, jm_eta+1); 

                // left/right west/east fluxes
                Eigen::Vector4d QW_L(Q[0](i-1,j), Q[1](i-1,j), Q[2](i-1,j), Q[3](i-1,j));
                Eigen::Vector4d QW_R(Q[0](i ,j), Q[1](i ,j), Q[2](i ,j), Q[3](i,j));
                Eigen::Vector4d QE_L = QW_R;
                Eigen::Vector4d QE_R(Q[0](i+1,j), Q[1](i+1,j), Q[2](i+1,j), Q[3](i+1,j));

                // left/right south/north fluxes
                Eigen::Vector4d QS_L(Q[0](i, j-1), Q[1](i, j-1), Q[2](i, j-1), Q[3](i, j-1));
                Eigen::Vector4d QS_R(Q[0](i, j), Q[1](i, j), Q[2](i, j), Q[3](i, j));
                Eigen::Vector4d QN_L = QS_R;
                Eigen::Vector4d QN_R(Q[0](i,j+1), Q[1](i,j+1), Q[2](i,j+1), Q[3](i,j+1));

                // compute west/east face fluxes
                auto FW = StegerWarming::computeFaceFlux(QW_L, QW_R, gamma, nW, maskW);
                auto FE = StegerWarming::computeFaceFlux(QE_L, QE_R, gamma, nE, maskE);

                // compute south/north face fluxes
                auto FS = StegerWarming::computeFaceFlux(QS_L, QS_R, gamma, nS, maskS);
                auto FN = StegerWarming::computeFaceFlux(QN_L, QN_R, gamma, nN, maskN);

                // scale fluxes
                FW *= Lw;
                FE *= Le;
                FS *= Ls;
                FN *= Ln;

                // volume check
                double vol = V(i, j);
                if (!std::isfinite(vol) || vol <= 0.0) {
                    std::cerr << "*** Invalid volume at (" << i << "," << j << "): V=" << vol << "\n";
                    std::exit(EXIT_FAILURE);
                }
                double invV = 1.0 / vol;

                //std::cout << "Residual update..\n";
                // residual update
                res_rho(i,j)    = invV * ( FE(0) - FW(0) + FN(0) - FS(0) );
                res_rhou(i,j)   = invV * ( FE(1) - FW(1) + FN(1) - FS(1) );
                res_rhov(i,j)   = invV * ( FE(2) - FW(2) + FN(2) - FS(2) );
                res_energy(i,j)= invV * ( FE(3) - FW(3) + FN(3) - FS(3) );

                if (i==145 && j==1) {
                    int maskW = xiMask(im_xi, jm_xi), maskE = xiMask(im_xi+1, jm_xi);
                    int maskS = etaMask(im_eta, jm_eta), maskN = etaMask(im_eta, jm_eta+1);
                    std::cout<<"\n--- debug at (" << i << ", " << j << ") ---\n";
                    std::cout<<" masks W,E,S,N = "
                             <<maskW<<","<<maskE<<","<<maskS<<","<<maskN<<"\n";
                    std::cout<<" FW="<<FW.transpose()<<"\n";
                    std::cout<<" FE="<<FE.transpose()<<"\n";
                    std::cout<<" FS="<<FS.transpose()<<"\n";
                    std::cout<<" FN="<<FN.transpose()<<"\n";
                }
            }
        }

        //std::cout << "Computing norms...\n";
        // compute norms
        l2   = computeL2Norm  (res_rho, res_rhou, res_rhov, res_energy);
        linf = computeLinfNorm(res_rho, res_rhou, res_rhov, res_energy);

        //std::cout << "Computing time step...\n";
        // compute dt and update
        init.computeTimeStep(CFL);
        double dt = init.getTimeStep();

        if (std::isnan(dt)) {
            std::cerr << "*** NaN in dt! ***\n";
            return 1;
        }

        Q[0].block(1,1, ni-3, nj-3) -= dt * res_rho.block(1,1, ni-3, nj-3);
        Q[1].block(1,1, ni-3, nj-3) -= dt * res_rhou.block(1,1, ni-3, nj-3);
        Q[2].block(1,1, ni-3, nj-3) -= dt * res_rhov.block(1,1, ni-3, nj-3);
        Q[3].block(1,1, ni-3, nj-3) -= dt * res_energy.block(1,1, ni-3, nj-3);

        init.enforceBoundaryConditions();

        for (int i = 1; i < ni-1; ++i) {
            for (int j = 1; j < nj-1; ++j) {
              double rho = Q[0](i,j);
              // compute pressure p = (gamma-1)*(rho*E - ρ·ek) or fetch from flux.getPressure()
              double p = flux.getPressure()(i,j);
              if (!std::isfinite(rho) || rho <= 0.0 || !std::isfinite(p) || p <= 0.0) {
                std::cerr<<"*** Bad state at ("<<i<<","<<j<<"): rho="<<rho<<" p="<<p<<"\n";
                std::abort();
              }
            }
          }


        /*
        if (++iter % 100 == 0) {
            std::cout << "iter=" << iter
                      << "  L2=" << l2
                      << "  Linf=" << linf << "\n";
        }
            */          
        std::cout << "iter="<< iter << "  L2="<< l2 << "  Linf="<< linf << "  dt = " << dt << '\n';
        iter++;
    }

    // final plot
    //auto P = flux.getPressure().block(1,1,ni-3,nj-3);
    //plotFieldContour(grid, P, "Pressure Contours");

    if (iter == max_iter) {
         std::cout << "maximum number of iterations reached\n";
    } else {
        std::cout << "convergence met at iteration: "  << iter << "\n";
    }

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
    // dimensions
    int ni = r_rho.rows();
    int nj = r_rho.cols();
    // interior block size
    int ni_int = ni - 2;
    int nj_int = nj - 2;

    // grab each field’s interior block
    auto A = r_rho  .block(1,1,ni_int,nj_int).array();
    auto B = r_rhou .block(1,1,ni_int,nj_int).array();
    auto C = r_rhov .block(1,1,ni_int,nj_int).array();
    auto D = r_energy.block(1,1,ni_int,nj_int).array();

    // sum of squares
    double S = A.square().sum()
             + B.square().sum()
             + C.square().sum()
             + D.square().sum();

    return std::sqrt(S);
}

double computeLinfNorm(
    const Eigen::MatrixXd &r_rho,
    const Eigen::MatrixXd &r_rhou,
    const Eigen::MatrixXd &r_rhov,
    const Eigen::MatrixXd &r_energy
) {
    int ni = r_rho.rows();
    int nj = r_rho.cols();
    int ni_int = ni - 2;
    int nj_int = nj - 2;

    // max abs of each field over the interior
    double m1 = r_rho  .block(1,1,ni_int,nj_int).cwiseAbs().maxCoeff();
    double m2 = r_rhou .block(1,1,ni_int,nj_int).cwiseAbs().maxCoeff();
    double m3 = r_rhov .block(1,1,ni_int,nj_int).cwiseAbs().maxCoeff();
    double m4 = r_energy.block(1,1,ni_int,nj_int).cwiseAbs().maxCoeff();

    // the overall L_inf
    return std::max( std::max(m1,m2), std::max(m3,m4) );
}