#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "GridHandler.h"

class Initialize {
public:
    Initialize(GridHandler &grid,
               double R,
               double gamma,
               double Cp);

    // fill interior cells with uniform P0, T0, M0
    void setInitialConditions(double P0,
                              double T0,
                              double M0);

    // impose BCs on the halo (ghost) cells
    void applyBoundaryConditions();

    // switch from primatives (P, u, v, T) to state vector Q
    void packToQ(
        const Eigen::MatrixXd &P,
        const Eigen::MatrixXd &u,
        const Eigen::MatrixXd &v,
        const Eigen::MatrixXd &T
    );

    // primative fields from Q
    Eigen::MatrixXd computePressure() const; // pressure
    Eigen::MatrixXd computeU_Velo() const; // u-component velocity
    Eigen::MatrixXd computeV_Velo() const; // v-component velocity
    Eigen::MatrixXd computeTemp() const; // temperature

    // getter methods
    const std::array<Eigen::MatrixXd,4>& getQ() const { return Q; }
    const Eigen::MatrixXd& getRho()    const { return Q[RHO]; }
    const Eigen::MatrixXd& getRhoU()   const { return Q[RHO_U]; }
    const Eigen::MatrixXd& getRhoV()   const { return Q[RHO_V]; }
    const Eigen::MatrixXd& getEnergy() const { return Q[ENERGY]; }

        // index enum for Q components
        enum Conserved { RHO = 0, RHO_U = 1, RHO_V = 2, ENERGY = 3 };

private:
    GridHandler &grid;
    double R;         // gas constant
    double gamma;     // ratio of heats
    double Cp;        // [J / kg*K] specific heat

    // individual BC helpers (called by applyBoundaryConditions)
    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();

    std::array<Eigen::MatrixXd, 4> Q; // state vector, [rho, rho*u, rho*v, rho*E]
};

#endif  // INITIALIZE_H
