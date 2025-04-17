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

    // individual BC helpers (called by applyBoundaryConditions)
    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();

    // getter methods
    const Eigen::MatrixXd& getRho()   const { return  rho; }
    const Eigen::MatrixXd& getRhoU()  const { return  rhoU; }
    const Eigen::MatrixXd& getRhoV()  const { return  rhoV; }
    const Eigen::MatrixXd& getEnergy()const { return    E;  }

private:
    GridHandler &grid;
    double R;         // gas constant
    double gamma;     // ratio of heats
    double Cp;        // [J / kg*K] specific heat

    Eigen::MatrixXd rho, rhoU, rhoV, E; // conserved variables, sized to grid with ghosts
    Eigen::MatrixXd Q; // flux vector
};

#endif  // INITIALIZE_H
