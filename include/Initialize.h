#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "ComputationalGrid.h"

class Initialize {
public:
    Initialize(ComputationalGrid &grid,
               double R,
               double gamma,
               double Cp);

    // fill interior cells with uniform P₀, T₀, M₀
    void setInitialConditions(double P0,
                              double T0,
                              double M0);

    // impose BCs on the halo (ghost) cells
    void applyBoundaryConditions();

    // individual BC helpers (called by applyBoundaryConditions)
    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();

    // getters for your solver to read/write
    const Eigen::MatrixXd& getRho()   const { return  rho; }
    const Eigen::MatrixXd& getRhoU()  const { return  rhoU; }
    const Eigen::MatrixXd& getRhoV()  const { return  rhoV; }
    const Eigen::MatrixXd& getEnergy()const { return    E;  }

private:
    ComputationalGrid &grid;
    double R;         // gas constant
    double gamma;     // ratio of heats
    double Cp;        // [J / kg*K] specific heat

    // conserved variables, sized to grid with ghosts
    Eigen::MatrixXd rho, rhoU, rhoV, E;
};

#endif  // INITIALIZE_H
