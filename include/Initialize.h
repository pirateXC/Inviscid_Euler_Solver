// include/Initialize.h

#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "GridHandler.h"
#include "State.h"

class Initialize {
public:
    // grid_: mesh + halo
    // state_: owns Q & conversions
    // R_, gamma_, Cp_: gas properties
    Initialize(GridHandler &grid,
               State       &state,
               double      R,
               double      gamma,
               double      Cp);

    // fill interior cells with uniform P0, T0, M0
    void setInitialConditions(double P0,
                              double T0,
                              double M0);

    // impose BCs on the halo cells
    void applyBoundaryConditions();

private:
    GridHandler &grid;
    State       &state;
    double       R;
    double       gamma;
    double       Cp;

    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();
};

#endif // INITIALIZE_H
