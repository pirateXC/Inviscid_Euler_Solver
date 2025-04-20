#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "GridHandler.h"
#include "FluxState.h"

class Initialize {
public:
    Initialize(GridHandler &grid,
               FluxState       &flux,
               double           R,
               double           gamma,
               double           Cp);

    void setInitialConditions(double P0,
                              double T0,
                              double M0);

    void applyBoundaryConditions();

private:
    GridHandler &grid;
    FluxState   &flux;
    double       R, gamma, Cp;

    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();
};

#endif // INITIALIZE_H