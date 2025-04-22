#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "GridHandler.h"
#include "FluxState.h"

class Initialize {
public:
    Initialize(GridHandler &grid, FluxState &flux, double R, double gamma, double Cp);

    // sets the initial conditions based on free-stream quantities
    void setInitialConditions(double P0, double T0, double M0);

    // applies the inlet, outlet, and wall BCs
    void applyBoundaryConditions();

    // calculates the time step
    void computeTimeStep(double CFL);
    
    // getter method
    const double &getTimeStep() const { return dt; }

private:
    GridHandler &grid;
    FluxState &flux;
    double R, gamma, Cp;

    double dt; // timestep
    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();
};

#endif // INITIALIZE_H