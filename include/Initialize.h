#ifndef INITIALIZE_H
#define INITIALIZE_H

#include <Eigen/Dense>
#include "GridHandler.h"
#include "FluxState.h"

class Initialize {
public:
    Initialize(GridHandler &grid, FluxState &flux, double R, double gamma, double Cp, double P_inf, double T_inf, double M_f);

    // sets the initial conditions based on free-stream quantities
    void setInitialConditions();

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
    double P_inf, T_inf, M_inf;
    double u_inf, v_inf;

    double dt; // timestep
    void setInletConditions();
    void setOutletConditions();
    void setWallConditions();
};

#endif // INITIALIZE_H