#include "ComputationalGrid.h"
#include "Initialize.h"
#include <matplot/matplot.h>

using namespace matplot;


int main() {
    // Thermodynamic and Transport Properties
    const double R = 287.0; // [J / kg*K] gas constant
    const double Cp = 1005.0; // [J / kg*K] specific heat
    const double gamma = 1.400; // ratio of specific heats

    // Initial Conditions and Boundary Conditions \\
    // ------------------ Inlet ------------------ 
    double P_i = 11664; // [Pa] inlet pressure
    double T_i = 216.7; // [K] inlet temperature
    double M_i = 3.000; // inlet Mach Number

    ComputationalGrid grid;
    Initialize init(grid, R, gamma, Cp);

    // Read in grid file
    if (!grid.readGridFile("data\\g641x065uf.dat")) {
        return 1;
    }

    // Compute Cell Metrics
    grid.computeCellMetrics();

    // Apply free stream conditions
    init.setInitialConditions(P_i, T_i, M_i);





    /* Psuedo Code...

        temporal loop
        for ... defined using CFL number

            for .... (sweep entire mesh) construct fluxes

            end

            for .... (sweep entire mesh) difference fluxes

            end

            if .... (convergence criteria is met)
                break loop, YAY solver worked
            else
                continue
            end

        end
    
    */










    return 0;
}