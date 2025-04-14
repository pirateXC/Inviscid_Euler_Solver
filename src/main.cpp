#include "ComputationalGrid.h"
#include <matplot/matplot.h>

using namespace matplot;

int main() {
    ComputationalGrid grid;
    
    // Thermodynamic and Transport Properties
    double R = 287.0; // [J / kg*K] gas constant
    double Cp = 1005.0; // [J / kg*K] specific heat
    double gamma = 1.400; // ratio of specific heats

    // Initial Conditions and Boundary Conditions \\
    // ------------------ Inlet ------------------ 
    double P_i = 11664; // [Pa] inlet pressure
    double T_I = 216.7; // [K] inlet temperature
    double M_i = 3.000; // inlet Mach Number



    if (!grid.readGridFile("data\\g641x065uf.dat")) {
        return 1;
    }

    // Plot the original grid.
    grid.plotGrid("2D Mesh");

    // Compute Cell Metrics
    grid.computeCellGeometry();

    grid.plotGrid("2D Mesh with Halo Cells");

    // Display all windows.
    matplot::show();
    



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