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




    if (!grid.readGridFile("g641x065uf.dat")) {
        return 1;
    }

    // Plot the original grid.
    grid.plotGrid("Original Grid");

    // Augment the grid with halo cells and plot the result.
    grid.haloCell();
    grid.plotGrid("Augmented Grid (with Halo Cells)");

    // Display all windows.
    matplot::show();
    
    return 0;
}