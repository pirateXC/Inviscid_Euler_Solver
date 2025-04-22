#include "GridHandler.h"
#include "FluxState.h"
#include "Initialize.h"
#include "PostProcess.h"

#include <matplot/matplot.h>
using namespace matplot;

int main() {
    // gas properties
    const double R     = 287.0;
    const double Cp    = 1005.0;
    const double gamma = 1.400;

    // freestream
    double P_i = 11664, T_i = 216.7, M_i = 3.0;

    GridHandler grid;
    if (!grid.readGridFile("data/g641x065uf.dat")) return 1;
    grid.computeCellMetrics();

    // now you know ni, nj
    int ni = grid.getNX();
    int nj = grid.getNY();

    // construct your State with the correct size
    FluxState flux(ni, nj);

    // only now can you build Initialize
    Initialize init(grid, flux, R, gamma, Cp);

    // apply ICs
    init.setInitialConditions(P_i, T_i, M_i);

    // apply BCs
    init.applyBoundaryConditions();

    // …rest of solver…

    return 0;
}
