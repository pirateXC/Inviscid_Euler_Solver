#include "ComputationalGrid.h"
#include <matplot/matplot.h>

using namespace matplot;

int main() {
    ComputationalGrid grid;
    
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