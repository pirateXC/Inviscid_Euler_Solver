#include "gnuplot-iostream.h"
#include <vector>
#include <cmath>
#include <iostream>

int main() {
    //   "\"C:\\Users\\jmmcl\\OneDrive\\Desktop\\Project\\C_Plus_Plus_Libraries\\gnuplot\\bin\\gnuplot.exe\""
    Gnuplot gp("\"C:\\Users\\jmmcl\\OneDrive\\Desktop\\Project\\C_Plus_Plus_Libraries\\gnuplot\\bin\\gnuplot.exe\"");

    std::vector<std::pair<double, double>> xy_points;
    for (double x = 0; x < 10; x += 0.1) {
        xy_points.emplace_back(x, std::sin(x));
    }

    gp << "set title 'Sine Wave'\n";
    gp << "plot '-' with lines title 'sin(x)'\n";
    gp.send1d(xy_points);

    std::cout << "Press Enter to exit...\n";
    std::cin.get();
    return 0;
}
