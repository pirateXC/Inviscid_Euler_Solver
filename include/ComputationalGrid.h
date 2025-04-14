#ifndef COMPUTATIONALGRID_H
#define COMPUTATIONALGRID_H

#include <Eigen/Dense>
#include <string>

class ComputationalGrid {
public:
    ComputationalGrid();
    
    // Reads the grid file, parses header information, and maps the data into matrices.
    bool readGridFile(const std::string &filename);
    
    // Augments the grid by adding halo cells on all four boundaries.
    void haloCell();
    
    // Creates a new figure window and plots the grid (using matplot).
    void plotGrid(const std::string &windowTitle);

private:
    int nx, ny;
    Eigen::MatrixXd x;
    Eigen::MatrixXd y;
};

#endif  // COMPUTATIONALGRID_H
