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

    // Calculates the cell center and volume for each cell
    void computeCellGeometry();

    // Creates a new figure window and plots the grid (using matplot).
    void plotGrid(const std::string &windowTitle);

private:
    int nx, ny; // total grid points for i and j componenet
    Eigen::MatrixXd x; // grid nodes in i-component
    Eigen::MatrixXd y; // grid nodes in j-component
    Eigen::MatrixXd xCenter; // cell centers for i-component
    Eigen::MatrixXd yCenter; // cell centers for j-component
    Eigen::MatrixXd cellVolume; // cell volumes
    Eigen::MatrixXd xArea_Eta; // face areas for i-component in eta direction
    Eigen::MatrixXd yArea_Eta; // face areas for j-component in eta direction
    Eigen::MatrixXd xArea_Xi;  // face areas for i-component in xi direction
    Eigen::MatrixXd yArea_Xi; // face areas for j-component in xi direction
};

#endif  // COMPUTATIONALGRID_H
