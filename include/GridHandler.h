#ifndef GRIDHANDLER_H
#define GRIDHANDLER_H

#include <Eigen/Dense>
#include <string>

class GridHandler {
public:
    GridHandler();
    
    // Reads the grid file, parses header information, and maps the data into matrices.
    bool readGridFile(const std::string &filename);
    
    // Augments the grid by adding halo cells on all four boundaries.
    void haloCell();

    // Calculates the cell center and volume for each cell and calculates the face areas in the xi and eta directions.
    void computeCellMetrics();

    // Creates a new figure window and plots the grid (using matplot).
    void plotGrid(const std::string &windowTitle);

    // getter methods
    const int &getNX() const { return nx; }
    const int &getNY() const { return ny; }
    const Eigen::MatrixXd &getX() const { return x; }
    const Eigen::MatrixXd &getY() const { return y; }
    const Eigen::MatrixXd &getXCenter() const { return xCenter; }
    const Eigen::MatrixXd &getYCenter() const { return yCenter; }
    const Eigen::MatrixXd &getCellVolume() const { return cellVolume; }
    const Eigen::MatrixXd &getXAreaXi() const { return xArea_Xi; }
    const Eigen::MatrixXd &getYAreaXi() const { return yArea_Xi; }
    const Eigen::MatrixXd &getXAreaEta() const { return xArea_Eta; }
    const Eigen::MatrixXd &getYAreaEta() const { return yArea_Eta; }

private:
    int nx, ny; // total grid points for i and j component
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
