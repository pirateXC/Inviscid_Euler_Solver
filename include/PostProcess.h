#ifndef POSTPROCESS_H
#define POSTPROCESS_H
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <string>
#include <Eigen/Dense>
#include "GridHandler.h"

// Plot a field as a filled contour with mesh overlay
void plotFieldContour(const GridHandler& grid,
                      const Eigen::MatrixXd& field,
                      const std::string& titleText);

// Plot the computational mesh alone (structured grid lines)
void plotOnlyMesh(const GridHandler& grid,
                  const std::string& titleText);

#endif // POSTPROCESS_H
