#include "PostProcess.h"
#include <matplot/matplot.h>
#include <limits>

using namespace matplot;

// Filled contour plot with mesh overlay
void plotFieldContour(const GridHandler& grid,
                      const Eigen::MatrixXd& field,
                      const std::string& titleText)
{
    const auto& x = grid.getXCenter();
    const auto& y = grid.getYCenter();

    int ni = x.rows();
    int nj = x.cols();

    std::vector<std::vector<double>> X(ni, std::vector<double>(nj));
    std::vector<std::vector<double>> Y(ni, std::vector<double>(nj));
    std::vector<std::vector<double>> F(ni, std::vector<double>(nj));

    for (int i = 0; i < ni; ++i) {
        for (int j = 0; j < nj; ++j) {
            X[i][j] = x(i, j);
            Y[i][j] = y(i, j);
            F[i][j] = field(i, j);
        }
    }

    figure();
    contourf(X, Y, F, 30);
    colorbar();
    hold(on);

    // Add structured grid overlay (same as GridHandler::plotGrid)
    const auto& fullX = grid.getX();
    const auto& fullY = grid.getY();
    int nx = fullX.rows();
    int ny = fullX.cols();

    std::vector<double> X_rows, Y_rows, X_cols, Y_cols;
    X_rows.reserve((nx + 1) * ny);
    Y_rows.reserve((nx + 1) * ny);
    X_cols.reserve((ny + 1) * nx);
    Y_cols.reserve((ny + 1) * nx);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            X_rows.push_back(fullX(i, j));
            Y_rows.push_back(fullY(i, j));
        }
        X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X_cols.push_back(fullX(i, j));
            Y_cols.push_back(fullY(i, j));
        }
        X_cols.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_cols.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    plot(X_rows, Y_rows, "k");
    plot(X_cols, Y_cols, "k");

    title(titleText);
    xlabel("x/L");
    ylabel("y/L");
    axis(equal);
    show();
}

// Grid-only plotting function (structured lines)
void plotOnlyMesh(const GridHandler& grid,
                  const std::string& titleText)
{
    const auto& x = grid.getX();
    const auto& y = grid.getY();

    int nx = x.rows();
    int ny = x.cols();

    std::vector<double> X_rows, Y_rows, X_cols, Y_cols;
    X_rows.reserve((nx + 1) * ny);
    Y_rows.reserve((nx + 1) * ny);
    X_cols.reserve((ny + 1) * nx);
    Y_cols.reserve((ny + 1) * nx);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            X_rows.push_back(x(i, j));
            Y_rows.push_back(y(i, j));
        }
        X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            X_cols.push_back(x(i, j));
            Y_cols.push_back(y(i, j));
        }
        X_cols.push_back(std::numeric_limits<double>::quiet_NaN());
        Y_cols.push_back(std::numeric_limits<double>::quiet_NaN());
    }

    figure();
    plot(X_rows, Y_rows, "k");
    plot(X_cols, Y_cols, "k");
    title(titleText);
    xlabel("x/L");
    ylabel("y/L");
    axis(equal);
    show();
}
