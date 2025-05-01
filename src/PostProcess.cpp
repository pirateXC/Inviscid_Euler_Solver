#include "PostProcess.h"
#include <matplot/matplot.h>
#include <limits>

using namespace matplot;

// Filled contour plot with mesh overlay
#include "PostProcess.h"
#include <matplot/matplot.h>

using namespace matplot;

// Filled contour plot with mesh overlay
void plotFieldContour(const GridHandler& grid,
    const Eigen::MatrixXd& field,
    const std::string& titleText)
{
// --- 1) build nested vectors for contourf
const auto& xc = grid.getXCenter();    // ni×nj
const auto& yc = grid.getYCenter();
int ni = xc.rows();
int nj = xc.cols();

std::vector<std::vector<double>> Xc(ni, std::vector<double>(nj));
std::vector<std::vector<double>> Yc(ni, std::vector<double>(nj));
std::vector<std::vector<double>> Zc(ni, std::vector<double>(nj));
for (int i = 0; i < ni; ++i) {
for (int j = 0; j < nj; ++j) {
Xc[i][j] = xc(i, j);
Yc[i][j] = yc(i, j);
Zc[i][j] = field(i, j);
}
}

// --- 2) draw filled contours
figure();
contourf(Xc, Yc, Zc, 30);
colormap(matplot::palette::viridis());

// --- 3) add & position the colorbar
colorbar();
// tweak these numbers to taste
gca()->cb_position({
0.92f, // left
0.1f,  // bottom
0.05f, // width
0.8f   // height
});

hold(on);

// --- 4) overlay your structured‐grid lines
const auto& fullX = grid.getX();
const auto& fullY = grid.getY();
int nx = fullX.rows();
int ny = fullX.cols();

std::vector<double> X_rows, Y_rows, X_cols, Y_cols;
X_rows.reserve(nx * (ny + 1));
Y_rows.reserve(nx * (ny + 1));
X_cols.reserve(ny * (nx + 1));
Y_cols.reserve(ny * (nx + 1));

// horizontal lines
for (int j = 0; j < ny; ++j) {
for (int i = 0; i < nx; ++i) {
X_rows.push_back(fullX(i, j));
Y_rows.push_back(fullY(i, j));
}
X_rows.push_back(std::numeric_limits<double>::quiet_NaN());
Y_rows.push_back(std::numeric_limits<double>::quiet_NaN());
}
// vertical lines
for (int i = 0; i < nx; ++i) {
for (int j = 0; j < ny; ++j) {
X_cols.push_back(fullX(i, j));
Y_cols.push_back(fullY(i, j));
}
X_cols.push_back(std::numeric_limits<double>::quiet_NaN());
Y_cols.push_back(std::numeric_limits<double>::quiet_NaN());
}

plot(X_rows, Y_rows, "k")->line_width(0.5);
plot(X_cols, Y_cols, "k")->line_width(0.5);

// --- 5) decorate and show
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
