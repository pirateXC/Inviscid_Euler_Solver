#ifndef STEGER_WARMING_H
#define STEGER_WARMING_H

#include <Eigen/Dense>

namespace StegerWarming {

  // compute AL*QL dissipation
  Eigen::Vector4d computeDissipLeft(
    const Eigen::Vector4d &QL,
    double gamma,
    const Eigen::Vector2d &normal);

  // compute AR*QR dissipation
  Eigen::Vector4d computeDissipRight(
    const Eigen::Vector4d &QR,
    double gamma,
    const Eigen::Vector2d &normal);

  // symmetric Steger–Warming flux
  // 'mask' = 0 disables dissipation (for wall faces), =1 applies it
  Eigen::Vector4d computeFaceFlux(
    const Eigen::Vector4d &QL,
    const Eigen::Vector4d &QR,
    double gamma,
    const Eigen::Vector2d &normal,
    int mask);

    Eigen::Vector4d fluxDotNorm(
        const Eigen::Vector4d &Q,
        double gamma,
        const Eigen::Vector2d &normal);


} // namespace StegerWarming

#endif // STEGER_WARMING_H
