#include "StegerWarming.h"
#include <cmath>

namespace StegerWarming {

// unpack one state vector into primitives & derived quantities
struct Prim { double rho,u,v,E,p,a,ho,ek; };

inline Prim unpackQ(const Eigen::Vector4d &Q, double gamma) {
    Prim Qv;
    Qv.rho = Q(0);
    Qv.u   = Q(1) / Qv.rho;
    Qv.v   = Q(2) / Qv.rho;
    Qv.E   = Q(3) / Qv.rho;
    Qv.ek  = 0.5 * (Qv.u*Qv.u + Qv.v*Qv.v);
    Qv.p   = (gamma - 1.0) * (Qv.rho*Qv.E - Qv.rho*Qv.ek);
    Qv.a   = std::sqrt(gamma * Qv.p / Qv.rho);
    Qv.ho  = Qv.E + Qv.p / Qv.rho;
    return Qv;
}

Eigen::Vector4d fluxDotNorm(
    const Eigen::Vector4d &Q,
    double gamma,
    const Eigen::Vector2d &normal)
{
    Prim Qv = unpackQ(Q, gamma);
    double un = Qv.u * normal(0) + Qv.v * normal(1);

    Eigen::Vector4d F;
    // mass flux
    F(0) = Qv.rho * un;
    // x‑momentum flux
    F(1) = Qv.rho * un * Qv.u
         + Qv.p * normal(0);
    // y‑momentum flux
    F(2) = Qv.rho * un * Qv.v
         + Qv.p * normal(1);
    // energy flux
    F(3) = (Qv.rho * Qv.E + Qv.p) * un;
    return F;
}

Eigen::Vector4d computeDissipLeft(
    const Eigen::Vector4d &QL,
    double gamma,
    const Eigen::Vector2d &normal)
{
    double eps = 0.2; // accounts for sonic points, may have to alter

    Prim Qv = unpackQ(QL, gamma);
    Eigen::Vector2d t(-normal.y(), + normal.x());

    double un = Qv.u*normal.x() + Qv.v*normal.y();

    // Characteristic speeds
    double lam1 = un - Qv.a, lam2 = un,
           lam3 = un + Qv.a, lam4 = un;

    double sq1 = std::sqrt(lam1*lam1 + eps*eps);
    double lam1L = 0.5*( lam1 + sq1 );

    double sq2 = std::sqrt(lam2*lam2 + eps*eps);
    double lam2L = 0.5*( lam2 + sq2 );

    double sq3 = std::sqrt(lam3*lam3 + eps*eps);
    double lam3L = 0.5*( lam3 + sq3 );

    double sq4 = std::sqrt(lam4*lam4 + eps*eps);
    double lam4L = 0.5*( lam4 + sq4 );

    // Left‐eigenvectors (rows ℓ_i of L)
    double denom = 2.0 * Qv.a * Qv.a;
    double L1_0 = ((gamma-1.0)*Qv.ek + Qv.a*un)/denom;
    double L1_1 = ((1.0-gamma)*Qv.u - Qv.a*normal.x())/denom;
    double L1_2 = ((1.0-gamma)*Qv.v - Qv.a*normal.y())/denom;
    double L1_3 = (gamma - 1.0)/denom;

    double L2_0 = (Qv.a*Qv.a - (gamma-1.0)*Qv.ek)/(Qv.a*Qv.a);
    double L2_1 = ((gamma-1.0)*Qv.u)/(Qv.a*Qv.a);
    double L2_2 = ((gamma-1.0)*Qv.v)/(Qv.a*Qv.a);
    double L2_3 = (1.0 - gamma)/(Qv.a*Qv.a);

    double L3_0 = ((gamma-1.0)*Qv.ek - Qv.a*un)/denom;
    double L3_1 = ((1.0-gamma)*Qv.u + Qv.a*normal.x())/denom;
    double L3_2 = ((1.0-gamma)*Qv.v + Qv.a*normal.y())/denom;
    double L3_3 = (gamma - 1.0)/denom;

    // shear‐mode in rotated frame (no singularities)
    double L4_0 = 0.0;
    double L4_1 =  t.x();
    double L4_2 =  t.y();
    double L4_3 = 0.0;

    // Project Q into characteristic space
    Eigen::Vector4d w;
    w(0) = L1_0*QL(0) + L1_1*QL(1) + L1_2*QL(2) + L1_3*QL(3);
    w(1) = L2_0*QL(0) + L2_1*QL(1) + L2_2*QL(2) + L2_3*QL(3);
    w(2) = L3_0*QL(0) + L3_1*QL(1) + L3_2*QL(2) + L3_3*QL(3);
    w(3) = L4_0*QL(0) + L4_1*QL(1) + L4_2*QL(2) + L4_3*QL(3);

    // Apply dissipation
    w(0) *= lam1L;
    w(1) *= lam2L;
    w(2) *= lam3L;
    w(3) *= lam4L;

    // Right‐eigenvectors (columns r_i of R)
    Eigen::Vector4d r1, r2, r3, r4;
    r1 << 1.0,
          (un - Qv.a)*normal.x(),
          (un - Qv.a)*normal.y(),
          Qv.ho - Qv.a*un;

    r2 << 1.0,
          Qv.u,
          Qv.v,
          Qv.ek;

    r3 << 1.0,
          (un + Qv.a)*normal.x(),
          (un + Qv.a)*normal.y(),
          Qv.ho + Qv.a*un;

    // shear‐mode in rotated frame
    r4 << 0.0,
        t.x(),
        t.y(),
        Qv.u * t.x() + Qv.v * t.y();

    return r1*w(0) + r2*w(1) + r3*w(2) + r4*w(3);
}

Eigen::Vector4d computeDissipRight(
    const Eigen::Vector4d &QR,
    double gamma,
    const Eigen::Vector2d &normal)
{
    double eps = 0.2; // accounts for sonic points, may have to alter

    Prim Qv = unpackQ(QR, gamma);
    Eigen::Vector2d t(-normal.y(), + normal.x());

    double un = Qv.u*normal.x() + Qv.v*normal.y();

    // Characteristic speeds
    double lam1 = un - Qv.a, lam2 = un,
           lam3 = un + Qv.a, lam4 = un;

    double sq1 = std::sqrt(lam1*lam1 + eps*eps);
    double lam1R = 0.5*( lam1 - sq1 );

    double sq2 = std::sqrt(lam2*lam2 + eps*eps);
    double lam2R = 0.5*( lam2 - sq2 );

    double sq3 = std::sqrt(lam3*lam3 + eps*eps);
    double lam3R = 0.5*( lam3 - sq3 );

    double sq4 = std::sqrt(lam4*lam4 + eps*eps);
    double lam4R = 0.5*( lam4 - sq4 );

    // Left‐eigenvectors (rows ℓ_i of L)
    double denom = 2.0 * Qv.a * Qv.a;
    double L1_0 = ((gamma-1.0)*Qv.ek + Qv.a*un)/denom;
    double L1_1 = ((1.0-gamma)*Qv.u - Qv.a*normal.x())/denom;
    double L1_2 = ((1.0-gamma)*Qv.v - Qv.a*normal.y())/denom;
    double L1_3 = (gamma - 1.0)/denom;

    double L2_0 = (Qv.a*Qv.a - (gamma-1.0)*Qv.ek)/(Qv.a*Qv.a);
    double L2_1 = ((gamma-1.0)*Qv.u)/(Qv.a*Qv.a);
    double L2_2 = ((gamma-1.0)*Qv.v)/(Qv.a*Qv.a);
    double L2_3 = (1.0 - gamma)/(Qv.a*Qv.a);

    double L3_0 = ((gamma-1.0)*Qv.ek - Qv.a*un)/denom;
    double L3_1 = ((1.0-gamma)*Qv.u + Qv.a*normal.x())/denom;
    double L3_2 = ((1.0-gamma)*Qv.v + Qv.a*normal.y())/denom;
    double L3_3 = (gamma - 1.0)/denom;

    // shear‐mode in rotated frame (no singularities)
    double L4_0 = 0.0;
    double L4_1 =  t.x();
    double L4_2 =  t.y();
    double L4_3 = 0.0;

    // Project Q into characteristic space
    Eigen::Vector4d w;
    w(0) = L1_0*QR(0) + L1_1*QR(1) + L1_2*QR(2) + L1_3*QR(3);
    w(1) = L2_0*QR(0) + L2_1*QR(1) + L2_2*QR(2) + L2_3*QR(3);
    w(2) = L3_0*QR(0) + L3_1*QR(1) + L3_2*QR(2) + L3_3*QR(3);
    w(3) = L4_0*QR(0) + L4_1*QR(1) + L4_2*QR(2) + L4_3*QR(3);

    // Apply dissipation
    w(0) *= lam1R;
    w(1) *= lam2R;
    w(2) *= lam3R;
    w(3) *= lam4R;

    // Right‐eigenvectors (columns r_i of R)
    Eigen::Vector4d r1, r2, r3, r4;
    r1 << 1.0,
          (un - Qv.a)*normal.x(),
          (un - Qv.a)*normal.y(),
          Qv.ho - Qv.a*un;

    r2 << 1.0,
          Qv.u,
          Qv.v,
          Qv.ek;

    r3 << 1.0,
          (un + Qv.a)*normal.x(),
          (un + Qv.a)*normal.y(),
          Qv.ho + Qv.a*un;

    // shear‐mode in rotated frame
    r4 << 0.0,
        t.x(),
        t.y(),
        Qv.u * t.x() + Qv.v * t.y();;

    return r1*w(0) + r2*w(1) + r3*w(2) + r4*w(3);
}

// Symmetric Steger–Warming flux with mask to disable dissipation at walls
Eigen::Vector4d computeFaceFlux(
    const Eigen::Vector4d &QL,
    const Eigen::Vector4d &QR,
    double gamma,
    const Eigen::Vector2d &normal,
    int mask)
{
    Eigen::Vector4d FL = fluxDotNorm(QL, gamma, normal);
    Eigen::Vector4d FR = fluxDotNorm(QR, gamma, normal);

    // Exact dissipation
    Eigen::Vector4d AQL = computeDissipLeft(QL, gamma, normal);
    Eigen::Vector4d AQR = computeDissipRight(QR, gamma, normal);


    Eigen::Vector4d Fc = 0.5*(FL + FR);
    Eigen::Vector4d D  = 0.5*(AQR - AQL);
    Eigen::Vector4d Dmask = mask * D;

    return Fc - Dmask;

    // Symmetric flux split
   // return 0.5*(FL + FR) - 0.5 * mask * (AQR - AQL);
}

} // namespace StegerWarming
