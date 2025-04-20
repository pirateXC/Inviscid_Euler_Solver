#ifndef FLUXSTATE_H
#define FLUXSTATE_H

#include <Eigen/Dense>
#include <array>

class FluxState {
public:
    enum Conserved { RHO = 0, RHO_U = 1, RHO_V = 2, ENERGY = 3 };

    FluxState(int ni_, int nj_)
      : ni(ni_), nj(nj_)
    {
        // initialize state vector
        for (auto &mat : Q) {
            mat.setZero(ni, nj);
        }

        // initialize primatives
        P   = Eigen::MatrixXd::Zero(ni, nj);
        T   = Eigen::MatrixXd::Zero(ni, nj);
        u   = Eigen::MatrixXd::Zero(ni, nj);
        v   = Eigen::MatrixXd::Zero(ni, nj);

    }

    void packToQ(double R, double gamma)
    {
        Eigen::ArrayXXd rho   = (P.array() / T.array()) / R;
        Eigen::ArrayXXd rho_u = rho * u.array();
        Eigen::ArrayXXd rho_v = rho * v.array();
        Eigen::ArrayXXd energy = (P.array()/(gamma - 1.0))
                              + 0.5 * rho * (u.array().square() + v.array().square());

        Q[RHO].block(1, 1, ni-2, nj-2) = rho.matrix().block(1, 1, ni-2, nj-2);
        Q[RHO_U].block(1, 1, ni-2, nj-2) = rho_u.matrix().block(1, 1, ni-2, nj-2);
        Q[RHO_V].block(1, 1, ni-2, nj-2) = rho_v.matrix().block(1, 1, ni-2, nj-2);
        Q[ENERGY].block(1, 1, ni-2, nj-2) = energy.matrix().block(1, 1, ni-2, nj-2);
    }

    // unpack state vector to primitive fields after update
    void unpackFromQ(double R, double gamma) {
        P = computePressure(R, gamma);
        T = computeTemp(R, gamma);
        u = (Q[RHO_U].array() / Q[RHO].array()).matrix();
        v = (Q[RHO_V].array() / Q[RHO].array()).matrix();
    }

    Eigen::MatrixXd computePressure(double R, double gamma) const {
        auto rho   = Q[RHO].array();
        auto rho_u = Q[RHO_U].array();
        auto rho_v = Q[RHO_V].array();
        auto energy= Q[ENERGY].array();

        return ((gamma - 1.0) *
                (rho * energy
                 - 0.5 * rho * ((rho_u/rho).square() + (rho_v/rho).square())))
               .matrix();
    }

    Eigen::MatrixXd computeTemp(double R, double gamma) const {
        return (computePressure(R, gamma).array() / Q[RHO].array() / R).matrix();
    }

    Eigen::MatrixXd computeVeloU() const {
        return (Q[RHO_U].array() / Q[RHO].array()).matrix();
    }

    Eigen::MatrixXd computeVeloV() const {
        return (Q[RHO_V].array() / Q[RHO].array()).matrix();
    }

    // getters
    std::array<Eigen::MatrixXd,4>& getQ() { return Q; }
    const std::array<Eigen::MatrixXd,4>& getQ() const { return Q; }
    Eigen::MatrixXd& pressure() { return P; }
    Eigen::MatrixXd& temperature() { return T; }
    Eigen::MatrixXd& velocityU() { return u; }
    Eigen::MatrixXd& velocityV() { return v; }

private:
    int ni; // num of i-components (including halo cells)
    int nj; // num of j-components (including halo cells)
    Eigen::MatrixXd P, T, u, v; // primatives pressure, temp, u-velo, v-velo
    std::array<Eigen::MatrixXd,4> Q;  // [rho, rho*u, rho*v, rho*E]
};

#endif // FLUXSTATE_H
