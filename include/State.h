// include/State.h

#ifndef STATE_H
#define STATE_H

#include <Eigen/Dense>
#include <array>

class State {
public:
    // index enum for Q components
    enum Conserved { RHO = 0, RHO_U = 1, RHO_V = 2, ENERGY = 3 };

    // build a zeroed Q‐field sized (ni × nj)
    State(int ni, int nj) {
        for (auto &mat : Q)
            mat.setZero(ni, nj);
    }

    // switch from primitives (P, u, v, T) to conserved Q
    void packToQ(const Eigen::MatrixXd &P,
                 const Eigen::MatrixXd &u,
                 const Eigen::MatrixXd &v,
                 const Eigen::MatrixXd &T,
                 double R,
                 double gamma)
    {
        Eigen::ArrayXXd rho   = (P.array() / T.array()) / R;
        Eigen::ArrayXXd rho_u = rho * u.array();
        Eigen::ArrayXXd rho_v = rho * v.array();
        Eigen::ArrayXXd energy = (P.array()/(gamma - 1.0))
                              + 0.5 * rho * (u.array().square() + v.array().square());

        Q[RHO]    = rho.matrix();
        Q[RHO_U]  = rho_u.matrix();
        Q[RHO_V]  = rho_v.matrix();
        Q[ENERGY] = energy.matrix();
    }

    // extract primitives from Q
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
        return (computePressure(R, gamma).array()
                / Q[RHO].array() / R)
               .matrix();
    }

    Eigen::MatrixXd computeVeloU() const {
        return (Q[RHO_U].array() / Q[RHO].array()).matrix();
    }

    Eigen::MatrixXd computeVeloV() const {
        return (Q[RHO_V].array() / Q[RHO].array()).matrix();
    }

    // mutable accessor so Initialize can write into Q
    std::array<Eigen::MatrixXd,4>& getQ() { return Q; }
    // const accessor for read‐only users
    const std::array<Eigen::MatrixXd,4>& getQ() const { return Q; }

private:
    std::array<Eigen::MatrixXd,4> Q;  // [rho, rho*u, rho*v, rho*E]
};

#endif // STATE_H
