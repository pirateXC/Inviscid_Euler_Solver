#ifndef FLUXSTATE_H
#define FLUXSTATE_H
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <Eigen/Dense>
#include <array>

#include <iostream>
class FluxState {
public:
    FluxState(int ni,int nj):ni(ni),nj(nj) {
        for(auto &m:Q)
            m.setZero(ni-1,nj-1);
            P.setZero(ni-1,nj-1);
            T.setZero(ni-1,nj-1);
            u.setZero(ni-1,nj-1);
            v.setZero(ni-1,nj-1);
    }

    void packToQ(double R,double gamma) {
        auto Tm = T.array();
        auto Pm = P.array();
        auto rho = Pm/(Tm*R);
        auto rho_u = rho*u.array();
        auto rho_v = rho*v.array();
        auto e = (Pm/(gamma-1.)) + 0.5*rho*(u.array().square()+v.array().square());

        Q[0]=rho.matrix(); 
        Q[1]=rho_u.matrix(); 
        Q[2]=rho_v.matrix(); 
        Q[3]=e.matrix();
    }

    void unpackFromQ(double R,double gamma) {
        auto rho = Q[0].array();

        P = computePressure(gamma);
        T = computeTemp(R,gamma);
        u = (Q[1].array()/rho).matrix();
        v = (Q[2].array()/rho).matrix();
    }

    Eigen::MatrixXd computePressure(double gamma) const {
        const auto rho  = Q[0].array();
        const auto rhoE = Q[3].array();
        const auto rhou = Q[1].array();
        const auto rhov = Q[2].array();
    
        const auto kin = 0.5 * (rhou.square() + rhov.square()) / rho;
        const auto p = (gamma - 1.0) * (rhoE - kin);
    
        return p.matrix();
    }

    Eigen::MatrixXd computeTemp(double R, double gamma) const {
        const auto p  = computePressure(gamma).array();
        return (p / Q[0].array() / R).matrix();
    }
    
    std::array<Eigen::MatrixXd,4>& getQ() { return Q; }
    const std::array<Eigen::MatrixXd,4>& getQ() const { return Q; }
    Eigen::MatrixXd& getPressure() { return P; }
    Eigen::MatrixXd& getTemp() { return T; }
    Eigen::MatrixXd& getVelo_U() { return u; }
    Eigen::MatrixXd& getVelo_V() { return v; }

private:
        int ni; // num of i-components (including halo cells)
        int nj; // num of j-components (including halo cells)
        Eigen::MatrixXd P, T, u, v; // primatives pressure, temp, u-velo, v-velo
        std::array<Eigen::MatrixXd,4> Q;
};
#endif