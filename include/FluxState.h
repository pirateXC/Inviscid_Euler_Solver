#ifndef FLUXSTATE_H
#define FLUXSTATE_H
#include <Eigen/Dense>
#include <array>

class FluxState {
public:
    enum Conserved{RHO,RHO_U,RHO_V,ENERGY};
    FluxState(int ni,int nj):ni(ni),nj(nj){for(auto &m:Q)m.setZero(ni,nj);P.setZero(ni,nj);T.setZero(ni,nj);u.setZero(ni,nj);v.setZero(ni,nj);}    
    void packToQ(double R,double gamma){
        auto Tm=T.array().max(1e-6), Pm=P.array().max(1e-6);
        auto rho=Pm/Tm/R;
        auto rho_u=rho*u.array(), rho_v=rho*v.array();
        auto e=(Pm/(gamma-1.))+0.5*rho*(u.array().square()+v.array().square());
        Q[RHO]=rho.matrix(); Q[RHO_U]=rho_u.matrix(); Q[RHO_V]=rho_v.matrix(); Q[ENERGY]=e.matrix();
    }
    void unpackFromQ(double R,double gamma){
        auto rho_arr=Q[RHO].array().max(1e-8);
        P=computePressure(R,gamma).array().max(1e-6).matrix();
        T=computeTemp(R,gamma).array().max(1e-6).matrix();
        u=(Q[RHO_U].array()/rho_arr).matrix();
        v=(Q[RHO_V].array()/rho_arr).matrix();
    }
    Eigen::MatrixXd computePressure(double R,double g)const{
        auto r=Q[RHO].array(), ru=Q[RHO_U].array(), rv=Q[RHO_V].array(), E=Q[ENERGY].array();
        return ((g-1)*(r*E-0.5*r*((ru/r).square()+(rv/r).square()))).matrix();
    }
    Eigen::MatrixXd computeTemp(double R,double g)const{
        return (computePressure(R,g).array()/Q[RHO].array()/R).matrix();
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