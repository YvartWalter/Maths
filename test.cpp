

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <time.h>
#include <cstring>

//STD C++
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <limits>
#include <complex>
#include <valarray>
#include <random>
#include <chrono>
#include <Eigen/Dense>

#include "B_functions.hpp"

template<class Engine=std::mt19937_64>
class RndU
{
  public:
    RndU(void):rd(), dist(0.0,1.0), gen( rd()) { };
    RndU(const double& min, const double& max):rd(), dist(min,max), gen( rd()) { };
    RndU(const std::string& token, const double& min, const double& max):rd(token), dist(min,max), gen( rd()) { };
    RndU(const RndU<Engine>& r):rd(r.rd), dist(r.dist.param()), gen(rd()) { };
    virtual ~RndU() { };
    double operator()(void) { return dist(gen); };
    void seed(uint_fast64_t sd){gen.seed(sd);};
  protected:
    std::random_device rd;
    std::uniform_real_distribution<double> dist;
    Engine gen;
};


using namespace std;
using namespace Eigen;

int main(int argc, char* argv[])
{
  
  RndU<std::mt19937_64> Uniform;
  //Uniform.seed( std::chrono::system_clock::now().time_since_epoch().count() );
  d_invpow Func;
  d_polynom Pln(3);
  
  unsigned int dim=25, degre = Pln; //degré réel + 1
  unsigned int i,j,k;
  VectorXd Y(dim), X(dim), Poly(degre), SY(dim), Z(dim), Y0, P2(degre);
  
  // dim = 25
  Y << 147.419 , 134.429 , 122.749 , 112.928 , 104.412 , 96.7632 , 89.8072 , 83.4222 , 77.5066 , 71.9985 , 66.8701 , 62.1014 , 57.6506 , 53.4605 , 49.4906, 45.7377 , 42.2168 , 38.9176 , 35.8143 , 32.8988 , 30.1586 , 27.5675 , 25.1187 , 22.9306 , 21.4945;
  Y0 = Y;
   
  
  // Func.Redefine(1.0, 1.0);
  for(i=0;i<dim;++i)
  {
    double v = Uniform(), w = Uniform();
    double I = (double)i;
    X(i) = 0.1+I/10.0;
    Y(i) += 20.0*(v-w);
    SY(i) = 1.0;
  }
  
//   bool L=Func;
  double A=0.0,B=0.0;
//   unsigned int n_par = Func;
  Pln.Parameters( Poly.data(), Poly.size() );
//   unsigned int n_par = Poly.size();
  cout<<"POLY:\n"<<Poly<<endl;
  P2=Poly;
//   Func.Parameters(A,B);
//   Poly.resize( n_par );
//   P2.resize( n_par );
//   for(unsigned int i=0;i<n_par;++i) Poly(i) = Func[i]; //paramètres initiaux
//   Func.Retrieve(X.data(), Y.data(), X.size());
//   A = LinFit<d_invpow, Eigen::MatrixXd, Eigen::VectorXd>( Func,X.size(), X.data(), Y.data(), n_par, Poly.data(), P2.data(), SY.data(), Z.data(), L);
//   B = Func.Retrieve(X.data(), Y.data(), X.size());
     B = Pln.Retrieve<Eigen::MatrixXd, Eigen::VectorXd>(X.data(), Y.data(), Z.data(), X.size());
     Pln.Parameters( P2.data(), P2.size() );
//   Func.Parameters(A,B);
  cout<<"i  X(i)  Y(i)  Sol  F_sol(X)(i)  Y0(i)"<<endl;
  for(i=0;i<dim;++i)
  {
    double I=(double)i;
//     cout<<i<<"  "<<X(i)<<"  "<<Y(i)<<"   "<<Func(X(i))<<"   "<<Z(i)<<"   "<<Y0(i)<<endl;
    cout<<i<<"  "<<X(i)<<"  "<<Y(i)<<"   "<<Pln(X(i))<<"   "<<Z(i)<<"   "<<Y0(i)<<endl;
//     cout<<i<<"  "<<X(i)<<"  "<<Y(i)<<"   "<<0<<"   "<<Z(i)<<"   "<<Y0(i)<<endl;
  }
  cout<<endl;
  
  return EXIT_SUCCESS;
}





