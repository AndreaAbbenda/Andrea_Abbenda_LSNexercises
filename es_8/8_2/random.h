/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>

#ifndef __Random__
#define __Random__

double err(double *av2, double *ave, int N);
double Minimum (double a, double b);
double Maximum (double a, double b);

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};

class WaveFunction {

private:
  double _mu, _sigma;

protected:

public:
  // constructors
  WaveFunction();
  WaveFunction(double mu, double sigma);
  // destructor
  ~WaveFunction();
  // methods
  double SetMu(double a){_mu=a;};
  double SetSigma(double b){_sigma=b;};
  double GetMu(){double a = _mu; return a;};
  double GetSigma(){double b = _sigma; return b;};
  double Eval(double x);
  double Eval2(double x);
  double D2Eval(double x);
  double H(double x);
  double AveH(int steps, Random rnd, double sigma);
  void Copy(WaveFunction phi);

};


#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
