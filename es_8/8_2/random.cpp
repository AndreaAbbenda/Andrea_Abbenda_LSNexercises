/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "random.h"


using namespace std;

double err(double *av2, double *ave, int N){
	if (N==0){return 0;}
	else {return sqrt(abs(av2[N]-pow(ave[N],2))/N);}
	}
	
double Maximum (double a, double b){ 
	if (a > b){return a; 
	}else if(b >= a){return b;
		}
	}

double Minimum (double a, double b){ 
	if (a > b){return b; 
	}else if(b >= a){return a;
		}
	}


Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

WaveFunction :: WaveFunction(){
	_mu=1;
	_sigma=1;
	}

WaveFunction :: WaveFunction(double mu, double sigma){
	_mu=mu;
	_sigma=sigma;
	}

	
WaveFunction :: ~WaveFunction(){}

void WaveFunction :: Copy(WaveFunction phi){
	this->SetMu(phi.GetMu());
	this->SetSigma(phi.GetSigma());
	}

double WaveFunction :: Eval(double x){
	double sigma2 = _sigma*_sigma;
	double a = exp(-pow(x-_mu,2)/(2*sigma2)) + exp(-pow(x+_mu,2)/(2*sigma2));
	double norm = sqrt(2 * _sigma * sqrt(M_PI)* (1+exp(-_mu*_mu/sigma2)));
	return a / norm;
	}
	
double WaveFunction :: Eval2(double x){
	double a = Eval(x);
	return a*a;
	}

double WaveFunction :: D2Eval(double x){
	double sigma2 = _sigma*_sigma;
	double sigma4 = _sigma*_sigma*_sigma*_sigma;
	double a = exp(-pow(x-_mu,2)/(2*sigma2)) + exp(-pow(x+_mu,2)/(2*sigma2));
	double b = exp(-pow(x+_mu,2)/(2*sigma2)) - exp(-pow(x-_mu,2)/(2*sigma2));
	double norm = sqrt(2 * _sigma * sqrt(M_PI)* (1+exp(-_mu*_mu/sigma2)));
	return (a * (x*x + _mu*_mu - sigma2)/sigma4 + b * 2*x*_mu/sigma4) / norm;
	}

double WaveFunction :: H(double x){
	double Kin = -D2Eval(x)/2;
	double Pot = x*x*x*x - 5*x*x/2;
	return Kin + Pot*Eval(x);
	}

double WaveFunction :: AveH(int steps, Random rnd, double sigma){
	
	double A, Q, test_random, _try, x_0, sum_h, h, p;
	x_0 = rnd.Rannyu();
	sum_h = 0;
	A = 0;
	
	for (int i=0; i<steps; i++){
		
	        test_random = rnd.Rannyu();
	        	        
	        while (A < test_random){	
			_try = rnd.Gauss(x_0, sigma);	 
			Q = Eval2(_try)/Eval2(x_0);
			A = Minimum (1, Q);
			test_random = rnd.Rannyu();
			}
		A = 0;
		x_0 = _try;
		
		h = H(_try);
		p = Eval(_try);
		
		sum_h += h/p;
	        }
	        
	return sum_h/(double)steps;
	}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
