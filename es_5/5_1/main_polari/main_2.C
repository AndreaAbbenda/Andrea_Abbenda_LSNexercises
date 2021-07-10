#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

double AEW(int n, double r, double theta){
	if (n == 1){ return pow(exp(-r)/sqrt(M_PI),2);}
	else if (n == 2){ return pow(sqrt(2/M_PI)/8*r*exp(-r/2)*cos(theta),2);}
	}
	
double EvalGauss(double mean, double sigma, double x){
	return (1/sqrt(2*M_PI*sigma*sigma))*exp(-pow((x-mean),2)/(2*sigma*sigma));
	}

double Minimum (double a, double b){ if (a > b){return b; }else if(b >= a){return a;}}


int main(int argc, char** argv){
   
   Random rnd;
   int seed[3];
   int p1, p2;
   ifstream Primes("Primes");   
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");  
   string property;            
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
   
	int nsteps = 10000;
	double sigma_r = 1;
	double r1_old = abs(rnd.Gauss(0,sigma_r));
	double theta1_old = 2*M_PI*rnd.Rannyu();
	double r2_old = abs(rnd.Gauss(0,sigma_r));
	double theta2_old = 2*M_PI*rnd.Rannyu();
	double A1=0, A2=0;
	int sumA1=0, sumA2=0;
	ofstream N1("es5_1_N1.dat");
	ofstream N2("es5_1_N2.dat");
	
	double r1_new = 0;
	double theta1_new = 0;
	double r2_new = 0;
	double theta2_new = 0;
	double test_random = rnd.Rannyu();
	
	for (int i=1; i <= nsteps; i++){
		if (i % 1000 == 0){ cout << "step n°"<< i <<endl; }
		
		while (A1 < test_random){
			r1_new = rnd.Gauss(r1_old,sigma_r);
			theta1_new = 2*M_PI*rnd.Rannyu();
				
			double Q1 = AEW(1, abs(r1_new), theta1_new)/AEW(1, r1_old, theta1_old);
			A1 = Minimum (1, Q1);
			test_random = rnd.Rannyu();
			sumA1 ++;
			}
		A1 = 0;
		double phi = rnd.Rannyu()*2*M_PI;
		N1 << r1_old << "	" << theta1_old << "	"<< phi << endl;
		r1_old = abs(r1_new);
		theta1_old = theta1_new;
	
		while (A2 < test_random){
			r2_new = rnd.Gauss(r2_old,sigma_r);
			theta2_new = 2*M_PI*rnd.Rannyu();
				
			double Q2 = AEW(2, abs(r2_new), theta2_new)/AEW(2, r2_old, theta2_old);
			A2 = Minimum (1, Q2);
			test_random = rnd.Rannyu();
			sumA2 ++;
			}
		A2 = 0;
		phi = rnd.Rannyu()*2*M_PI;
		N2 << r2_old << "	" << theta2_old << "	"<< phi << endl;
		r2_old = abs(r2_new);
		theta2_old = theta2_new;
		
	}
	cout << "Probability of accepting a step (quantic number n = 1): " << (double)nsteps/sumA1 << endl;
	cout << "Probability of accepting a step (quantic number n = 2): " << (double)nsteps/sumA2 << endl;
	N1.close();
	N2.close();	
		
return 0;
}
