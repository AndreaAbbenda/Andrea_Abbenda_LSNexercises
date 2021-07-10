#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"


using namespace std;

double Boltzmann(double Hnew, double Hold, double t){
	return exp(-(Hnew-Hold)/t);
	}

int main(int argc, char** argv){
   
   Random rnd;
   int seed[4];
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
  
 
 	TApplication app("app", &argc,argv);
	TH2D *isto= new TH2D("a","b",200,0,2,200,0,2);

 
  
  double range_sigma = 0.5;
  double range_mu = 0.5;
  int H_steps = 5000;
  double H_sigma = 1.5;
  
  
  double Tmax = 1;
  int n_t = 50;
  int t_steps = 1000;
  
  double sigma_0 = 1;
  double mu_0 = 1;
  
  WaveFunction Phi_test;  
  double H_old, H_new;
  H_old = pow(Phi_test.AveH(H_steps, rnd, H_sigma) + 40,2);
  
  double sigma_test, mu_test;
  double A, Q, test;
  
  int attempted = 0;
  int accepted = 0;
  
  
  int trig =0;
  ofstream sm("sigmamu.dat");
  
  for (int i=0; i<n_t; i++){
  	
	double t = Tmax - i*(Tmax/(double)n_t);
	
	for (int j=0; j<t_steps; j++){
		
		if (trig == 0){
	//		sigma_test = Maximum(0.1, abs(sigma_0 + range_sigma*(rnd.Rannyu()-0.5)));
			sigma_test = abs(sigma_0 + range_sigma*(rnd.Rannyu()-0.5));
			mu_test = mu_0;
			Phi_test.SetSigma(sigma_test);
			trig = 1;
		}else if(trig == 1){
			mu_test = abs(mu_0 + range_mu*(rnd.Rannyu()-0.5));
			sigma_test = sigma_0;
			Phi_test.SetMu(mu_test);
			trig = 0;
			}
		
			
		H_new = Phi_test.AveH(H_steps, rnd, H_sigma);
		Q = Boltzmann(H_new, H_old, t);
		A = Minimum(1, Q);
		test = rnd.Rannyu();
		
		if(A > test){
			sigma_0 = sigma_test;
			mu_0 = mu_test;
			H_old = H_new;
			accepted++;
			}
		attempted++;
		isto->Fill(mu_0, sigma_0);
		sm << mu_0 << "          " << sigma_0 << endl;
	}
	
	cout << "Probability of accepting a move at t = "<<t<<" : " << (double)accepted/(double)attempted << endl;
	cout<< "energia: "<< H_old<<endl;
	cout<< "mu: "<< mu_0<<endl;
	cout<< "sigma: "<< sigma_0<<endl<<endl;
	accepted = 0;
	attempted = 0;
  }
  
  cout<< " best mu: "<< mu_0<< endl;
  cout<< " best sigma: "<< sigma_0<< endl;
  sm.close();
  
	int nsteps = 100000;
	int blocks = 100;
	int steps = nsteps/blocks;
	
	WaveFunction Phi(mu_0,sigma_0);
	double sigma = 1.5;
	double x_0 = rnd.Rannyu();
	double _try;
	
	attempted = 0;
	accepted = 0;
	A = 0;
	double test_random = rnd.Rannyu();
	
	double sum_h = 0;
	double ave_h[blocks];              //errori
	double ave2_h[blocks];
	double err_h[blocks];	
	
	ofstream sam("sampling_best.dat");
	ofstream ham("hamilton_best.dat");
	
		
	for (int i=1; i <= nsteps; i++){ 
		while (A < test_random){	
			_try = rnd.Gauss(x_0, sigma);	 
			Q = Phi.Eval2(_try)/Phi.Eval2(x_0);
			A = Minimum (1, Q);
			test_random = rnd.Rannyu();
			attempted ++;
			}
		x_0 = _try;	
		A = 0;
		accepted ++;
		double h = Phi.H(_try);
		double p = Phi.Eval(_try);
		sam << _try << "	" << p << "	"<< h/p << endl;

			
		sum_h += h/p;		
		if (i % steps == 0){
			int j = i/steps;
			if (j == 1){
				ave_h[0] = sum_h/steps;
				ave2_h[0] = pow(ave_h[0],2);
				err_h[0] = 0;
			}else {
				ave_h[j-1] = (ave_h[j-2]*(j-1) + sum_h/steps)/j;
				ave2_h[j-1] = (ave2_h[j-2]*(j-1) + pow(sum_h/steps,2))/j;
				err_h[j-1] = sqrt((ave2_h[j-1]-pow(ave_h[j-1],2))/(j-1));
				}
			ham<<ave_h[j-1]<<"	"<<err_h[j-1]<<endl;
			sum_h = 0;
			if (j==blocks){cout<<"energia : "<<ave_h[j-1]<<endl;}
			}						
	}
	cout << "Probability of accepting a step, multivariate normal transition probability : " << (double)accepted/(double)attempted << endl;
	ham.close();
	sam.close();
	
	TCanvas *c1=new TCanvas("c1","plot");
	c1->cd();
	isto->Draw();
	app.Run();

			
return 0;
}
