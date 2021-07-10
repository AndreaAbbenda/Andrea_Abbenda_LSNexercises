#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"


using namespace std;


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
  
	int nsteps = 100000;
	int blocks = 100;
	int steps = nsteps/blocks;
	double mu_in, sigma_in;
	cout<<"mu?"<<endl;
	cin>>mu_in;
	cout<<"sigma?"<<endl;
	cin>>sigma_in;
	WaveFunction Phi(mu_in,sigma_in);
	
	double sigma = 1.0;
	double x_0 = rnd.Rannyu();
	
	int attempted = 0;
	int accepted = 0;
	double A = 0;
	double test_random = rnd.Rannyu();
	double _try, Q;
	
	double sum_h = 0;
	double ave_h[blocks];              //errori
	double ave2_h[blocks];
	double err_h[blocks];	
	
	ofstream sam("sampling.dat");
	ofstream ham("hamilton.dat");
	
	
	TApplication app("app", &argc,argv);
	double range = Phi.GetMu()+3*Phi.GetSigma();
	TH1D *isto= new TH1D("a","b",200,-range,range);

	
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
		double h = Phi.Ham(_try);
		double p = Phi.Eval(_try);
		double hp = Phi.H_p(_try);
		sam << _try << "          " << p << "          "<<hp<<"          "<< h/p << endl;
	
		isto->Fill(_try);	
			
		sum_h += h/p;		
		if (i % steps == 0){                                 //calcolo degli errori statistici
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
