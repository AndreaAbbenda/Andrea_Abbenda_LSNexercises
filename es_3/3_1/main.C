#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

 
int main(int argc, char** argv){
   
   Random rnd;
   int seed[5];
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


	int Massets = 10000;          
	int Nblocks = 100;
	
	double S0 = 100;                   // definisco i parametri necessari per la stima degli asset prices
	double Tdel = 1.0;
	double K = 100;
	double r = 0.1;
	double sigma = 0.25;
	
	double ave[4][Nblocks];            // definisco le variabili per il data-blocking
	double av2[4][Nblocks];
	double err_prog[4][Nblocks];
	double sum[4];
	for (int i=0; i<4; i++){ave[i][0]=0; av2[i][0]=0; err_prog[i][0]=0; sum[i]=0;}
	
	
	ofstream o_call("es3_1_call.txt");
	ofstream o_put("es3_1_put.txt");
	
	
	for (int i=0; i<Nblocks; i++){  //ciclo sui blocchi
		for (int j=0; j<4; j++){sum[j]=0;}
		cout<<"Block "<<i+1<<" of "<<Nblocks<<endl;
		for (int j=0; j<Massets; j++){  //ciclo sugli asset prices
		
		//calcolo il valore a t=T degli asset prizes col metodo direct e cumulative
		
			double direct=0;
			double cumulative=S0;
		
			double W = rnd.Gauss(0,Tdel);
			direct = S0*(exp((r-0.5*pow(sigma,2))*Tdel+sigma*W));   //direct method
			
			double Ti=Tdel/100;
			for (int i=0; i<100; i++){   //cumulative method
				cumulative=cumulative*exp((r-0.5*pow(sigma,2))*Ti+sigma*rnd.Gauss(0,1)*sqrt(Ti));
				}
			
			double callDir = exp(-r*Tdel)*Max(0, direct-K);           //stima delle opzioni put e call
			double callCum = exp(-r*Tdel)*Max(0, cumulative-K);
			double putDir = exp(-r*Tdel)*Max(0, K-direct);
			double putCum = exp(-r*Tdel)*Max(0, K-cumulative);
			
			sum[0] = sum[0] + callDir;
			sum[1] = sum[1] + callCum;
			sum[2] = sum[2] + putDir;
			sum[3] = sum[3] + putCum;			
			}
		
		for(int k=0; k<4; k++){
		
		//stima degli errori statistici  (data-blocking)
		
			if (i==0){ ave[k][0]=sum[k]/Massets; av2[k][0]=pow(ave[k][0],2); err_prog[k][0]=0;
			}else{
				ave[k][i]= (ave[k][i-1]*i + sum[k]/Massets)/(i+1);
				av2[k][i]= (av2[k][i-1]*i + pow(sum[k]/Massets,2))/(i+1);
				err_prog[k][i]= err(av2[k], ave[k], i);
				}
			}
		
		o_call<<ave[0][i]<<"	"<<err_prog[0][i]<<"		"<<ave[1][i]<<"	"<<err_prog[1][i]<<endl;
		o_put<<ave[2][i]<<"	"<<err_prog[2][i]<<"		"<<ave[3][i]<<"	"<<err_prog[3][i]<<endl;
		
	}
	
	o_call.close();
	o_put.close();
		
return 0;
}
