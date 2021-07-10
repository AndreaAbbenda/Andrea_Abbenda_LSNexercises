#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>


using namespace std;

	
double cosPi(double x){
	return 0.5*M_PI*cos(0.5*M_PI*x);}
		
//la densità di probabilità scelta per migliorare l'integrale di M_PI/2*cos(M_PI*x) è d(x) = 2(1-x). il metodo per generare numeri casuali secondo la distribuzione scelta è stato implementato come funzione della classe Random nei file .h e .cpp col nome di "d_lin()"
 
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

	ofstream dist_rand("es2_1_random.txt");
//	ofstream dist_prob("es2_1_dprob.txt");
	ofstream dist_lin("es2_1_dlin.txt");
	
	int Mthrows = 100000;
	int Nblocks = 100;
	double sumrand = 0;
//	double sumprob = 0;
	double sumlin = 0;
	double Irand = 0;
//	double Iprob = 0;
	double Ilin = 0;
	double ave[3][Nblocks];
	double av2[3][Nblocks];
	double sum_prog[3][Nblocks];
	double su2_prog[3][Nblocks];
	double err_prog[3][Nblocks];
	double y;
	
	//calcolo gli integrali Irand e Iprob della funzione data generando numeri casuali seondo le due distribuzioni random e d(x)=3/2(1-x^2) (d_prob)

	
	for (int i=0; i<Nblocks; i++){
		cout<<"Block "<<i+1<<" of "<<Nblocks<<endl;
		sumrand = 0;
//		sumprob = 0;
		sumlin = 0;
		for (int j=0; j<Mthrows; j++){
		
			sumrand = sumrand + cosPi(rnd.Rannyu());
			
//			double x = rnd.d_prob();	
//			sumprob = sumprob + cosPi(x)/(1.5*(1-pow(x,2))); 
			
			double y = rnd.d_lin();
			sumlin = sumlin + cosPi(y)/(2*(1-y)); 
			
			}
		Irand = sumrand/Mthrows;
//		Iprob = sumprob/Mthrows;
		Ilin = sumlin/Mthrows;
		
		ave[0][i] = Irand;
		av2[0][i] = pow(Irand,2);
//		ave[1][i] = Iprob;
//		av2[1][i] = pow(Iprob,2);
		ave[2][i] = Ilin;
		av2[2][i] = pow(Ilin,2);

		sum_prog[0][i] = 0;
		su2_prog[0][i] = 0;
//		sum_prog[1][i] = 0;
//		su2_prog[1][i] = 0;
		sum_prog[2][i] = 0;
		su2_prog[2][i] = 0;
		}
		
	//calcolo le incertezze statistiche 
	
	for (int i=0; i<Nblocks; i++){
		for (int j=0; j<i+1; j++){
			sum_prog[0][i] = sum_prog[0][i] + ave[0][j];
			su2_prog[0][i] = su2_prog[0][i] + av2[0][j];
//			sum_prog[1][i] = sum_prog[1][i] + ave[1][j];
//			su2_prog[1][i] = su2_prog[1][i] + av2[1][j];
			sum_prog[2][i] = sum_prog[2][i] + ave[2][j];
			su2_prog[2][i] = su2_prog[2][i] + av2[2][j];
		}
		sum_prog[0][i] = sum_prog[0][i]/(i+1.);
		su2_prog[0][i] = su2_prog[0][i]/(i+1.);
		err_prog[0][i]=err(su2_prog[0],sum_prog[0],i);
//		sum_prog[1][i] = sum_prog[1][i]/(i+1.);
//		su2_prog[1][i] = su2_prog[1][i]/(i+1.);
//		err_prog[1][i]=err(su2_prog[1],sum_prog[1],i);
		sum_prog[2][i] = sum_prog[2][i]/(i+1.);
		su2_prog[2][i] = su2_prog[2][i]/(i+1.);
		err_prog[2][i]=err(su2_prog[2],sum_prog[2],i);

		dist_rand<<sum_prog[0][i]<<"		"<<err_prog[0][i]<<endl;
//		dist_prob<<sum_prog[1][i]<<"		"<<err_prog[1][i]<<endl;
		dist_lin<<sum_prog[2][i]<<"		"<<err_prog[2][i]<<endl;

		}
		
	dist_rand.close();	
				
return 0;
}
