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
  
	int nsteps = 100000;
	int blocks = 100;
	int steps = nsteps/blocks;
	
	double range1 = 2.6;
	double range2 = 3.5;
	double sigma1 = 1.4;	
	double sigma2 = 2;
	
	Point Start1 (rnd.Gauss(30,3), rnd.Gauss(30,3), rnd.Gauss(30,3));  // nei file .h e .cpp ho inizializzato una classe punto
	Point Start2 (rnd.Gauss(30,5), rnd.Gauss(30,5), rnd.Gauss(30,5));  // i cui datamembri sono x,y,z
	
	Point Old1(0,0,0), Old2(0,0,0);                         // 'Old' è il punto di partenza dell'algoritmo,
	Old1.Copy(Start1);                                      // preso gaussianamente intorno all'origine
	Old2.Copy(Start2);
	
	Point New1 (0,0,0);                                           
	Point New2 (0,0,0);                                           
	double A1=0, A2=0;
	int sumA1=0, sumA2=0;
		
	double ave_N1[blocks];              //errori
	double ave2_N1[blocks];
	double err_N1[blocks];
	double ave_N2[blocks];
	double ave2_N2[blocks];
	double err_N2[blocks];
	double sum_N1 = 0;
	double sum_N2 = 0;

	ofstream N1("es5_1_N1_unif_30.dat");         // output files per scatterplot
	ofstream N2("es5_1_N2_unif_30.dat");
	ofstream R1("es5_1_R1_unif_30.dat");         // output files per analisi <r>
	ofstream R2("es5_1_R2_unif_30.dat");	
	
	
		
	double test_random = rnd.Rannyu();            // genero un numero casuale tra 0 e 1 e lo confronto con A1, ovvero
	                                              // la probabilità di transizione tra il punto Old ed il nuovo punto proposto New.

	for (int i=1; i <= nsteps; i++){              // se test_random è minore di A1 il punto è accettato, altrimenti vengono 
		while (A1 < test_random){	      // generati un nuovi punti New finché ciò non accade
			Point Test1 (Old1, range1, rnd);      //a partire dal punto Old, genero un nuovo punto New. le coordinate di New
			New1.Copy(Test1); 		      //sono estratte uniformemente attorno alle coordinate di Old
			double Q1 = AEW(1, Test1)/AEW(1, Old1);	     //calcolo Q1 ed A1. AEW (Atomic Electron Wavefunction) prende in input
			A1 = Minimum (1, Q1);                        //un oggetto punto ed  e restituisce il modulo
			test_random = rnd.Rannyu();                  //quadro della funzione d'onda valutato nel punto scelto
			sumA1 ++;
			}
			
		sum_N1 += Old1.Radius();		
		if (i % steps == 0){                                 //calcolo degli errori statistici
			int j = i/steps;
			if (j == 1){
				ave_N1[0] = sum_N1/steps;
				ave2_N1[0] = pow(ave_N1[0],2);
				err_N1[0] = 0;
			}else {
				ave_N1[j-1] = (ave_N1[j-2]*(j-1) + sum_N1/steps)/j;
				ave2_N1[j-1] = (ave2_N1[j-2]*(j-1) + pow(sum_N1/steps,2))/j;
				err_N1[j-1] = sqrt((ave2_N1[j-1]-pow(ave_N1[j-1],2))/(j-1));
				}
			R1<<ave_N1[j-1]<<"	"<<err_N1[j-1]<<endl;
			sum_N1 = 0;
			}				
		A1 = 0;
		N1 << Old1.GetX() << "	" << Old1.GetY() << "	"<< Old1.GetZ() << endl;
		Old1.Copy(New1);
	
//		sigma = 4;
	
		while (A2 < test_random){                       //lo stesso procedimento è ripetuto per entrambe le funzioni d'onda
			Point Test2 (Old2, range2, rnd);
			New2.Copy(Test2);
			double Q2 = AEW(2, Test2)/AEW(2, Old2);
			A2 = Minimum (1, Q2);
			test_random = rnd.Rannyu();
			sumA2 ++;
			}
		sum_N2 += Old2.Radius();		
		if (i % steps == 0){
			int j = i/steps;
			if (j == 1){
				ave_N2[0] = sum_N2/steps;
				ave2_N2[0] = pow(ave_N2[0],2);
				err_N2[0] = 0;
			}else {
				ave_N2[j-1] = (ave_N2[j-2]*(j-1) + sum_N2/steps)/j;
				ave2_N2[j-1] = (ave2_N2[j-2]*(j-1) + pow(sum_N2/steps,2))/j;
				err_N2[j-1] = sqrt((ave2_N2[j-1]-pow(ave_N2[j-1],2))/(j-1));
				}
			R2<<ave_N2[j-1]<<"	"<<err_N2[j-1]<<endl;
			sum_N2 = 0;
			}							
		A2 = 0;
		N2 << Old2.GetX() << "	" << Old2.GetY() << "	"<< Old2.GetZ() << endl;
		Old2.Copy(New2);
		
	}
	cout << "Probability of accepting a step (quantic number n = 1), multivariate uniform transition probability : " << (double)nsteps/sumA1 << endl;
	cout << "Probability of accepting a step (quantic number n = 2): multivariate uniform transition probability " << (double)nsteps/sumA2 << endl<<endl;
	N1.close();
	N2.close();
	R1.close();
	R2.close();

	Old1.Copy(Start1);
	Old2.Copy(Start2);

	N1.open("es5_1_N1_norm_30.dat");         // output files per scatterplot
	N2.open("es5_1_N2_norm_30.dat");
	R1.open("es5_1_R1_norm_30.dat");         // output files per analisi <r>
	R2.open("es5_1_R2_norm_30.dat");	
	
	sumA1=0, sumA2=0;
	test_random = rnd.Rannyu();                   
	                                             
	for (int i=1; i <= nsteps; i++){              
		while (A1 < test_random){	     
			Point Test1 (Old1, 0, sigma1, rnd);   // ripeto lo stesso procedimento di prima variando usando una probabilità 
			New1.Copy(Test1); 		      // di transizione multivariata distribuita normalmente 
			double Q1 = AEW(1, Test1)/AEW(1, Old1);	     
			A1 = Minimum (1, Q1);                      
			test_random = rnd.Rannyu();   
			sumA1 ++;
			}
			
		sum_N1 += Old1.Radius();		
		if (i % steps == 0){                                 
			int j = i/steps;
			if (j == 1){
				ave_N1[0] = sum_N1/steps;
				ave2_N1[0] = pow(ave_N1[0],2);
				err_N1[0] = 0;
			}else {
				ave_N1[j-1] = (ave_N1[j-2]*(j-1) + sum_N1/steps)/j;
				ave2_N1[j-1] = (ave2_N1[j-2]*(j-1) + pow(sum_N1/steps,2))/j;
				err_N1[j-1] = sqrt((ave2_N1[j-1]-pow(ave_N1[j-1],2))/(j-1));
				}
			R1<<ave_N1[j-1]<<"	"<<err_N1[j-1]<<endl;
			sum_N1 = 0;
			}				
		A1 = 0;
		N1 << Old1.GetX() << "	" << Old1.GetY() << "	"<< Old1.GetZ() << endl;
		Old1.Copy(New1);
	
		while (A2 < test_random){                       
			Point Test2 (Old2, 0, sigma2, rnd);
			New2.Copy(Test2);
			double Q2 = AEW(2, Test2)/AEW(2, Old2);
			A2 = Minimum (1, Q2);
			test_random = rnd.Rannyu();
			sumA2 ++;
			}
		sum_N2 += Old2.Radius();		
		if (i % steps == 0){
			int j = i/steps;
			if (j == 1){
				ave_N2[0] = sum_N2/steps;
				ave2_N2[0] = pow(ave_N2[0],2);
				err_N2[0] = 0;
			}else {
				ave_N2[j-1] = (ave_N2[j-2]*(j-1) + sum_N2/steps)/j;
				ave2_N2[j-1] = (ave2_N2[j-2]*(j-1) + pow(sum_N2/steps,2))/j;
				err_N2[j-1] = sqrt((ave2_N2[j-1]-pow(ave_N2[j-1],2))/(j-1));
				}
			R2<<ave_N2[j-1]<<"	"<<err_N2[j-1]<<endl;
			sum_N2 = 0;
			}							
		A2 = 0;
		N2 << Old2.GetX() << "	" << Old2.GetY() << "	"<< Old2.GetZ() << endl;
		Old2.Copy(New2);
		
	}
	cout << "Probability of accepting a step (quantic number n = 1), multivariate normal transition probability : " << (double)nsteps/sumA1 << endl;
	cout << "Probability of accepting a step (quantic number n = 2), multivariate normal transition probability : " << (double)nsteps/sumA2 << endl;
	N1.close();
	N2.close();
	R1.close();
	R2.close();	
		
return 0;
}
