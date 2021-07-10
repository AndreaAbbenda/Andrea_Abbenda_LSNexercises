#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

//funzione signum(x): restituisce il segno di x

int signum(double x){
	if (x >= 0){return 1;} else if(x < 0){return -1;}
	}
	
int main(int argc, char** argv){
   
   Random rnd;
   int seed[7];
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

	ofstream discrete("es2_2_discrete.txt");
	ofstream continuum("es2_2_continuum.txt");
	
	
	int Mblocks = 10000;
	int Nsteps = 100;
	double ave[2][Nsteps];
	double av2[2][Nsteps];
	for (int i=0; i<Nsteps; i++){ave[0][i]=0; av2[0][i]=0; ave[1][i]=0; av2[1][i]=0;}

	double distance = 0;
	double L = 1;   //lunghezza dello step 
	
	//i random walk sono implementati come matrici 3XNsteps. partono sempre da 0,0,0
	
	double rw_disc[3][Nsteps+1];	rw_disc[0][0]=0; rw_disc[1][0]=0; rw_disc[2][0]=0;
	double rw_cont[3][Nsteps+1];	rw_cont[0][0]=0; rw_cont[1][0]=0; rw_cont[2][0]=0;
	
	for (int j=0; j<Mblocks; j++){	
		for (int i=1; i<=Nsteps; i++){
		
			//discrete random walk
		
			int dim = 3*rnd.Rannyu();  //scelgo la dimensione in cui compiere il psso
			int direction = signum(rnd.Rannyu() - 0.5);  //scelgo la direzione del passo
			for (int d=0; d<3; d++){
				if (d==dim){rw_disc[dim][i] = rw_disc[dim][i-1] + L*direction;}
				else {rw_disc[d][i] = rw_disc[d][i-1];}
				}
			distance = sqrt(pow(rw_disc[0][i],2)+pow(rw_disc[1][i],2)+pow(rw_disc[2][i],2));
			ave[0][i-1] = ave[0][i-1] + distance/Mblocks;
			av2[0][i-1] = av2[0][i-1] + pow(distance,2)/Mblocks;
				
			//continue random walk
			
			double Theta = rnd.Rannyu()*2*M_PI;  //genero l'angolo theta sul piano XY, theta in [0,2*M_PI]
			double Phi = rnd.Rannyu()*M_PI;   //genero l'angolo phi rispetto all'asse Z, phi in [0,M_PI]
			rw_cont[0][i] = rw_cont[0][i-1] + L*cos(Theta)*sin(Phi);
			rw_cont[1][i] = rw_cont[1][i-1] + L*sin(Theta)*sin(Phi);
			rw_cont[2][i] = rw_cont[2][i-1] + L*cos(Phi);
			distance = sqrt(pow(rw_cont[0][i],2)+pow(rw_cont[1][i],2)+pow(rw_cont[2][i],2));
			ave[1][i-1] = ave[1][i-1] + distance/Mblocks;
			av2[1][i-1] = av2[1][i-1] + pow(distance,2)/Mblocks;		
		}
	}
	
	// calcolo media e deviazione standard al variare del passo
	
	double Raverage = 0;
	double Rerror = 0;
	for (int i=0; i<Nsteps; i++){
		Raverage = ave[0][i];
		Rerror = sqrt(abs(av2[0][i]-pow(ave[0][i],2))/Nsteps);
		discrete<<Raverage<<"	"<<Rerror<<endl;
		
		Raverage = ave[1][i];
		Rerror = sqrt(abs(av2[1][i]-pow(ave[1][i],2))/Nsteps);
		continuum<<Raverage<<"	"<<Rerror<<endl;
	}
	
	discrete.close();
	continuum.close();
	
return 0;
}
