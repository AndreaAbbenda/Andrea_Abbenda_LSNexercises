#include "random.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;
 
int main(int argc, char** argv){
   
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");   //imposto = a p1 e p2 prendendo i primi 2 numeri dal file Primes
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");  //carico in un vettore di 4 elementi i numeri 0000, 0000, 0000, 0001 e costruisco un oggetto di classe random
   string property;            //utilizzando i quattro numeri e p1 e p2
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

	double matrix[12][10000];
	
	for (int i=0; i<12; i++){
		for (int j=0; j<10000; j++){
			matrix[i][j]=0;
			}
		}
				
	for (int i=1; i<10000; i++){
		matrix[0][i]=rnd.Rannyu();
		matrix[4][i]=rnd.Expo(1);
		matrix[8][i]=rnd.Lorentz(0,0.5);
		for (int j=0; j<2; j++){
			matrix[1][i]+=rnd.Rannyu()/2;
			matrix[5][i]+=rnd.Expo(1)/2;
			matrix[9][i]+=rnd.Lorentz(0,1)/2;
			}
		for (int j=0; j<10; j++){
			matrix[2][i]+=rnd.Rannyu()/10;
			matrix[6][i]+=rnd.Expo(1)/10;
			matrix[10][i]+=rnd.Lorentz(0,1)/10;
			}
		for (int j=0; j<100; j++){
			matrix[3][i]+=rnd.Rannyu()/100;
			matrix[7][i]+=rnd.Expo(1)/100;
			matrix[11][i]+=rnd.Lorentz(0,1)/100;
			}
		}
	
	ofstream std("es1_2_std.txt");
	ofstream lor("es1_2_lor.txt");
	ofstream exp("es1_2_exp.txt");
	
	
	for (int i=0; i<10000;i++){
		std<<matrix[0][i]<<"	"<<matrix[1][i]<<"	"<<matrix[2][i]<<"	"<<matrix[3][i]<<endl;
		exp<<matrix[4][i]<<"	"<<matrix[5][i]<<"	"<<matrix[6][i]<<"	"<<matrix[7][i]<<endl;
		lor<<matrix[8][i]<<"	"<<matrix[9][i]<<"	"<<matrix[10][i]<<"	"<<matrix[11][i]<<endl;
		}

	std.close();
	exp.close();
	lor.close();

	TApplication app("app", &argc,argv);
	
	
	TH1D *standard1= new TH1D("a","b",100,-3,3);
	TH1D *standard2= new TH1D("a","b",100,-3,3);
	TH1D *standard10= new TH1D("a","b",100,-3,3);
	TH1D *standard100= new TH1D("a","b",100,-3,3);
	
	for (int i=0; i<10000; i++){
		standard1->Fill(matrix[4][i]);
		standard2->Fill(matrix[5][i]);
		standard10->Fill(matrix[6][i]);
		standard100->Fill(matrix[7][i]);
		}
		
	TCanvas *c1=new TCanvas("c1","plot");
	c1->Divide(2,2);	
	c1->cd(1);
	standard1->Draw("ALP");
	c1->cd(2);
	standard2->Draw("ALP");
	c1->cd(3);
	standard10->Draw("ALP");
	c1->cd(4);
	standard100->Draw("ALP");

app.Run();
return 0;
}
