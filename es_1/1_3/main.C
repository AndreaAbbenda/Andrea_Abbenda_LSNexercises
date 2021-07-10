#include "random.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

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



	int Nblocks=100;
	int Mthrows=1000;
	int Nhit=0;
	int Nthr=0;
	int d=10;
	int L=7;
	double x1, x2, x3 =0;
	double y1, y2, y3 =0;
	double pi=0;
	double ave[Nblocks];
	double sum_prog[Nblocks];
	double su2_prog[Nblocks];
	double err_prog[Nblocks];
	
	ofstream PI("es1_3_pi.txt");
	TApplication app("app", &argc,argv);
	
	TH1D *angle= new TH1D("a","b",200,-M_PI,M_PI);
	TH1D *coso= new TH1D("a","b",200,-1,1);
	TGraph *scatter= new TGraph();
	

	for (int j=0; j<Nblocks; j++){	
		cout<<j+1<<" di "<<Nblocks<<endl;
		for (int i=0; i<Mthrows; i++){
		
//considero una porzione di piano "A" di dimensione D*D, con D=100*d (d=spacing tra 2 linee). Genero il primo punto P1(x1,y1) casuale all'interno di A: esso sarà il primo estremo del mio 'ago'
			
			x1=100*d*(2*rnd.Rannyu()-1);
			y1=100*d*(2*rnd.Rannyu()-1);
	
//genero un secondo punto P2(x2,y2) all'interno di una circonferenza di raggio unitario centrata nell'origine. costruisco un terzo punto P3(x3,y3) a distanza L da P1, la cui direzione è pari a quella di P2 rispetto a O
			
			x2=2*rnd.Rannyu()-1;
			y2=2*rnd.Rannyu()-1;
			double distance=sqrt(pow(x2,2)+pow(y2,2));
	
			while (distance>1){
				x2=2*rnd.Rannyu()-1;
				y2=2*rnd.Rannyu()-1;
				distance=sqrt(pow(x2,2)+pow(y2,2));
				}	
			scatter->SetPoint(scatter->GetN(), x2,y2);	
			double Costheta=x2/distance;
			double Sintheta=y2/distance;
			x3=x1+L*Costheta;
			y3=y1+L*Sintheta;
			angle->Fill(acos(Costheta));	
			coso->Fill(x2);	
			
		
//verifico se il segmento P1P2 interseca la griglia sul piano orizzontale. Le linee sono prese ortogonali all'asse delle ascisse
			Nthr=Nthr+1;
			if (int (x3/d) != int (x1/d)){Nhit=Nhit+1;}
			else if (x3*x1 < 0){Nhit=Nhit+1;}
			}
			
//alla fine della j-esima simulazione, calcolo valor medio ed errore progressivo
			
		if (Nhit==0){pi=2.*L*Nthr/d;
		}else{pi=2.*L*Nthr/(Nhit*d);}
		
		Nhit = 0;
		Nthr = 0;
				
		ave[j]=pi;
		if (j==0){
			sum_prog[j]=pi;
			su2_prog[j]=pow(pi,2);
		}else {
			sum_prog[j]=(sum_prog[j-1]*j+pi)/(j+1.);
			su2_prog[j]=(su2_prog[j-1]*j+pow(pi,2))/(j+1.);
		}
		ave[j]=sum_prog[j];
		err_prog[j]=err(su2_prog,sum_prog,j);
		PI<<ave[j]<<"		"<<err_prog[j]<<endl;
	}
		
	PI.close();
	cout<<ave[Nblocks-1]<<endl;	
	TCanvas *c1=new TCanvas("c1","plot");
	c1->Divide(2);
	c1->cd(1);
	angle->Draw();
	c1->cd(2);
	coso->Draw();
	app.Run();
return 0;
}
