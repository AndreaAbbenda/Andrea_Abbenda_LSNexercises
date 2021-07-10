#include "random.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"

#include <fstream>
#include <iostream>
#include <string>

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

	
	
	ofstream o_elementi("es1_1_mean.txt");

	TApplication app("app", &argc,argv);
	
	int M=100000;
	int N=100;
	int L=M/N;
	
	double r[M];
	double x[N];
	double ave[N];
	double av2[N];
	double sum_prog[N];
	double su2_prog[N];
	double err_prog[N];
	double err_x[N];
	for (int i=0; i<M; i++){r[i]=rnd.Rannyu();}
	for (int i=0; i<N; i++){
		x[i]=i;
		ave[i]=0;
		av2[i]=0;
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
		err_x[i]=0;}
		
	for (int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){sum=sum+r[i*L+j];}
		ave[i]=sum/L;
		av2[i]=pow(ave[i],2);
		}
	
	for (int i=0; i<N; i++){
		for (int j=0; j<i+1; j++){
			sum_prog[i]=sum_prog[i]+ave[j];
			su2_prog[i]=su2_prog[i]+av2[j];}
		sum_prog[i]=sum_prog[i]/(i+1);
		su2_prog[i]=su2_prog[i]/(i+1);
		err_prog[i]=err(su2_prog,sum_prog,i);
		o_elementi<<sum_prog[i]-0.5<<"	"<<err_prog[i]<<endl;
		}
	o_elementi.close(); 
	TGraphErrors *graf1= new TGraphErrors(N,x,sum_prog,err_x,err_prog);
	graf1->SetTitle("mean;#throws;<r>");

	ofstream o_elementi2("es1_1_err.txt");
	
	for (int i=0; i<N; i++){
		ave[i]=0;
		av2[i]=0;
		sum_prog[i]=0;
		su2_prog[i]=0;
		err_prog[i]=0;
		err_x[i]=0;}
	
	for (int i=0; i<N; i++){
		double sum=0;
		for (int j=0; j<L; j++){sum=sum+pow(r[i*L+j]-0.5,2);}
		ave[i]=sum/L;
		av2[i]=pow(ave[i],2);
		}
	
	for (int i=0; i<N; i++){
		for (int j=0; j<i+1; j++){
			sum_prog[i]=sum_prog[i]+ave[j];
			su2_prog[i]=su2_prog[i]+av2[j];}
		sum_prog[i]=sum_prog[i]/(i+1);
		su2_prog[i]=su2_prog[i]/(i+1);
		err_prog[i]=err(su2_prog,sum_prog,i);
		o_elementi2<<sum_prog[i]-1./12.<<"	"<<err_prog[i]<<endl;
		}
	o_elementi2.close(); 
	TGraphErrors *graf2= new TGraphErrors(N,x,sum_prog,err_x,err_prog);
	graf2->SetTitle("error;#throws;<r-0.5>^2");

	ofstream o_elementi3("es1_1_chi.txt");

	double chi[100];
	double Xchi[100];
	for (int i=0; i<100;i++){Xchi[i]=i;}
	int n=10000;
	M=100;
	double s[n];
	int exp[M];
	double chi2=0;
	for (int j=0; j<100; j++){	
		int a=0;
		for (int i=0; i<M; i++){exp[i]=0;}
		for (int i=0; i<n; i++){
			s[i]=rnd.Rannyu();
			a=int(M*s[i]);
			exp[a]++;}
		chi2=0;
		for (int i=0; i<M; i++){
			chi2+=pow(exp[i]-n/M,2)/(n/M);
			}
		chi[j]=chi2;
		o_elementi3<<chi[j]<<endl;
		}
	o_elementi3.close();
	
	TGraph *graf3= new TGraph(100,Xchi,chi);
	graf3->SetTitle("chi^2;#throws;chi^2");
		
	TCanvas *c1=new TCanvas("c1","plot");
	c1->Divide(3);
	
	c1->cd(1);
	graf1->Draw("ALP");
	c1->cd(2);
	graf2->Draw("ALP");
	c1->cd(3);
	graf3->Draw("ALP");

	
	

app.Run();
return 0;
}
