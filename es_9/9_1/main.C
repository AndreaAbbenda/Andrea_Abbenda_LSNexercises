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
   
   Random rnd2;
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
            rnd2.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  

  
  int shape;                                  // definisco i parametri iniziali
  int box_side = 100;
  int n_chromosomes = 100;
  int n_mutations = 30000;
  int norm = 1;
  

  cout << "Select cities' disposition [s/c]"  << endl;          //seleziono la disposizione delle città
  string s_c;
  cin >> s_c;
  if (s_c == "s"){ shape = 1;}
  else if (s_c == "c"){ shape = 0;}


  CreatePattern(box_side, shape, rnd2);
  
  Population POP(n_chromosomes, rnd2);
  
  Initialize(POP, shape);  
  	
  	
  for(int i=1; i<=n_mutations; i++){
  	
  	if (i % 500 == 0) {cout<<"Mutation number: "<<i<<endl<<"best pattern, L1 norm :" << POP.Fitness(0,1)<<endl; cout<<endl; }
  	  	
  	Mutate(POP, norm);
  	
  	POP.Order();
  	POP.Check();
  	POP.Accumulate(shape,norm);
  	
  }
  
  POP.PrintBest(shape,norm);	
   
  /*
  ofstream multi1, multi2;
  if (shape ==0){
  	multi1.open("progressive1_circle.dat");
  	multi2.open("progressive2_circle.dat");
  }else if (shape==1){
  	multi1.open("progressive1_square.dat");
  	multi2.open("progressive2_square.dat");
  }
  
  n_mutations = 10000; 
  
  for (int i=1; i<=2; i++){
  norm = i;
  cout<< "norma L^"<<i<<endl<<endl;
  
  	for (int j=1; j<= 15; j++){
  		n_chromosomes = j*10;
  		cout<<"n chromosomes: "<<n_chromosomes<<endl;
  		Population POP2(n_chromosomes,rnd2);
  		Initialize(POP2, shape); 
  		
  		for(int k=1; k<=n_mutations; k++){
  		  	Mutate(POP2, norm);
   		    	POP2.Order();
  			POP2.Check();
  		}
  		double a=POP2.Fitness(0,1);
  		cout<<a<<endl;
  		if (i==1){multi1<<a<<endl;} else if(i==2){multi2<<a<<endl;}
  	}
  }
  
  multi1.close(); multi2.close();
			
*/

		
return 0;
}








void CreatePattern(int a, int shape, Random &rnd){
	if(shape == 0){
		ofstream Pos("Positions_circle.dat");
		for (int i=0; i<Cities; i++){
			double x = a*cos(2*M_PI/Cities *i);
			double y = a*sin(2*M_PI/Cities *i);
			Pos<< x<<"          "<<y<<endl;
			}
		Pos.close();
	}else {
		ofstream Pos("Positions_square.dat");
		for (int i=0; i<Cities; i++){
			double x = a*rnd.Rannyu();
			double y = a*rnd.Rannyu();
			Pos<< x<<"          "<<y<<endl;
			}
		Pos.close();
	}		
}

void Initialize(Population &POP, int shape){
	POP.Start(shape);
	POP.Order();
	POP.Check();	
}

void appo(Chromosome &c, int a){
	c.SetComponent(a,100);
	}

void Crossover(Chromosome &x1, Chromosome &x2, double random){

	int N = x1.GetN();
	int a = double(N-1) * random + 1;
	
	Chromosome save1(N);
	Chromosome save2(N);
	
	save1.Copy(x1);
	save2.Copy(x2);
	
	int reference =0;
	
	for (int i=a; i<N; i++){
		
		int trigger1 = 0;
		int component = i-a;
		while (trigger1 ==0){
			
			int trovata = 0;			
			reference = save2.GetComponent(component);
			
			for (int j=0; j<i; j++){
				if (x1.GetComponent(j)==reference){trovata = 1;  break;}					
			}

			if (trovata ==0){
				x1.SetComponent(i, reference);
				trigger1 = 1;
			}else if(trovata == 1){
				component++;
			}		
		}
				
		trigger1 = 0;
		component = i-a;
		while (trigger1 ==0){
			
			int trovata = 0;			
			reference = save1.GetComponent(component);
			
			for (int j=0; j<i; j++){
				if (x2.GetComponent(j)==reference){trovata = 1;  break;}					
			}

			if (trovata ==0){
				x2.SetComponent(i, reference);
				trigger1 = 1;
			}else if(trovata == 1){
				component++;
			}		
		}
	}
}


void Mutate(Population &POP, int norm){
	
	Population support(POP);
	int n = POP.GetN()/2;
	Chromosome C1(POP.GetChromosome(0));
	Chromosome C2(POP.GetChromosome(0));
	double x;
	
	for (int i=0; i<n; i++){
		
		C1.Copy(POP.Select(norm));
		C2.Copy(POP.Select(norm));
		double prob =POP.GetRandom();
		if (prob>0.4){Crossover(C1, C2, POP.GetRandom());}
		
		x = POP.GetRandom(); if(x<0.1){C1.Switch();}
		x = POP.GetRandom(); if(x<0.1){C1.Shift();}
		x = POP.GetRandom(); if(x<0.1){C1.Permute();}
		x = POP.GetRandom(); if(x<0.1){C1.Invert();}
		
		x = POP.GetRandom(); if(x<0.1){C2.Switch();}
		x = POP.GetRandom(); if(x<0.1){C2.Shift();}
		x = POP.GetRandom(); if(x<0.1){C2.Permute();}
		x = POP.GetRandom(); if(x<0.1){C2.Invert();}
		
		support.SetChromosome(C1, 2*i);
		support.SetChromosome(C2, 2*i+1);
	}
	
	for (int i=0; i<POP.GetN(); i++) POP.SetChromosome(support.GetChromosome(i), i);
}

	
	





























