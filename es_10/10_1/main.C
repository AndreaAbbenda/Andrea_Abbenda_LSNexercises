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
  
  

  
  int shape;
  int box_side = 100;
  int n_chromosomes = 100;
  int norm = 1;
  
  int Tmax =7500;               // imposto le temperature per il SA
  int n_t = 15;
  int t_steps = 2000;
  int t;
  
  int sum =0; int acc_att =0;
  
  
  cout << "Select cities' disposition [s/c]"  << endl;
  string s_c;
  cin >> s_c;
  if (s_c == "s"){ shape = 1;}
  else if (s_c == "c"){ shape = 0;}


  CreatePattern(box_side, shape, rnd2);

  Population POP(n_chromosomes,rnd2);
  
  Initialize(POP, shape);  
 
    	
  for (int j=0; j<n_t; j++){
     t = Tmax * (1 - (double)j /(double)n_t);
     	
     for(int i=1; i<=t_steps; i++){
  	
  	if (i % t_steps == 0) {
  	    double sum1=0;                  // calcolo la lunghezza media dei migliori N/2 cammini
  	    int N = POP.GetN()/2;
		for(int i=0; i<N; i++){
		sum1 += POP.Fitness(i,1);
		}

   	    if(shape ==0){cout<<"Ordered disposition on a circumference"<<endl;}
  	    else if(shape ==1){cout<<"Casual distribution in squared box"<<endl;} 	
  	    cout<<"Temperature: "<<t<<endl;
  	    cout<<"Average pattern,L1 norm :"<< sum1/(double)N<<endl;
  	    cout<<"best pattern, L1 norm :" << POP.Fitness(0,1)<<endl; 
  	    cout<<"acceptance rate: "<<(double)sum/n_chromosomes/(double)t_steps<<endl<<endl; 
  	    sum=0;
  	}
  	  	
  	acc_att = Mutate(POP, norm, t);
  	sum += acc_att;
  	
  	POP.Order();
  	POP.Check();
  	POP.Accumulate(shape,norm);
  	
     }
  }
  
  POP.PrintBest(shape,norm);	
   
 
		
return 0;
}








void CreatePattern(int a, int shape, Random rnd){
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


int Mutate(Population &POP, int norm, int t){
	
	Population support(POP);
	int n = POP.GetN()/2;
	Chromosome C1(POP.GetChromosome(0));
	Chromosome C2(POP.GetChromosome(0));
	int index1, index2; int ret=0;
	double x, Boltz1, Boltz2;
	
	for (int i=0; i<n; i++){
		
		index1 = POP.SelectInt(norm);
		index2 = POP.SelectInt(norm); //cout<<index1<<"     "<<index2<<endl;
		C1.Copy(POP.GetChromosome(index1)); 
		C2.Copy(POP.GetChromosome(index2));
		
		double prob =POP.GetRandom();
		if (prob>0.4){Crossover(C1, C2, POP.GetRandom());}
		
		x = POP.GetRandom(); if(x<0.1){C1.Switch();}; 
		x = POP.GetRandom(); if(x<0.1){C1.Shift();}
		x = POP.GetRandom(); if(x<0.1){C1.Permute();}
		x = POP.GetRandom(); if(x<0.1){C1.Invert();}
		
		x = POP.GetRandom(); if(x<0.1){C2.Switch();}
		x = POP.GetRandom(); if(x<0.1){C2.Shift();}
		x = POP.GetRandom(); if(x<0.1){C2.Permute();}
		x = POP.GetRandom(); if(x<0.1){C2.Invert();}
		
		support.SetChromosome(C1, 2*i);
		support.SetChromosome(C2, 2*i+1);		
		
		Boltz1 = Minimum(1, exp(-(support.Fitness(2*i, norm) - POP.Fitness(index1, norm))/(double)t));
		Boltz2 = Minimum(1, exp(-(support.Fitness(2*i+1, norm) - POP.Fitness(index2, norm))/(double)t));
		
//		cout<<Boltz1<<"        "<<Boltz2<<endl;
		
		if(Boltz1 < POP.GetRandom()){support.SetChromosome(POP.GetChromosome(index1), 2*i);}else{ret ++;}
		if(Boltz2 < POP.GetRandom()){support.SetChromosome(POP.GetChromosome(index2), 2*i+1);}else{ret ++;}		
	}
	
	for (int i=0; i<POP.GetN(); i++) POP.SetChromosome(support.GetChromosome(i), i);
	
	return ret;
	
}

	
	





























