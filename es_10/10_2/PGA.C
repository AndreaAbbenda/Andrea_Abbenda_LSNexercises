#include "random.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>

#include "mpi.h"

using namespace std;


int main(int argc, char** argv){
  
  
   
  int size, rank;
  int shape = 1;



   Random rnd2;                                     //inizializzo generatore di numeri casuali
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
  
  
  
  int box_side = 100;              //definisco i parametri comuni di ogni core
  int n_chromosomes = 100;
  int n_mutations = 30000;
  int norm = 1;
  
  CreatePattern(box_side, shape, rnd2);     //creo la disposizione delle città
  
  
  
  MPI_Init(&argc,&argv);                    //parallelizzo
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat;
  int itag=1;



   int p[2][4];                                      //imposto il seed dell'oggetto random a seconda del core
   Primes.open("Primes");   
   if (Primes.is_open()){
      for (int i=0; i<4; i++) {Primes >> p[0][i] >> p[1][i];}
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

  rnd2.SetRandom(seed,p[0][rank],p[1][rank]);



  Population POP(n_chromosomes, rnd2);           
                                                   // genero ed inizializzo la Population
  Initialize(POP, shape);   

  
  int support[Cities];

  for(int i=1; i<=n_mutations; i++){
  	
  	Mutate(POP, norm);
  	
  	POP.Order();
  	POP.Check();
  	POP.Accumulate(shape,rank);
  	                                                   //ogni 2000 mutazioni scambio i best pattern
  	if (i % 2000 == 0 and i != n_mutations){
  		if (rank ==0) cout<<"shape "<<shape<<", mutation n° "<<i<<" of "<<n_mutations<<", exchanging best patterns"<<endl;
  		
  		for (int j=0; j<Cities; j++) support[j] = POP.GetChromosome(0).GetComponent(j);   //copio il best pattern in support[]
  		
  		if (rank < size-1) MPI_Send(&support,Cities,MPI_INTEGER,rank+1,itag,MPI_COMM_WORLD);  //invio support[]
		else if (rank == size-1) MPI_Send(&support,Cities,MPI_INTEGER,0,itag,MPI_COMM_WORLD);
	
		if (rank > 0) MPI_Recv(&support,Cities,MPI_INTEGER,rank-1,itag,MPI_COMM_WORLD,&stat);  //ricevo support[]
		else if(rank == 0) MPI_Recv(&support,Cities,MPI_INTEGER,size-1,itag,MPI_COMM_WORLD,&stat);
		
		Chromosome Csupp(support, Cities);      //creo un Chromosome a partire da support[] ricevuto e lo imposto in Population
		POP.SetChromosome(Csupp, 0);
	}  	
  }


  
  double best_norm[size];  for(int i=0; i<size; i++) best_norm[i]=0;
  int best_path[size][Cities];
  
  double send_norm = POP.Fitness(0, norm);
  int send_path[Cities]; for(int i=0; i<Cities; i++) send_path[i] = POP.GetChromosome(0).GetComponent(i);
  
  int index;
  
  MPI_Gather(&send_norm,1,MPI_REAL8,best_norm,1,MPI_REAL8,0,MPI_COMM_WORLD);                 //raccolgo la fitness dei migliori cammini
  MPI_Gather(&send_path,Cities,MPI_INTEGER,best_path,Cities,MPI_INTEGER,0,MPI_COMM_WORLD);   //raccolgo i migliori cammini
  
  if (rank == 0){ 
  	index = MinIndex(best_norm, size); cout<<endl<<"index of best path: "<<index<<endl;   //trovo il valore minore
  	for(int i=0; i<size; i++) cout<<best_norm[i]<<endl;
  	
  	Chromosome Csupp(best_path[index], Cities);
	POP.SetChromosome(Csupp, 0);
  	POP.PrintBest(shape,norm);                     //stampo il cammino migliore
	
	Filter(shape, index);       		
  	
  }
  

  MPI_Finalize();	


		
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

	
int MinIndex(double vals[], int len){
	int index = 0;
	double ref = vals[0];
	for (int i=0; i<len; i++) {
		if (vals[i] < ref) {
			index = i;
			ref = vals[i];
		}
	}
	return index;
}	

void Filter(int shape, int index){

  	ifstream L_all, AveL_all;
  	ofstream L_best, AveL_best;
  	
	if (shape == 0){
		L_all.open("L1_circle_all.dat");
		AveL_all.open("AveL1_circle_all.dat");
		L_best.open("L1_circle_best.dat");
		AveL_best.open("AveL1_circle_best.dat");
	} else {
		L_all.open("L1_square_all.dat");
		AveL_all.open("AveL1_square_all.dat");
		L_best.open("L1_square_best.dat");
		AveL_best.open("AveL1_square_best.dat");
	}
	
	while (!L_all.eof()){
        	int x;
        	double appo1, appo2;
        	L_all >> x >> appo1 >> appo2;
        	if (x == index){
        		L_best<< appo1 << "         " << appo2 << endl;
        	}
        }

	while (!AveL_all.eof()){
        	int x;
        	double appo1, appo2;
        	AveL_all >> x >> appo1 >> appo2;
        	if (x == index){
        		AveL_best<< appo1 << "         " << appo2 << endl;
        	}
        }
        
        L_all.close(); AveL_all.close();
  	L_best.close(); AveL_best.close();

}




























