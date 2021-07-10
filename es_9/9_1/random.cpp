/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "random.h"

using namespace std;

double err(double *av2, double *ave, int N){
	if (N==0){return 0;}
	else {return sqrt(abs(av2[N]-pow(ave[N],2))/N);}
	}
double Minimum (double a, double b){ 
	if (a > b){return b; 
	}else if(b >= a){return a;
		}
	}

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}





//---------------------------------------------------------------------------------------------------------------------

Chromosome::Chromosome(){
	N=Cities; 
	C=new int[N]; 
	for(int i=0; i<Cities; i++){C[i]=i;}
	
	int seed[4] = {0000, 0000, 0000, 0001};
	int p1 = 2892;
	int p2 = 2707;
	rnd.SetRandom(seed,p1,p2);
}




Chromosome::Chromosome(unsigned int n){
	N=n;
	C=new int[N];
	for(int i=0; i<Cities; i++){
		C[i]=i;
	};
	
	int seed[4] = {0000, 0000, 0000, 0001};
	int p1 = 2892;
	int p2 = 2707;
	rnd.SetRandom(seed,p1,p2);
};


Chromosome::Chromosome(const Chromosome & c){
	N=c.N;
	C=new int[N];
	for (int i=0;i<N;i++) C[i]=c.C[i];
	
	int seed[4] = {0000, 0000, 0000, 0001};
	int p1 = 2892;
	int p2 = 2707;
	rnd.SetRandom(seed,p1,p2);
};


Chromosome::~Chromosome(){
	delete [] C;
};

void Chromosome::Copy(const Chromosome c){
	for (int i=0;i<N;i++) C[i]=c.C[i];
	
	int seed[4] = {0000, 0000, 0000, 0001};
	int p1 = 2892;
	int p2 = 2707;
	rnd.SetRandom(seed,p1,p2);
};


void Chromosome :: SetComponent(int a, int b){
	C[a]=b;
};

int Chromosome :: GetComponent(int a){
	return C[a];
};

int Chromosome :: GetN(){
	return N;
};

void Chromosome :: Switch(){
	double x = abs(rnd.Rannyu());
	int a = (double)(N-2) * x + 1;
	int b = C[a];
	C[a]=C[a+1];
	C[a+1]=b;
};

void Chromosome :: Shift(){
	int lotto = (double)(N-2) * abs(rnd.Rannyu()) + 1;
	int start = (N-1-lotto) * abs(rnd.Rannyu()) + 1;
	int shift = (N-1-start-lotto)* abs(rnd.Rannyu()) + 1;
	for (int i=0; i<lotto; i++){
		for (int j=0; j<shift; j++){
			int appo = C[start +lotto-  i + j -1];
			C[start +lotto- i + j-1] = C[start + lotto-i + j ];
			C[start + lotto-i + j ] = appo;
			if(start + lotto-i + j  ==N){cout<<"err"<<endl; break;};
			}
		}
};

void Chromosome :: Permute(){
//	cout<<abs(rnd.Rannyu())<<endl;
	int lotto = double((N-1)/2) * abs(rnd.Rannyu()) + 1;
	int start = (N-2*lotto) * abs(rnd.Rannyu()) + 1;
	for (int i=0; i<lotto; i++){
		int appo = C[start + i];
		C[start + i] = C[start + i + lotto];
		C[start + i + lotto] = appo;
		if(start + i + lotto==N){cout<<"err"<<endl; break;};
	}
};

void Chromosome :: Invert(){
	int lotto = (double)((N-1)/2  ) * abs(rnd.Rannyu()) + 1;
	int start = (double)(N-2*lotto) * abs(rnd.Rannyu()) + 1;
	for (int i=0; i<lotto; i++){
		int appo = C[start + i];
		C[start + i] = C[start + 2*lotto -i -1];
		C[start + 2*lotto -i -1] = appo;
		if(start + 2*lotto -i -1 >=N){cout<<"err"<<endl; break;};
		if(start + i <=0){cout<<"err"<<endl; break;};
	}
};


void Chromosome :: Print(){
	for (int i=0; i<N; i++){
		cout<<C[i]<<"  ";
		}
	cout<<endl;
};

int Chromosome :: CheckChromo(){
	int sum =0;
	for (int i=0; i<N; i++){
		sum += C[i];
	}
	int check;
	if (sum == 496 and C[0] == 0){check =0;}
	else {
		check = 1;
		if (sum != 496){ cout<<"ERR: somma indici errata: "<<sum<< endl;}
		if (C[0] != 0){ cout<<"ERR: primo indice errato: "<<C[0]<<endl;}
	}
	return check;
};

















//-------------------------------------------------------------------------------------------------------


Population :: Population(Population &POP){
	N=POP.N;
	P=new Chromosome[N];
	for(int i=0; i<N; i++){P[i].Copy(POP.GetChromosome(i));}
};

Population::Population(unsigned int n){
	N=n;
	P=new Chromosome[N];
};

Population::Population(unsigned int n, Random &r){
	N=n;
	P=new Chromosome[N];
	rnd = r;
};

Population::~Population(){
	delete [] P;
};

void Population :: SetChromosome(Chromosome c, int a){
	P[a].Copy(c);
};

Chromosome Population :: GetChromosome(int a){
	return P[a];
};

int Population :: GetN(){
	return N;
};





void Population :: Start(int shape){              // riempie la matrice Positions di Population con le posizioni da Position.dat
                                                  // inizializza gli N cromosomi
	ifstream Pos;
	if (shape == 0){Pos.open("Positions_circle.dat");
	}else { Pos.open("Positions_square.dat");}
	if (Pos.is_open()){
		for (int i=0; i<Cities; i++){
			Pos >> Positions[0][i] >> Positions[1][i];
		}	
	} else cerr << "PROBLEM: Unable to open Positions" << endl;
	Pos.close();
	
	Chromosome alpha(Cities);
	for (int i=0; i<N; i++){alpha.Switch(); alpha.Shift(); alpha.Invert(); alpha.Permute();}
	for (int i=0; i<N; i++){
		P[i].Copy(alpha);	
		alpha.Switch(); alpha.Shift(); alpha.Invert(); alpha.Permute();
	};
};

void Population :: Order(){
	double fit[N];
	for (int i=0; i<N; i++){
		fit[i] = Fitness(i,1);
	}
	Chromosome appo(Cities);
	double appo2;

	for (int i=0; i<N; i++){		
		double val=fit[i];
		int support = i;
		for(int j=i; j<N; j++){
			if (fit[j] <= val){
				val = fit[j];
				support = j;
			};
		};
		
		appo.Copy(P[i]);
		P[i].Copy(P[support]);
		P[support].Copy(appo);
		
		appo2=fit[i];
		fit[i] = fit[support];
		fit[support]=appo2;
	}		
};

double Population :: Fitness(int chromo, int L){
	double sum = 0;
	double x,y;
	int comp,succ;
	if (L==2){
		for (int i=0; i<Cities-1; i++){
			comp= P[chromo].GetComponent(i);
			succ= P[chromo].GetComponent(i+1);
			x = pow(Positions[0][comp] - Positions[0][succ],2);
			y = pow(Positions[1][comp] - Positions[1][succ],2);
			sum += x; sum +=y;
		}
		comp= P[chromo].GetComponent(Cities-1);
		succ= P[chromo].GetComponent(0);
		x = pow(Positions[0][comp] - Positions[0][succ],2);
		y = pow(Positions[1][comp] - Positions[1][succ],2);
		sum += x; sum +=y;
	}else if(L==1){
		for (int i=0; i<Cities-1; i++){
			comp= P[chromo].GetComponent(i);
			succ= P[chromo].GetComponent(i+1);
			x = pow(Positions[0][comp] - Positions[0][succ],2);
			y = pow(Positions[1][comp] - Positions[1][succ],2);
			sum += sqrt(x+y);
		}
		comp= P[chromo].GetComponent(Cities-1);
		succ= P[chromo].GetComponent(0);
		x = pow(Positions[0][comp] - Positions[0][succ],2);
		y = pow(Positions[1][comp] - Positions[1][succ],2);
		sum += sqrt(x+y);
	}
	return sum;
};

void Population :: Check(){
//	cout<< "Checking ...  ";
	for (int i=0; i<N; i++){
		int check = P[i].CheckChromo();
		if (check == 1){cout<< "ERRORE: problema cromosoma "<<i<<endl; break;}
	};
//	cout<<" OK"<<endl;
	
};

Chromosome Population :: Select(int norm){

    double trigger = abs(rnd.Rannyu());
    int index = 0;
    if (trigger > 0.5){

	double sum =0;
	double prob[N];
	for (int i=0; i<N; i++){
		prob[i] = 1.0/Fitness(i,norm);
		sum += prob[i];		
	}	
	
	for (int i=0; i<N; i++){
		prob[i] = prob[i]/sum;
		}
	
	double test = abs(rnd.Rannyu());
	test = test - prob[index];
	while (test > 0){
		index ++;
		if (index >= N){cout<<"SELECTION ERROR"<<endl; break;}
		test = test - prob[index];
	}
	
    }else if(trigger <= 0.5){
	
	double test = abs(rnd.Rannyu());
	index = int(N * pow(test,2));
    }
    
    Chromosome selected(P[index].GetN());
    selected.Copy(P[index]);
    
    return selected;
}

void Population::Print(){
	for (int i=0; i<N; i++){
		P[i].Print();
	}
	cout<<endl;
}

double Population :: GetRandom(){
	return rnd.Rannyu();
}

void Population::PrintBest(int shape, int norm){
	ofstream out;
	if (shape == 0){
		if(norm == 1){
			out.open("best1_circle.dat",ios::app);
		}else if(norm == 2){
			out.open("best2_circle.dat",ios::app);
		}
	} else {
		if(norm == 1){
			out.open("best1_square.dat",ios::app);
		}else if(norm == 2){
			out.open("best2_square.dat",ios::app);
		}
	}
	
	for (int i=0; i<Cities; i++){
		out<< Positions[0][P[0].GetComponent(i)] << "          "<< Positions[1][P[0].GetComponent(i)]<< endl;
	}
	
	out.close();
}

void Population::Accumulate(int shape, int norm){
	
	ofstream L, AveL;
	if (shape == 0){
		if(norm == 1){
			L.open("L1_circle.dat",ios::app);
			AveL.open("AveL1_circle.dat",ios::app);
		}else if(norm == 2){
			L.open("L2_circle.dat",ios::app);
			AveL.open("AveL2_circle.dat",ios::app);
		}
	} else {
		if(norm == 1){
			L.open("L1_square.dat",ios::app);
			AveL.open("AveL1_square.dat",ios::app);
		}else if(norm == 2){
			L.open("L2_square.dat",ios::app);
			AveL.open("AveL2_square.dat",ios::app);
		}
	}
	
	L << Fitness(0,1)<<"          "<<Fitness(0,2)<<endl;
	
	int n = N/2;
	double sum =0;
	double sum2 =0;
	for (int i=0; i<n; i++){
		double a = Fitness(i,1);
		sum += a;
		sum2 += a*a;
	}
	double ave = sum/(double)n;
	double err = sqrt((sum2/(double)n-ave*ave)/(double)n);
	AveL << ave <<"          "<<err<<endl;
	
	
	L.close(); AveL.close();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
