/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <cmath>

#ifndef __Random__
#define __Random__

double err(double *av2, double *ave, int N);
double Minimum (double a, double b);

const int Cities = 32;

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};

void CreatePattern(int, int, Random);



class Chromosome {

private:

protected:
   unsigned int N;
   int* C;
   Random rnd;

public:
  // constructors
  Chromosome();
  Chromosome(unsigned int);
  Chromosome(const Chromosome&);
  // destructor
  ~Chromosome();
  // methods
  void Copy(const Chromosome c);
  void SetComponent(int a, int b);
  int GetComponent(int);
  int GetN();
  
  void Switch();
  void Shift();
  void Permute();
  void Invert();
 
  void Print();
  int CheckChromo();

};

void Crossover(Chromosome &, Chromosome &, double);
void appo(Chromosome &, int);


class Population {

private:
  
protected:
  Chromosome* P;
  unsigned int N;
  Random rnd;
  double Positions[2][Cities];

public:
  // constructors
  Population(){N=20;};
  Population(Population&);
  Population(unsigned int);
  Population(unsigned int, Random);
  // destructor
  ~Population();
  // methods
  int GetN();
  Chromosome GetChromosome(int);
  void SetChromosome(Chromosome,int);
  void Start(int);
  void Order();
  double Fitness(int,int);
  void Check();
  Chromosome Select(int);
  int SelectInt(int);
  void Print();
  void PrintBest(int,int);
  double GetRandom();
  void Accumulate(int,int);
};

void Initialize(Population&, int);
int Mutate(Population&, int,int);

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
