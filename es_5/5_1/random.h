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

class Point {

private:
  double _x, _y, _z;

protected:

public:
  // constructors
  Point();
  Point(double x, double y, double z);
  Point(Point p, double range, Random rnd);
  Point(Point p, double mean, double sigma, Random rnd);
  // destructor
  ~Point();
  // methods
  double SetX(double a){_x=a;};
  double SetY(double b){_y=b;};
  double SetZ(double c){_z=c;};
  double GetX(){double a = _x; return a;};
  double GetY(){double b = _y; return b;};
  double GetZ(){double c = _z; return c;};
  double Radius(){ return sqrt(_x*_x + _y*_y + _z*_z);};
  void Copy(Point p);
  void Print();

};

double AEW(int n, Point p);

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
