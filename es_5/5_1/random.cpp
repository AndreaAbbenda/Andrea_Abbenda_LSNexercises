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

Point :: Point(){
	_x=0;
	_y=0;
	_z=0;
	}

Point :: Point(double x, double y, double z){
	_x=x;
	_y=y;
	_z=z;
	}

Point :: Point(Point p, double mean, double sigma, Random rnd){
	_x= rnd.Gauss(p.GetX() + mean, sigma);
	_y= rnd.Gauss(p.GetY() + mean, sigma);
	_z= rnd.Gauss(p.GetZ() + mean, sigma);
	}


Point :: Point(Point p, double range, Random rnd){
	_x= p.GetX() + (rnd.Rannyu()*2 - 1)*range;
	_y= p.GetY() + (rnd.Rannyu()*2 - 1)*range;
	_z= p.GetZ() + (rnd.Rannyu()*2 - 1)*range;
	}

	
Point :: ~Point(){}

void Point :: Copy(Point p){
	this->SetX(p.GetX());
	this->SetY(p.GetY());
	this->SetZ(p.GetZ());
	}

void Point :: Print(){
	cout<<"x: "<<this->GetX()<<endl;
	cout<<"y: "<<this->GetY()<<endl;
	cout<<"z: "<<this->GetZ()<<endl<<endl;
	}

double AEW(int n, Point p){
	double x = p.GetX();
	double y = p.GetY();
	double z = p.GetZ();
	if (n == 1){ return pow(exp(-sqrt(x*x + y*y + z*z)/sqrt(M_PI)),2);}
	else if (n == 2){ return pow(sqrt(2/M_PI)/8*sqrt(x*x + y*y + z*z)*exp(-sqrt(x*x + y*y + z*z)/2)*z/sqrt(x*x + y*y + z*z),2);}
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
