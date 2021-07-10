/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __NVT__
#define __NVT__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,igofr;
double stima_pot, stima_kin, stima_etot, stima_temp,stima_gofr,bin_size,nbins;
double walker[m_props];
double com;

// averages
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double err_pot,err_kin,err_etot,err_temp,err_gofr;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part],xorigin[m_part],yorigin[m_part],zorigin[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, Nstep, nblk, iprint;
double delta;

// errors
double sumetot = 0, sumepot = 0, sumekin = 0, sumtemp = 0;
double sum2etot = 0, sum2epot = 0, sum2ekin = 0, sum2temp = 0;
//functions
void Input(void);
void Input(const std::string& filename);
void Move(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double,double,int);

void PrintConf(const std::string& filename);

#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
