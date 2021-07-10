/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

#include <string>

using namespace std;

int main(){ 
  com=0;
  ifstream test("old.final");
  if (test.is_open()){
    cout << "Continue previous simulation? [y/n]"  << endl;
    string yes_no;
    cin >> yes_no;
    if (yes_no == "y"){            //Inizialization
       Input("old.final");                               // ho sovrascritto la funzione Input() affinchè il sistema possa ripartire
       }else if (yes_no == "n"){                         // con le ultime due configurazioni della simulazione precedente. 
       Input();                                          // Input(filename) prende in input la configurazione precedente a 
       };                                                // config.final
  }else{
    Input();
  }
  test.close();
          
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    com=0;
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
      if (iblk == nblk and istep == nstep-1){
     	PrintConf("old.final");
      }

    }
    cout<<iblk<<"        "<<com<<endl;
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
//  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses S.I. units " << endl<<endl;
  cout << "Simulating gas phase " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
  ReadInput.open("input.gas"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho/pow(0.00000000034,3)<< " particles per m^3"<< endl;
  vol = (double)npart/(rho);
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol * pow(0.00000000034,3) << " m^3"<< endl;
  cout << "Edge of the simulation box = " << box * 0.00000000034 << " m"<< endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> Nstep;
  ReadInput >> nblk;
  nstep = (int)(Nstep/nblk);
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta * 0.34/1000000000*sqrt(1.38*1.66/1200000)<<" sec"<< endl;
  cout << "Number of total steps = " << Nstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps per block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void Input(const std::string& filename){    //Prepare all stuff for the simulation using a configuration from a given file.dat
  ifstream ReadInput,ReadConf;
//  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses S.I. units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
  ReadInput.open("input.gas"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho/pow(0.00000000034,3)<< " particles per m^3"<< endl;
  vol = (double)npart/(rho);
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol * pow(0.00000000034,3) << " m^3"<< endl;
  cout << "Edge of the simulation box = " << box * 0.00000000034 << " m"<< endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> Nstep;
  ReadInput >> nblk;
  nstep = (int)(Nstep/nblk);
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta * 0.34/1000000000*sqrt(1.38*1.66/1200000)<<" sec"<< endl;
  cout << "Number of total steps = " << Nstep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps per block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  
  //measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

  //Initialising penultimate configuration from old.final
  cout << "Read old configuration from file "<<filename<< endl << endl;
  ReadConf.open(filename);
  for (int i=0; i<npart; ++i){
    ReadConf >> xorigin[i] >> yorigin[i] >> zorigin[i];
    xorigin[i] = xorigin[i] * box;
    yorigin[i] = yorigin[i] * box;
    zorigin[i] = zorigin[i] * box;
    xold[i] = xorigin[i]; yold[i] = yorigin[i]; zold[i] = zorigin[i]; 
  }
  ReadConf.close();

//Initialising initial configuration from config.final
  cout << "Read initial configuration from file config.final"<< endl << endl;
  ReadConf.open("config.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  
  double fs = 0;
  int j = 0;
  
  while (fs < 0.99 or fs > 1.01){              // al fine di termalizzare il sistema alla giusta temperatura, l'inizializzazione
                                               // del sistema è ripetuta più volte finché il fattore di scala della 
                                               // temperatura sta tra 0.95 e 1.05
   Move();                                    
                                              
   double sumv[3] = {0.0, 0.0, 0.0};          
   for (int i=0; i<npart; ++i){               // a partire dalle configurazioni old.final (-1) e cinfig.final (0), muovo le 
     sumv[0] += vx[i];                        // particelle con l'algoritmo di Verlet, ottenendo la configurazione (1). Con
     sumv[1] += vy[i];                        // (-1) e (1) clcolo le velocità in (0), estrapolo il fattore di scala delle 
     sumv[2] += vz[i];                        // velocità e le riscalo secondo la temperatura desiderata
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   cout<<"velocity scale factor after restart number "<< j + 1 << ": "<< fs << endl << endl;
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xorigin[i] = Pbc(xold[i] - vx[i] * delta);
     yorigin[i] = Pbc(yold[i] - vy[i] * delta);
     zorigin[i] = Pbc(zold[i] - vz[i] * delta);
     
     x[i] = xold[i];
     y[i] = yold[i];
     z[i] = zold[i];
     
     xold[i] = xorigin[i];
     yold[i] = yorigin[i];
     zold[i] = zorigin[i];
     }
                                                           // una volta estrapolate le configurazioni (-1_final) e (0), il sistema compie
                                                           // Nstep/10 steps al fine di termalizzarsi, dopodichè la procedura di
     for(int istep=1; istep <= (int)(Nstep/10); ++istep){  // inizializzazione è fatta ripartire fino a che la temperatura del sistema 
     Move();                                               // è pressoché uguale a quella desiderata
     
//     tt = 0;
//     for (int i=0; i<npart; ++i) tt += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
//     thermalgraph << (2.0 / 3.0) * tt/(double)npart << endl;
     }
     
     
   j++;
   }
//   thermalgraph.close();                                                       
                                                        
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
    
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
//  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  
//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
//update of the histogram of g(r)
     walker[igofr + (int)(dr/bin_size)] += 2;
     com += dr;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv] = 120*v/(double)npart; //Potential energy per particle      //modifico le unità di misura affinché i risultati non siano
    walker[ik] = 120*t/(double)npart; //Kinetic energy per particle        //in unità ridotte ma in S.I. nel caso di Argon
    walker[it] = 120*(2.0 / 3.0) * t/(double)npart; //Temperature         //(m = 39.948 amu, sigma = 0.34 nm, epsilon/k = 120 K)
    walker[ie] = 120*(t+v)/(double)npart; //Total energy per particle

    return;
}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
   
//   double r, gdir;
   ofstream Gofr, Epot, Etot, Ekin, Temp;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    
    Epot.open("output.epot.gas",ios::app);
    Etot.open("output.etot.gas",ios::app);
    Ekin.open("output.ekin.gas",ios::app);
    Temp.open("output.temp.gas",ios::app);
    Gofr.open("output.gofr.gas",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
    stima_kin = blk_av[ik]/blk_norm; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_kin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
    
    stima_temp = blk_av[it]/blk_norm; //Potential energy
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    
    for (int k=igofr; k<igofr+nbins; ++k){   //G(r)
      stima_gofr = (blk_av[k]/blk_norm) / (rho * (double)npart * 4 * M_PI * (pow((k-igofr+1)*bin_size,3)-pow((k-igofr)*bin_size,3))/3);
      glob_av[k] += stima_gofr;
      glob_av2[k] += stima_gofr*stima_gofr;
      err_gofr=Error(glob_av[k],glob_av2[k],iblk);
      if (iblk == nblk){
        Gofr << (k - igofr) * bin_size * 0.34<< "		"<< glob_av[k]/(double)iblk << "	"<< err_gofr << endl;
        }
      }
      

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_pot << endl;
//Kinetic energy per particle
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_kin << endl;
//Total energy per particle
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
    
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Epot.close();
    Temp.close();
    Gofr.close();

}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void PrintConf(const std::string& filename){ //Write any chosen configuration
  ofstream WriteConf;

  WriteConf.open(filename);

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
