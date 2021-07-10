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
#include "MolDyn_NVE.h"

#include <string>

using namespace std;

int main(){ 
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
  
// errors
  const int blocks = (int) nstep/100;
  double ave_etot[blocks], ave_epot[blocks], ave_ekin[blocks], ave_temp[blocks];
  double ave2_etot[blocks], ave2_epot[blocks], ave2_ekin[blocks], ave2_temp[blocks];
  double err_etot[blocks], err_epot[blocks], err_ekin[blocks], err_temp[blocks];
  ofstream errEpot, errEkin, errEtot, errTemp;

  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
        
	sumetot += stima_etot;                            //da riga 51 a 103 il codice calcola le incertezze relative alle grandezze 
	sumepot += stima_pot;                             //fisiche studiate con ilmetodo del data-blocking
	sumekin += stima_kin;
	sumtemp += stima_temp;
    
	if(istep%100 == 0){
	   errEpot.open("output_err_epot.dat",ios::app);
  	   errEkin.open("output_err_ekin.dat",ios::app);
  	   errTemp.open("output_err_temp.dat",ios::app);
  	   errEtot.open("output_err_etot.dat",ios::app);
    	
    	   int j = (int)(istep/100);
    	   if (j == 1){
    	     ave_etot[0] = sumetot/10;
    	     ave2_etot[0] = pow((sumetot/10),2);
    	     
    	     ave_epot[0] = sumepot/10;
    	     ave2_epot[0] = pow((sumepot/10),2);
    	     
    	     ave_ekin[0] = sumekin/10;
    	     ave2_ekin[0] = pow((sumekin/10),2);
    	     
    	     ave_temp[0] = sumtemp/10;
    	     ave2_temp[0] = pow((sumtemp/10),2);
    	     
    	   }else {    	     
    	     ave_etot[j - 1] = (ave_etot[j-2] * (j-1) + sumetot/10) / j;
    	     ave2_etot[j - 1] = (ave2_etot[j-2] * (j-1) + pow((sumetot/10),2)) / j;
    	     err_etot[j-1] = sqrt((ave2_etot[j - 1] - pow(ave_etot[j - 1],2))/(j-1));
    	     
    	     ave_epot[j - 1] = (ave_epot[j-2] * (j-1) + sumepot/10) / j;
    	     ave2_epot[j - 1] = (ave2_epot[j-2] * (j-1) + pow((sumepot/10),2)) / j;
    	     err_epot[j-1] = sqrt((ave2_epot[j - 1] - pow(ave_epot[j - 1],2))/(j-1));
    	     
    	     ave_ekin[j - 1] = (ave_ekin[j-2] * (j-1) + sumekin/10) / j;
    	     ave2_ekin[j - 1] = (ave2_ekin[j-2] * (j-1) + pow((sumekin/10),2)) / j;
    	     err_ekin[j-1] = sqrt((ave2_ekin[j - 1] - pow(ave_ekin[j - 1],2))/(j-1));
    	     
    	     ave_temp[j - 1] = (ave_temp[j-2] * (j-1) + sumtemp/10) / j;
    	     ave2_temp[j - 1] = (ave2_temp[j-2] * (j-1) + pow((sumtemp/10),2)) / j;
    	     err_temp[j-1] = sqrt((ave2_temp[j - 1] - pow(ave_temp[j - 1],2))/(j-1));
    	   }
    	   errEtot << ave_etot[j - 1] << "		" << err_etot[j-1] << endl;
    	   errEpot << ave_epot[j - 1] << "		" << err_epot[j-1] << endl;
    	   errEkin << ave_ekin[j - 1] << "		" << err_ekin[j-1] << endl;
    	   errTemp << ave_temp[j - 1] << "		" << err_temp[j-1] << endl;
    	
    	   sumetot = 0, sumepot = 0, sumekin = 0, sumtemp = 0;
    	   errEpot.close();
	   errEkin.close();
	   errTemp.close();
	   errEtot.close();
       }

//        ConfXYZ(nconf);   //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
     
     if (istep == nstep-1){
     	PrintConf("old.final");
     }
  }
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

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
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.solid"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
  
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
  double tt;
  int j = 0;
//  ofstream thermalgraph("es4_1_thermal.dat");
  
  while (fs < 0.95 or fs > 1.05){              // al fine di termalizzare il sistema alla giusta temperatura, l'inizializzazione
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
     int nconf = 1;                                        // nstep/10 steps al fine di termalizzarsi, dopodichè la procedura di
     for(int istep=1; istep <= (int)(nstep/10); ++istep){  // inizializzazione è fatta ripartire fino a che la temperatura del sistema 
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
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

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

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
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
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
