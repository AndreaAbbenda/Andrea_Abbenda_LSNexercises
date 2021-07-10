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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  true_false = true;
  T = -1;
  Input(T);
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation                            //La struttura del main è stata mantenuta pressoché
  {                                                                             //inalterata: l'unica differenza ho aggiunto 
    Reset(iblk);   //Reset block averages                                       //una parte di codice per poter campionare
    for(int istep=1; istep <= nstep; ++istep)                                   //il sistema a temperatura variabile. A tal proposito
    {                                                                           //ho modificato la funzione Input affinchè fosse 
      Move(metro);                                                              //possibile variare la temperatura del sistema
      Measure();
      Accumulate(); //Update block averages               //la variabile booleana true_false è usata come semplice 'interruttore':
    }                                                     //nello studiare il sistema a t variabile acquisisco certi dati, nel caso
    Averages(iblk);   //Print results for current block   //invece di t costante ne acquisisco altri
  }
  
  ConfFinal(); //Write final configuration
  
  true_false = false;
  for (int t = 5; t <=20; t++)
  {  
    T = (double)t/10;
    cout << "System temperature: " << T << endl;

    Input(T);
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages                                 //dall'analisi dei risultati è possibile notare come il sistema
      for(int istep=1; istep <= nstep; ++istep)                             //equilibri dopo circa 40 blocchi. per questo motivo ho deciso
      {                                                                     // di acquisire misure a T variabile solo dopo il 40° blocco
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      if (iblk >= 40){Averages(iblk);}  //Print results for current block
    }
  }


  return 0;
}


void Input(double T)
{
  ifstream ReadInput, ReadConf;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  if (T < 0) {ReadInput >> temp;}
  else {ReadInput >> temp; temp = T;}
  beta = 1.0/temp;
  ReadInput >> nspin;

  ReadInput >> J;

  ReadInput >> h;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
  
  if(true_false == true){
    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    cout << "Temperature = " << temp << endl;
    cout << "Number of spins = " << nspin << endl;
    cout << "Exchange interaction = " << J << endl;
    cout << "External field = " << h << endl << endl;
    if(metro==1) cout << "The program perform Metropolis moves" << endl;
    else cout << "The program perform Gibbs moves" << endl;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();
  }

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  ReadConf.open("config.final");
  string yes_no;
  if (ReadConf.is_open()){
    if(true_false == true){
      cout << "Continue previous simulation? [y/n]"  << endl;     
      cin >> yes_no;
    }else if(true_false == false){yes_no = "n";}
    if (yes_no == "y"){            //Inizialization     // ho modificato la funzione Input(T) affinchè il sistema possa ripartire 
         for (int i=0; i<nspin; ++i){                   // con l'ultima configuraziona della simulazione precedente. 
           ReadConf >> s[i];
         }
                                                     
    }else if (yes_no == "n"){                     
         for (int i=0; i<nspin; ++i){
            if(rnd.Rannyu() >= 0.5) s[i] = 1;
            else s[i] = -1;
         }
    };
  }else{
     for (int i=0; i<nspin; ++i){
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
     }
  }
  ReadConf.close();  
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  if(true_false == true){cout << "Initial energy = " << walker[iu]/(double)nspin << endl;}
}

//--------------------------------------------------------------------------------------------------

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
//  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted += 1;
    if(metro==1) //Metropolis
    {
      double test = rnd.Rannyu();
      sm = s[o] * -1;      
      energy_new = exp(-beta*Boltzmann(sm, o));
      energy_old = exp(-beta*Boltzmann(s[o], o));
      p = Min(1, energy_new/energy_old);
      if (p > test){s[o] = sm; accepted += 1;}
    }
    else //Gibbs sampling
    {
      double test = rnd.Rannyu();      
      sm = s[o]*(-1);
      energy_new = Boltzmann(sm, o);
      energy_old = Boltzmann(s[o],o);
      p = exp(-beta*energy_new)/(exp(-beta*energy_new) + exp(-beta*energy_old));
      if (p > test){s[o] = sm; accepted += 1;}
    }
  
  }
}

//--------------------------------------------------------------------------------------------------

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

double Min(const double a, const double b)
{
  double ret;
  if(a > b){ret = b;}
  else if (b >= a){ret = a;}
  return ret;
}

//--------------------------------------------------------------------------------------------------

void Measure(void)
{
//  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }  
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}

//--------------------------------------------------------------------------------------------------

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
           glob_av_temp[i] = 0;
           glob_av2_temp[i] = 0;

       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

//--------------------------------------------------------------------------------------------------

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

//--------------------------------------------------------------------------------------------------

void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=25;
    
    if(true_false == true){
      cout << "Block number " << iblk << endl;
      cout << "Acceptance rate " << accepted/attempted << endl << endl;
    }
    
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    if (true_false == true){                         //misure a temperatura fissa
      Ene.open("output.ene.0.metro",ios::app);
      Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
      Ene.close();     
    }else if(true_false == false){                   //misure a temperatura varibile. Acquisisco dal 40-esimo blocco in poi
     
        glob_av_temp[iu]  += stima_u;
        glob_av2_temp[iu] += stima_u*stima_u;
	if (iblk == nblk){
	     err_u=Error(glob_av_temp[iu],glob_av2_temp[iu],iblk-39);
             Ene.open("output.ene.temp.metro",ios::app);
             Ene << setw(wd) << T << setw(wd) << stima_u << setw(wd) << glob_av_temp[iu]/(double)(iblk-39) << setw(wd) << err_u << endl;
             Ene.close();
             }       
    }

    stima_c = beta*beta*(blk_av[ic]/blk_norm/(double)nspin - stima_u*stima_u*(double)nspin); //Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    if (true_false == true){ 
      Heat.open("output.heat.0.metro",ios::app);
      Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
      Heat.close();     
    }else if(true_false == false){     

        glob_av_temp[ic]  += stima_c;
        glob_av2_temp[ic] += stima_c*stima_c;
	if (iblk == nblk){
	     err_c=Error(glob_av_temp[ic],glob_av2_temp[ic],iblk-39);
             Heat.open("output.heat.temp.metro",ios::app);
             Heat << setw(wd) << T << setw(wd) << stima_c << setw(wd) << glob_av_temp[ic]/(double)(iblk-39) << setw(wd) << err_c << endl;
             Heat.close();
             }       
    }
    
    stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    if (true_false == true){ 
      Mag.open("output.mag.002.metro",ios::app);
      Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
      Mag.close();     
    }else if(true_false == false){
      
        glob_av_temp[im]  += stima_m;
        glob_av2_temp[im] += stima_m*stima_m;
	if (iblk == nblk){
	     err_m=Error(glob_av_temp[im],glob_av2_temp[im],iblk-39);
             Mag.open("output.mag.002temp.metro",ios::app);
             Mag << setw(wd) << T << setw(wd) << stima_m << setw(wd) << glob_av_temp[im]/(double)(iblk-39) << setw(wd) << err_m << endl;
             Mag.close();
        }
    }

    stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    if (true_false == true){ 
      Chi.open("output.chi.0.metro",ios::app);
      Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
      Chi.close();
      cout << "----------------------------" << endl << endl;
    }else if(true_false == false){
     
        glob_av_temp[ix]  += stima_x;
        glob_av2_temp[ix] += stima_x*stima_x;
	if (iblk == nblk){
	     err_x=Error(glob_av_temp[ix],glob_av2_temp[ix],iblk-39);
             Chi.open("output.chi.temp.metro",ios::app);
             Chi << setw(wd) << T << setw(wd) << stima_x << setw(wd) << glob_av_temp[ix]/(double)(iblk-39) << setw(wd) << err_x << endl;
             Chi.close();
             }       
    }
}

//--------------------------------------------------------------------------------------------------

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

//--------------------------------------------------------------------------------------------------

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

//--------------------------------------------------------------------------------------------------

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
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
