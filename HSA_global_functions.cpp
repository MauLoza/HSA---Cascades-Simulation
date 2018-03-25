#include "HSA.h"
#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include <vector>
#include <fstream>
#include <stdio.h>

using namespace std;

//Define global variables--------------------------------------|
double mass[4] = {134.9766*10e6, 139.57018*10e6, 139.57018*10e6, 105.6583*10e6};   
// neutral, positive, negative, muon (MeV/c^2)
double mean_time = 2.55e-8;                                                               
// lifetime of pion in lab frame (s)
double c = 3.0e5;  // km/s

//Define global functions--------------------------------------|
double km_to_gcm(double h){ return 1030*exp(-h/6.4); }
double gcm_to_km(double g){ return 6.4*log(1030/g);}
double lorentz_factor(double E_0){ return (E_0/mass[1])+1; } 
//Mass already includes /c^2

double exp_inv(double lambda){ 
  double U=rand_unif_real(0,1);          //Random number  
  return (-1/lambda)*log(1-U); 
}

double exp_lambda(double lambda){
  
  double l;
  random_device rand_dev;
  mt19937 generator(rand_dev());
  exponential_distribution<double> distr(lambda);
  
  l = distr(generator);
  return l;
}

double rand_unif_real(double E_0, double E_f){

  std::random_device              rand_dev;
  std::mt19937                        generator(rand_dev());
  
//Declare uniform distribution generator
  std::uniform_real_distribution<double> distr(E_0, E_f);
 
  double energy = distr(generator);
  return energy;
}

vector<double> Split_En_Unif_2N(double E_0, int N){

  vector<double> energy;  
  double i_energy;
  double E_initial;
  double E_sum;

  E_initial = E_0;

//Split energy
  for(int i=0; i<(pow(2,N)-1); i++){
  
  i_energy = rand_unif_real(0.0, E_0);
  energy.push_back(i_energy);
  
  E_sum += i_energy;  
  E_0 -= i_energy;
    
  }

  energy.push_back(E_initial-E_sum);
  return energy;
}

int pion_chooser(){

  int mu;
  std::random_device              rand_dev;
  std::mt19937                        generator(rand_dev());
  
//Declare uniform distribution generator
  std::uniform_int_distribution<int> distr(0, 2);
  
  mu = distr(generator);
  return mu; 

//0 equals neutral pion, 1 equals positive pion
//2 equals negative pion, 3 muon and 4 proton  
}

