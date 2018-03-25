#ifndef HSA_ALGORITHM_H
#define HSA_ALGORITHM_H

//Includes all functions that work with struct Particle

#include <iostream>
#include <chrono>
#include <random>
#include <math.h>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;

struct Particle{
  double energy;  //Initial energy
  int type; //Kind of particle: 0 neutral, 1 positive, 2 negative, 4 proton
  double distance; //g/cm^2
  string status;       //Particle interacts or decays
  
  Particle(double E, int kind, double d){
  energy = E;
  type = kind;  
  distance = d;
  status = "i";
  }  
};

//Complete split of energy into secondaries (pions). HSA ALGORITHM

void Writte(Particle p, ofstream& file){  file << p.energy << " " << p.type << 
" " << p.distance << " " << p.status << endl; }

void HSA_secondaries(Particle& p, vector<double> branch_energy, 
double E_th, int N, ofstream& energy_type, string name){

  vector<double> temporal_energy;
  string file = name.append(".dat");
  energy_type.open(file, ios_base::app);  //Append in file 
  
  temporal_energy = Split_En_Unif_2N(p.energy, 1); 
  //1. split energy of incoming particle
  p.energy = temporal_energy[0];  
  //2. rewrite energy of incoming particle
  Writte (p, energy_type);
  //3. split energy in 2^N branches  
  branch_energy = Split_En_Unif_2N(temporal_energy[1], N);       
  
//GENERATE SECONDARIES - pions
   
  for(int i=0; i<pow(2,N); i++){
    //4. divide energy of branch [kin_energy, remaining energy]
    temporal_energy = Split_En_Unif_2N(branch_energy[i], 1);       
    
    do{
      
      Particle s(temporal_energy[0], pion_chooser(), p.distance);       
      //5. assign kin_energy for pion, temporal[0]
      Writte(s, energy_type);  //Save info of secondary     
      temporal_energy[1] -= mass[s.type];                                            
      //6.1 Substract pion mass from remaining energy
      temporal_energy = Split_En_Unif_2N(temporal_energy[1], 1);  
      //6.2 Split the remaining energy again
    
    }while (temporal_energy[1]>E_th);   
  }
  
  energy_type.close();
}

//Read secondaries and select only hadronic components 
//Vector selector_hadronic is a vector of particles                                                                                

void selector_hadronic(vector<Particle>& p_vec, 
ifstream& particles_now, string name){

  double temporal[3];
  char c_temp;
  string file =  name.append(".dat");
  particles_now.open(file);     
  //Opens file. Write as energy, type, distance, status

  while(particles_now >>temporal[0] >>temporal[1] >>temporal[2] >> c_temp){       

    if (temporal[1]==0) continue; 
    //Ignore electromagnetic components
    Particle p(temporal[0], temporal[1], temporal[2]);
    p.status = c_temp;
    p_vec.push_back(p);
  }
  particles_now.close();
} 

//Lambda by particle

double lambda_int_proton(double E_0){
  if(E_0 <= 0.1e12) return 87.0;
  else return (80.8 - 2.78*log(E_0/1e12));
}

double lambda_int_pion(double E_0){
  if(E_0 <= 0.1e12) return 116;
  else return (105 - 4.23*log(E_0/1e12));
}

double beta(Particle p){
double factor;
factor = p.energy/mass[p.type]+1;
return  sqrt(1-pow(1/factor,2));
}

double lambda_dy_pion(Particle p){
  return (mean_time*c*lorentz_factor(p.energy));
}

//Energy split for decays of charged pions and selector of decay lenght

vector<double> Split_Decay_Pion(Particle p){

  vector<double> energy;  
  double i_energy;
  double E_in, E_fin;
  double E_sum;
  double E_0;           //initial value for flat distribution

  E_in = p.energy;
  E_0 = (mass[3]*mass[3])/(mass[1]*mass[1])*p.energy;

//Split energy
  
  i_energy = rand_unif_real(E_0, p.energy);
  energy.push_back(i_energy);
  
  E_sum += i_energy;  
  E_0 -= i_energy;

  energy.push_back(E_in-E_sum);
  return energy;
}
double Distance_Dy_Pion(Particle p){

  double h_1, h_2;  //Distance in km: h_1 initial and h_2 final
  double X_2;          //Final distance in g/cm^2
  double l_d;            //lenght of decay km

  h_1 = gcm_to_km(p.distance);  //Original height in g/cm^2
  l_d = exp_lambda(lambda_dy_pion(p)); 
  h_2 = h_1-l_d;   //Decay lenght in km
  X_2 = km_to_gcm(h_2);   //Decay length in g/cm^2
  return X_2 - p.distance; //Final decay lenght in g/cm^2
}

//Final function that determines the new particle's distance
void distance(Particle& p){
  
  double d_dy, d_int;
  double l;
  
  if (p.type == 4){
     if (p.energy>10e9) p.distance += exp_lambda(lambda_int_proton(p.energy));  
    else p.distance += 1030;  
    //Proton reaches sea level at energies below 10GeV
  }
  
  else if ((p.type==1)||(p.type==2)){

    if(p.energy>1e12)  p.distance += exp_lambda(lambda_int_pion(p.energy));
    else if(p.energy<10e9) {
      p.status = "d"; 
      p.distance += Distance_Dy_Pion(p);
    }
    else if(10e9<p.energy<1e12){ 
         d_int = exp_lambda(lambda_int_pion(p.energy));
         l = km_to_gcm(lambda_dy_pion(p));
         d_dy = exp_lambda(l);  

      //Return shortest distance for charged pions
      if(d_int>d_dy){ 
        p.status = "d";
        p.distance += Distance_Dy_Pion(p);
      }
      else if (d_int<d_dy) p.distance += d_int;
    }
  }
}

void move_particles(vector<Particle>& p_vec){
  for (int i = 0; i<p_vec.size(); i++) { 
    
    distance(p_vec[i]);
    if(p_vec[i].status=="d"){ 
    p_vec[i].type=5;  //If pion decays turn it into muon
    p_vec[i].energy= Split_Decay_Pion(p_vec[i])[0]-mass[4];  
    //Rewrite its energy according to function, assign energy of muon.
    p_vec[i].distance = 1030;    //Displace it til sea level, as muon
    }
  }
}

//Write vector of particles
void Write_p_vec(vector<Particle>& p_vec, 
ofstream& particles_after, string name){
  string file = name.append(".dat");
  particles_after.open(file, ios_base::app);
  for (int i = 0; i<p_vec.size(); i++)  Writte(p_vec[i], particles_after);
  particles_after.close();
}

//Remove for energy below threshold & particles with bigger density than sea
//except muons
void Erase_below_Eth(vector<Particle>& p_vec, double E_th){
    p_vec.erase( 
      remove_if(p_vec.begin(), p_vec.end(), [&](Particle const & p){
//      if(p.type==5) return false;
      return p.energy<E_th; }), 
      p_vec.end() ); 
}

void Erase_below_sea(vector<Particle>& p_vec){
    p_vec.erase( 
      remove_if(p_vec.begin(), p_vec.end(), [&](Particle const & p){
      if(p.type==5) return false;
      return p.distance>1030; }), 
      p_vec.end() ); 
}

void Erase_muons(vector<Particle>& p_vec){
    p_vec.erase( 
      remove_if(p_vec.begin(), p_vec.end(), [&](Particle const & p){
        return p.type==5; }), 
      p_vec.end() ); 
}

void Erase_pions(vector<Particle>& p_vec){
    p_vec.erase( 
      remove_if(p_vec.begin(), p_vec.end(), [&](Particle const & p){
        return p.type!=5; }), 
      p_vec.end() ); 
}

void Erase_muons_Th(vector<Particle>& p_vec, double E_th_mu){
    p_vec.erase( 
      remove_if(p_vec.begin(), p_vec.end(), [&](Particle const & p){
        return p.energy<E_th_mu; }), 
      p_vec.end() ); 
}

//Muon saver

void save_delete_muons(vector<Particle>& p_vec, ofstream& particles_after, 
int& counter, double E_th_mu){
 
  string name = "muons";
  vector<Particle> p_muons = p_vec;     //Copy p_vec
  Erase_pions(p_muons);
  Erase_muons(p_vec);
  Erase_muons_Th(p_muons);
  counter += p_muons.size();
  Write_p_vec(p_muons, particles_after, name);

}

#endif
