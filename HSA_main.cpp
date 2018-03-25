#include "HSA.h"
#include "HSA_algorithm.h"

using namespace std;

int main()
{
  int points = 4;
  double E_0;       //Unity is eV
  int N[4]= {0, 1, 2, 3};  //2^N branches
  vector<double> branch_energy;
  double E_th_HSA = 250e6;  //Threshold
  double E_th_mu = 1e3;         //Muon threshold
  ifstream particles_now;         //File to read .txt
  ofstream particles_after;       //File to write .txt
  vector<Particle> p_vec;        //Vector that will include all secondaries
  double height=20.0;              //In km
  string name="secondaries";
  char number[4][12] = {"zero30.dat", "one30.dat", "two30.dat", "three30.dat"}; 
  //name of out in function of the branch J.dat

//Remove any possible repetition of our data
  remove(number[1]);
  remove(number[2]);
  remove(number[3]);
  remove(number[4]);
  remove("muons.dat");

for(int k=0; k<4; k++){
particles_after.open(number[k]); //write info for file number[k]
for(int i=2; i<4; i++){

  E_0=pow(10,i)*1e12;

for (int j=0; j<1; j++){
  //FIRST COLLISION
  int counter=0;
  Particle p(E_0,  4,  km_to_gcm(height));   //4 stands for proton at height in km
  p_vec.push_back(p);                                 //p_vec saves proton p

  //SUBSEQUENT COLLISIONS
  //Develop secondaries-All particles left have status i, 
  //and muons are excluded with if(p_vec[i].type !=5)
  //Only muons have status d, because pions decaying were turned into muons.


do
  {
  remove("secondaries.dat");
  move_particles(p_vec); 
  //Moved secondaries and decay some of them
  Erase_below_Eth(p_vec, E_th_HSA); 
  //Delete secondaries below minimum energy, except muons
  Erase_below_sea(p_vec); 
  //Delete secondaries below sea level, except muons
  save_delete_muons(p_vec, particles_after, counter, E_th_mu); 
  //Delete muons from p_vec and write file with their info

  for(int m= 0; m<p_vec.size(); m++)  if(p_vec[m].type !=5){ 
  HSA_secondaries(p_vec[m], branch_energy, E_th_HSA, 
  N[k], particles_after, name);
  }
  //All remaining particles generate secondaries except muons.
  //HSA produces file with secondaries. 
  //Inside the loop keeps writing all the remaining particles 
  p_vec.clear(); 
  //Remove info from original branches, they are already included in name.dat
  selector_hadronic(p_vec, particles_now, name); 
  //add elements from HSA to vector with info of secondaries (no neutral pions)
  Erase_below_Eth(p_vec, E_th_HSA);

  }while(p_vec.size()>1);

  p_vec.clear();

  particles_after.open(number[k], ios_base::app); 
  //write info of muons reaching sea level
  particles_after  << counter << " ";
  }

  particles_after << endl;
  particles_after.close();
  }
  }
}




