#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <ctime>
#include <random>

using namespace std;

#define J 1              /* Magnetic dipole spin */

#define SEED 201991748        /* Random generator */


// ================ Define global parameter ================== //


static double TMAX, T, dT;                   // Temperature loop //
static double NORM;                          // Normalization at the ensemble average //
static int N, NSPINS, MCS;                   // Half number of points, N*N, Monte Carlo steps


random_device rd;                          // Random Device //
mt19937 engine;                            // Engine of random device //
uniform_real_distribution<double> realDist;       // Produce a random real number between 0,1 //
uniform_int_distribution<int> intDist;           // poduce an integer number between 0,1 //


const string UP   = "\u21BF";   // spin - up
const string DOWN = "\u21C2";   // spin - down

// ============= Create lookup table of exp(-DeltaE / T) ===============
double exp_table[9];                                   // It can store 10 numbers //




/* ====== A function to Print out the lattice using arrows to indicate spin ========= */
/* ====== param the 2D array representing the lattice ======== */

void printSpins(int **lat)             // The input is a pointer to pointer 
{
  for (int x = 0; x < N; x++) 
    {
      for (int y = 0; y < N; y++) 
	{
	  string s;
	  switch (lat[x][y]) 
	    {
	    case 1:
	      s = UP;
	      break;
	    case -1:
	      s = DOWN;
	      break;
	    default:
	      s = "ERR";
            }
	  cout << s << " ";
        }
      cout << endl;
    }
}


/* ==================================================== */
/// The bond energy at point (x,y)
///
/// The bond energy between site i and neighbor j is defined as
///  -J * s_i * s_j
///  For all neighbors:
///  H_i = -J * s_i * sum(s_j)
///
///  @param lat is a 2D array of integers representing the lattice of spins
///  @param x is the x coordinate of the point in the lattice
//   @param y is the y coordinate of the point in the lattice
//   @result the bond energy of point (x,y) in the lattice
//   The important point the priodic boundary condition

int bondEnergy(int **lat, int x, int y) 
{
    /* Bond energy between site i and neighbor j is
     *          -J * s_i * s_j
     *
     * For all neighbors:
     *          H_i = -J * s_i * sum(s_j)
     *
     * Neighbor Key
     *
     *      o a o
     *      d i b
     *      o c o
     */
  //**************************************************//
  // This function accepts a 2D array and the position of each point.//
  // It should return the energy for any point in the lattice. lat[x][y] means //
  // to access a point in lattice for computing the energy for this point.//
  // The syntax (x+1)%N which means (x+1) mod N apply periodic boundary condition from the east of a point [x][y].//
  // The syntax (y+1)%N which means (y+1) mod N apply periodic boundary condition from the north of a point [x][y].//
  // The syntax x?(x-1):(N-1) which means if x is non zero consider x-1 and otherwise consider consider N-1 apply periodic boundary condition from the west of a point [x][y].//
  // The syntax y?(y-1):(N-1) which means if y is non zero consider y-1 and otherwise consider consider N-1 apply periodic boundary condition from the south of a point [x][y].//
  //**********************************************//.
  return (-J * lat[x][y] * ( lat[(x+1)%N][y] +  lat[x?(x-1):(N-1)][y] + lat[x][(y+1)%N] + lat[x][y?(y-1):(N-1)]   ));
}





/// Calculate and return the total energy of the lattice
///
/// The total bond energy for the lattice is sum(H_i)
/// The bond energy between site i and neighbor j is defined as
///  -J * s_i * s_j
///  For all neighbors:
///  H_i = -J * s_i * sum(s_j)
///
///  @param lat is a 2D array of integers representing the lattice of spins
///  @result the total energy of the lattice

int energy(int **lat) 
{
  int H = 0;
  /// Total bond energy for lattice is sum(H_i)
  // ***************************************//
  
  // We should move on every point of the lattice and calculate the energy point for each spin, then sum over all of them.//

  // ************************************//
  for(int x = 0; x < N; x++)
    {
      for(int y = 0; y < N; y++)
	{
	  H = H + 0.5*bondEnergy(lat, x, y);
	}
    }
  return (H);
  // ************************************//
  
  // We need a factor 0.5 to prevent two times calculate the energy.//
  
  // *********************************//
}
 


/// Calculate and return the total magnetization of a lattice
///
/// The total magnetization of a lattice is the sum of all lattice points
///
/// @param lat is a 2D array of integers representing the lattice of spins
/// @result is the net magnetization of the lattice

int magnetization(int **lat)
{
  // ************************************************************** //
  
  // We should move on every point of the lattice and calculate the magnetization over all of them.//
  // The initial value should be zero for sum//

  // ************************************************************* //
  int m = 0;
  for(int x = 0; x < N; x++)
    {
      for(int y = 0; y < N; y++)
	{
	  m = m + lat[x][y];
	}
    }
  return (m);
  // **************************************//
  
  // It is only a sum over all spins. //

  // **************************************//
  
	
}



/// Calculate whether or not to flip the spin of a lattice position
///
/// Return true if the spin at (x,y) of the lattice lat should be
/// flipped based on the Metropolis criterion at temperature temp
/// The energy change of the spin flip is
/// dE = 2 * -H_i = J * s_i * sum(s_j)
/// you can use the BondEnergy function to calculate the change in spin
/// will be flipped
/// Note: total change is twice the bond energy
///
/// @param lat is a 2D array of integers representing the lattice of spins
/// @param x is the x position in the lattice of the spin to be flipped
/// @param y is the y position in the lattice of the spin to be flipped
/// @param dE is the change in energy corresponding to the spin flip
///        it is passed by reference and should be updated for the main
///        function
/// @param temp is the temperature
/// @result true if the spin at (x,y) should be flipped and otherwise false

bool testFlip(int **lat, int x, int y, int &dE, double temperature) 
{
  // CMSC6920: calculate the change in energy corresponding to the given
  // spin flip and apply the Metropolis criterion to it to determine if
  // the spin should be flipped
  /* Energy change upon reversal of spin is
   *          dE = 2 * -H_i = J * s_i * sum(s_j)
   *
   * Note: total change is twice the bond energy
   */
  dE = -2 * bondEnergy(lat, x, y);
  double zeta= 1.0 * rand()/RAND_MAX;               // store the random number in this variable
  double Pa= 0.0;                                   // store acceptance probablity here
  bool acceptance=false;

  // *************************************** //
  // Check the condition : //
  if(dE <= 0 )
    {
      Pa = 1; // means accept
    }
  if(dE > 0 && zeta <= exp(-dE/T))
    {
      Pa = 1; // means accept
    }
  if(Pa == 1)
    {
      acceptance  = true; // define accept
    }
    {
  // set the breakpoint on this line
      return(acceptance);
    }
}



// energies tested can only hold values of 0, 2, 4, 8
// for efficiency, the lookup table should be configured so that exp_table[dE] directly returns exp(-dE/T)
// i.e., exp_table[i=2] should access exp( -(dE=2) / T)
// note that negative values of dE are always accepted, so they should not be tested by exp_table

void initExpLookup(double T)
{
  for(int dE = 0; dE < 9; dE = dE + 1)
    {
      if(dE > 0) // Saman : Check the condition // 
	{
	  exp_table[dE] = exp(-dE/T);
	}
      
    }
  
}


/* Now, we can run a simulation */


int main(int argc, char *argv[]) 
{
  N=4;                                                // This program is only for 4*4 lattice
  
  engine=mt19937(SEED);                              // random engine
  uniform_int_distribution<> intDistLattice(0, 1);                   // pick up integer from zero and one
  uniform_int_distribution<> intLattice(0, N);                       // pick up integer between zero and N
  
  int **lattice;	                                             // A pointer to pointer for lattice



   // allocate and initialize NxN lattice and assign each point a spin of J

  // *****************************************************//

  // We need to initialize the pointer in order to access to the memory //
  lattice = new int* [N];
  
  for (int k=0; k<N; k++)  // Start to move on lattice //
	{
	  lattice[k]=new int[N]; //  Initialize the pointer //
	}
  
  //  Make the lattice that is requested //
  for(int x = 0; x < N; x++)
    {
      for(int y = 0; y < N; y++)
	{
	  lattice[x][y]  = J;
	}
    }
    
  // Observables 
  //   M - magnetization
  //   E - energy

  
  
  
  int M;
  double E; 



  // test case 1
  cout << "Tests of energy and magnetization functions" << endl << endl;
  cout << "Test Case 1" << endl;
  cout << "-----------" << endl;
  printSpins(lattice);
  E = energy(lattice);
  M = magnetization(lattice);       
  
  cout << "Energy: " << E << " (should be -32)" << endl;
  cout << "Magnetization: " << M << " (should be 16)" << endl;
  cout << endl;
  
  // test case 2
  cout << "Test Case 2" << endl;
  cout << "-----------" << endl;
  
  for (int x = 0; x < N; x++)
    {
      for (int y = 0; y < N; y++)
	{
	  lattice[x][y]=-J;
	}
    }
    
  printSpins(lattice);
  E = energy(lattice);
  M = magnetization(lattice);       
  
  cout << "Energy: " << E << " (should be -32)" << endl;
  cout << "Magnetization: " << M << " (should be -16)" << endl;
  cout << endl;
  
  // test case 3
  cout << "Test Case 3" << endl;
  cout << "-----------" << endl;
  
  for (int x = 0; x < N; x++)
    {
      for (int y = 0; y < N; y++)
	{
	  if( (x+y) % 2 ==0)
	    {
	      lattice[x][y]=J;
	    }
	  else
	    {
	      lattice[x][y]=-J;
	    }
	}
    }
  
  printSpins(lattice);
  E = energy(lattice);
  M = magnetization(lattice);       
  
  cout << "Energy: " << E << " (should be 32)" << endl;
  cout << "Magnetization: " << M << " (should be 0)" << endl;
  cout << endl << endl;
  
  cout << "Tests of bondEnergy" << endl << endl;
  
  for (int x = 0; x < N; x++)
    {
      lattice[x] = new int[N];
	
      for (int y = 0; y < N; y++)
	{
	  lattice[x][y]=J;
	}
    }
  printSpins(lattice);
  cout << "Bond energy at (0,0) = " << bondEnergy(lattice, 0, 0) << " (should be -4)" << endl;
  cout << "Bond energy at (1,1) = " << bondEnergy(lattice, 1, 1) << " (should be -4)" << endl;
  cout << "Bond energy at (3,3) = " << bondEnergy(lattice, 3, 3) << " (should be -4)" << endl;
    
  cout << "-----------" << endl;
  
  for (int x = 0; x < N; x++)
    {
      for (int y = 0; y < N; y++)
	{
	  if( (x+y) % 2 ==0)
	    {
		lattice[x][y]=J;
	    }
	  else
	    {
	      lattice[x][y]=-J;
	    }
	}
    }
  
  printSpins(lattice);
  
  cout << "Bond energy at (0,0) = " << bondEnergy(lattice, 0, 0) << " (should be 4)" << endl;
  cout << "Bond energy at (1,1) = " << bondEnergy(lattice, 1, 1) << " (should be 4)" << endl;
  cout << "Bond energy at (3,3) = " << bondEnergy(lattice, 3, 3) << " (should be 4)" << endl;
  cout << "-----------" << endl;
  
  for (int x = 0; x < N; x++)
    {
      for (int y = 0; y < N; y++)
	{
	  if( (x+y) % 2 ==0)
	    {
	      lattice[x][y]=J;
	    }
	  else
	    {
	      lattice[x][y]=-J;
	    }
	}
    }
  printSpins(lattice);
  
  cout << "Bond energy at (0,0) = " << bondEnergy(lattice, 0, 0) << " (should be 4)" << endl;
  cout << "Bond energy at (1,1) = " << bondEnergy(lattice, 1, 1) << " (should be 4)" << endl;
  cout << "Bond energy at (3,3) = " << bondEnergy(lattice, 3, 3) << " (should be 4)" << endl;

  cout << "---------" << endl;
  for (int x = 0; x < N; x++)
    {
      for (int y = 0; y < 3; y++)
	  {
	    lattice[x][y]=J;
	  }
    }
  
  for (int x = 0; x < N; x++)
    {
      for (int y = 3; y < N; y++)
	{
	  lattice[x][y]=-J;
	}
    }
  
  printSpins(lattice);
  int dE=4;
  double temperature=100.0;
  
  cout << "Bond energy at (0,0) = " << bondEnergy(lattice, 0, 0) << " (should be -2)" << endl;
  cout << "Bond energy at (1,1) = " << bondEnergy(lattice, 1, 1) << " (should be -4)" << endl;
  cout << "Bond energy at (2,2) = " << bondEnergy(lattice, 2, 2) << " (should be -2)" << endl;
  cout << "Bond energy at (3,3) = " << bondEnergy(lattice, 3, 3) << " (should be 0)" << endl;
  
  initExpLookup(5);
  cout << "exp lookup table" << endl;
  for(int i = 0 ;i<9;++i)
    {
      cout << "i = " << i << ", E = " << i <<  ", exp(-E/T) = " << exp_table[i] << endl;
    }
  testFlip(lattice, intLattice(engine), intLattice(engine), dE, temperature);


  // We want to deallocate memory //

  for(int l = 0; l < N; l++)                 // We should free the memory for new size of lattice
    {
      delete[] lattice[l];
    }
  delete[] lattice;
      
  return(0);
}
