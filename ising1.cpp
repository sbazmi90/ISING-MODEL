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

bool testFlip(int **lat, int x, int y, int &dE, double temp) 
{
  // CMSC6920: calculate the change in energy corresponding to the given
  // spin flip and apply the Metropolis criterion to it to determine if
  // the spin should be flipped
  dE = -2 * bondEnergy(lat, x, y);
  
  if (dE <= 0)
    {
      return true;
    }
  else
    {
      double xi = realDist(engine);
      if (xi < exp_table[dE])
	{
	  return true;
	}
    }
  return(false);
}

void initExpLookup(double T)
{
  for(int i=0;i<9;++i)
    {
      exp_table[i]=exp(-i/T);
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









// Run a simulation

int main(int argc, char *argv[]) 
{
  // Initialize global parameter variables from input file Exit with
  // error code 1 if input file has errors
  
  if (!argv[1]) 
    {
      cerr << "Specify an input file" << endl;
      return(1);
    }
  
  if (read_input(argv[1])) 
    {
      return(1);
    }
  
  if (!argv[2]) 
    {
      cerr << "Specify an output file" << endl;
      return 1;
    }
  
  ofstream ofile(argv[2]);
  
  if (!ofile) 
    {
      cerr << "Could not open output file! " << argv[2] << endl;
      return(1);
    }
  
  ofile << "#N " << N << endl;
  ofile << "#EQMCS " << EQMCS << endl;
  ofile << "#PMCS " << PMCS << endl;
  ofile << "#T " << T << endl;
  ofile << "#dT " << dT << endl;
  ofile << "#TMAX " << TMAX << endl;
  
  // Keep record of seed 
  ofile << "#SEED " << SEED << endl;
  ofile << "#" << endl;
  
  engine=mt19937(SEED);
  
  uniform_int_distribution<> intLattice(0, N-1);
  uniform_real_distribution<> realDist(0, 1);
  uniform_int_distribution<> intDistLattice(0, 1);
  
  NSPINS = N*N;                   // Number of spins
  NORM = 1.0 / (PMCS * NSPINS);    // Normalization constant
  
  int **lattice;
  
  // We need to initialize the pointer in order to access to the memory //
  lattice = new int* [N];
  
  for (int k=0; k<N; k++)  // Saman : Start to move on lattice //
	{
	  lattice[k]=new int[N]; // Saman : Initialize the pointer //
	}
  

  for (int x = 0; x < N; x++)           // Move on x direction 
    {
      for (int y = 0; y < N; y++)      // Move on y direction
	{
	  if( rd()/RAND_MAX < 0.5 )     //Determine the spin in lattice based on random generation
	    {
	      lattice[x][y] = J;
	    }
	  else
	    lattice[x][y] = -J;
	}
    }
  /************************ Saman ***********************/
  /* Observables 
   * M - magnetization
   * E - energy */
  
  int M = 0;
  double E = 0;        // for a given state
  double M_avg, E_avg;    // ensemble averages
  double M_tot, E_tot;    // sum for all states at a given temperature 

  // Column headers for data file 
  char s[68];
  sprintf(s, "%-7s%12s%12s%12s\n", "#TEMP", "<E>", "<E>^2", "|<M>|");
  ofile << s;
  
  // loop over temperatures from T to TMAX
  // and compute E_ave and M_ave at each temperature
  // by performing a Monte Carlo simulation at each 
  // temperature
