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

#define SEED srand(time(0))        /* Random generator */


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
