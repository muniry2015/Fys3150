#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;

ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {//<<<{refresh meaning of inline}>>>
    //modulo operator, represented by a percentage sign (%), gives the remainder of a division of two values
    //  x = 11 % 3;
   // results in variable x containing the value 2, since dividing 11 by 3 results in 3, with a remainder of 2.

  return (i+limit+add) % (limit);
}
// Function to read in data from screen
void read_input(int&, int&, double&, double&, double&);
// Function to initialise energy and magnetization
void initialize(int, double, int **, double&, double&, bool random=false);
// The metropolis algorithm
void Metropolis(int, long&, int **, double&, double&, double *, int&);
// prints to file the results of the calculations
void output(int, int, double, double *, int, double);

int main(int argc, char* argv[])
{
  //At start state all spins =+1, in case we want random start uncomment random.
  char *outfilename;
  long idum;
  int **spin_matrix, n_spins, mcs, flips;
  double w[17], average[5], initial_temp, final_temp, E, M, temp_step;

  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] <<
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename);
  //    Read in initial values such as size of lattice, temp and cycles
  read_input(n_spins, mcs, initial_temp, final_temp, temp_step);
  spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));


  idum = -1; // random starting point<<{this is used to initialize the method ran1 found in lib.cpp}
  for ( double temperature = initial_temp; temperature <= final_temp; temperature+=temp_step){

    flips=0;
    // setup array for possible energy changes
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temperature);
    // initialise array for expectation values
    for( int i = 0; i < 5;i++) average[i] = 0.;//{it is used to stor the following values E,E^2, M,M^2, |M|,flips}
    //the initialization calculate E and M we start with giving a start temperature and the size of the grid.
    //we send the info n_spins, temp, the spin:matrix and the adresse of E and M to updated the initial
    //values of the E & M of the  system.
    if(temperature == initial_temp)
    {
        //    initialise energy and magnetization
        E = M = 0.;
        initialize(n_spins, temperature, spin_matrix, E, M);
    }
    // start Monte Carlo computation
      //int cycles ;
    for (int cycles = 1; cycles <= mcs; cycles++)
    {
      Metropolis(n_spins, idum, spin_matrix, E, M, w, flips);
      // update expectation values
      //NB:for each Metro cycle we get a possible state, note the temp is unchanged so those are just possible states
      //for the one giving temp. The simulation should bring the system to the most possible
      //  E and M giving spesified temperature as it stablise.


      average[0] += E;    average[1] += E*E;
      average[2] += M;    average[3] += M*M; average[4] += fabs(M);

     //output(n_spins, cycles, temperature ,average, flips,E); //this line is used in case we want to store the averages as function of mcs
    }
    // print results

    cout << "Your flips prosent is "<< (double)flips/n_spins/n_spins/mcs <<"\n";
    output(n_spins, mcs, temperature,average, flips,E);//This line in case u want to save the data as function of tempirature
  }
  free_matrix((void **) spin_matrix); // free memory
  ofile.close();  // close output file
  return 0;
}


// read in input data
void read_input(int& n_spins, int& mcs, double& initial_temp,
        double& final_temp, double& temp_step)
{
  cout << "Number of Monte Carlo trials =";
  cin >> mcs;
  cout << "Lattice size or number of spins (x and y equal) =";
  cin >> n_spins;
  cout << "Initial temperature with dimension energy=";
  cin >> initial_temp;
  cout << "Final temperature with dimension energy=";
  cin >> final_temp;
  cout << "Temperature step with dimension energy=";
  cin >> temp_step;
} // end of function read_input




// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, double temperature, int **spin_matrix,
        double& E, double& M, bool random)
{



    if (random==true)
    {
        //a negative number indecate the initialitation of the random function.
        long idme=-1;
        for(int y =0; y < n_spins; y++)
        {
          for (int x= 0; x < n_spins; x++)
          {

              if ( ran1(&idme) <= 0.5 )
              {
            spin_matrix[y][x] = -1;
            //cout<< spin_matrix[y][x] ;// this was used as part of seeing the matrix to be sure the the method is working
              }
              else {
                 spin_matrix[y][x] = 1; // spin orientation for the ground state
                // cout<< spin_matrix[y][x] ;// this was used as part of seeing the matrix to be sure the the method is working
              }
          if (idme==-1)idme=1;
          }
           //cout<< "\n" ;// this was used as part of seeing the matrix to be sure the the method is working
        }


       cout<<"Our friend Random says Hello. ;)  \n";//random setup is not the normal run so I want a notice when it's on.
    }


  // setup spin matrix and intial magnetization
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      spin_matrix[y][x] = 1; // spin orientation for the ground state
      M +=  (double) spin_matrix[y][x];
    }
  }
  // setup initial energy
  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      E -=  (double) spin_matrix[y][x]*
    (spin_matrix[periodic(y,n_spins,-1)][x] +
     spin_matrix[y][periodic(x,n_spins,-1)]);//This ()return avalue only if we are in the border of the
      //lattice else it returns 0
    }
  }
}// end function initialise

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, int &flips)
{
  // loop over all spins

  for(int y =0; y < n_spins; y++) {
    for (int x= 0; x < n_spins; x++){
      int ix = (int) (ran1(&idum)*(double)n_spins);
      int iy = (int) (ran1(&idum)*(double)n_spins);
      int deltaE =  2*spin_matrix[iy][ix]*
    (spin_matrix[iy][periodic(ix,n_spins,-1)]+
     spin_matrix[periodic(iy,n_spins,-1)][ix] +
     spin_matrix[iy][periodic(ix,n_spins,1)] +
     spin_matrix[periodic(iy,n_spins,1)][ix]);
      if ( ran1(&idum) <= w[deltaE+8] ) {
    spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config

    //we need to update the M to have 2 multi the new spin
    //look at it as the following, one for erasing the previous value and one add the new spin value.
        M += (double) 2*spin_matrix[iy][ix];
        E += (double) deltaE;
        flips++;
      }
    }
  }
} // end of Metropolis sampling over spins




void output(int n_spins, int mcs, double temperature, double *average, int flips, double E )
{
  double norm = 1/((double) (mcs));  // divided by total number of cycles
  double Eaverage = average[0]*norm;
  double E2average = average[1]*norm;
  double Maverage = average[2]*norm;
  double M2average = average[3]*norm;
  double Mabsaverage = average[4]*norm;
  // all expectation values are per spin, divide by 1/n_spins/n_spins
  double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
  double Mvariance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << mcs;
  ofile << setw(15) << setprecision(8) << temperature;
  ofile << setw(15) << setprecision(8) << Eaverage/n_spins/n_spins;
  //to devided by n_spins twice is just as deviding by n_spinsxn_spins
  //so the energy avr stored is of each dipole
  ofile << setw(15) << setprecision(8) << Evariance/temperature/temperature;
  ofile << setw(15) << setprecision(8) << Maverage/n_spins/n_spins;
  ofile << setw(15) << setprecision(8) << Mvariance/temperature;
  ofile << setw(15) << setprecision(8) << Mabsaverage/n_spins/n_spins ;
  ofile << setw(15) << setprecision(8) << E ;
  ofile << setw(15) << setprecision(8) << flips << endl;

} // end output function



