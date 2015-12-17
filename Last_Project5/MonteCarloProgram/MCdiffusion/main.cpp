/*
  1-dim random walk program which computes the 
  probability for being in a given position after
  a certain amount of time steps.
*/

//Trial is the number of particles that will be in the system.

#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
#include "gaussiandeviate.cpp"
using namespace  std;

// Function to read in data from screen, note call by reference
void initialise(int&, int&, double&, double &D_x) ;
// The Mc sampling for random walks
void  mc_sampling(int, int, double, int *, int *, int *, int*, int , int*, int);
// prints to screen the results of the calculations 
void  output(int, int, int *, int *, int *, int , int *, int , int *);

void zero(int * array, int len)
{
    for (int i = 0; i < len; i++)
    {
     array[i] = 0;
    }
}

int main()
{
  int max_trials, number_walks;
  double move_probability, D_x;
  int *density;

  initialise(max_trials, number_walks, move_probability, D_x) ;

  // Read in data
  int length=1/D_x;
  density = new int[length];
  int g_lengde=1*length;
  int * gausdensity= new int[g_lengde];
     zero(gausdensity, g_lengde);
  //initialze the density a
  for(int i=0;i< 1/D_x; i++)
  {
      density[i]=0;
  }

  int *walk_cumulative = new int [number_walks+1];
  int *walk2_cumulative = new int [number_walks+1];
  int *probability = new int [2*(number_walks+1)];
  for (int walks = 1; walks <= number_walks; walks++)
  {   
    walk_cumulative[walks] = walk2_cumulative[walks] = 0;
  }
  for (int walks = 0; walks <= 2*number_walks; walks++)
  {   
    probability[walks] = 0;
  } // end initialization of vectors
   // cout<< " before mc_sampling"<< endl;
  // Do the mc sampling  
  mc_sampling(max_trials, number_walks, move_probability, walk_cumulative, walk2_cumulative, probability, density, length, gausdensity, g_lengde);
  // Print out results 
 // cout<< " before output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl;
  output(max_trials, number_walks, walk_cumulative, walk2_cumulative, probability, length, density , g_lengde, gausdensity);
  //cout<< " after output------------------------------------------------"<< endl;
  delete [] walk_cumulative; // free memory
  delete [] walk2_cumulative; delete [] probability;
 cout<<"-END and RETURN::::::::::::::::::::::::::::::::::::::::::::::"<<endl;
  return 0; 
} // end main function


//  This is the sampling part--------------------------------------
void mc_sampling(int max_trials, int number_walks,
                 double move_prob, int *walk_cumulative,
                 int *walk2_cumulative, int *probability, int *density, int lengde,int *gausdensity, int g_lengde)
{
  double move_probability;

  long idum;
  //int position;
  idum=-1;  // initialise random number generator
  //New Code
  //Start position allways full nr of particles.
//int lengde=int(1/D_x);
       //cout<<"I m in mc_sampling 1: "<<endl;
  int * nxtdensity= new int[lengde];

  int * nxtgausdensity= new int[g_lengde];
   zero(nxtgausdensity, g_lengde);

   gausdensity[(g_lengde)-1]=0;
   gausdensity[0]=max_trials;
  int gaus_ran;
  //zerodensity=density; // density should be an array of zeros
  zero(nxtdensity, lengde);
  density[(lengde)-1]=0;
  density[0]=max_trials;
  gausdensity[0]=max_trials;
    for (int walks = 1; walks <= number_walks; walks++)
    {//gaussian method:
        for(int p=0; p<g_lengde;p++)
        {
            for(int particles=0; particles< gausdensity[p];particles++)
            {  //cout<<"#particles= "<< particles << endl;
                double real_gaus=gaussian_deviate(&idum);
               gaus_ran=int( real_gaus*1);
                // cout<< "real gaus: "<< real_gaus<< endl;
                //cout<< "gaus: "<< gaus_ran<< endl;
                //cout<<">>";
               if(gaus_ran>0)
                 {
                   //positive= positive+1;
                   //cout<< "gaus"<< gaus_ran<< endl;
                   //cout<<"positive: "<<positive<<endl;
                    //cout<<"Step -1- "<<positive<<endl;

                      if(p+gaus_ran<g_lengde)
                       { // cout<<"Step -3- "<<positive<<endl;
                          nxtgausdensity[p+gaus_ran] += 1;
                       }

                  }

               if(gaus_ran<0)
                   {
                   //cout<<"negative: "<<negative<<endl;

                   if(p+gaus_ran>0)
                       {  nxtgausdensity[p+gaus_ran] += 1;
                       }

                   }
               if(gaus_ran==0)
                   {
              nxtgausdensity[p] += 1;

                   }

            }
        }
        for (int i=0; i< g_lengde; i++)
        {gausdensity[i]=nxtgausdensity[i];
           // cout<<"density: "<<nxtdensity[i]<<endl;
        }
      //  cout<<"------------------------------------------"<<endl;
        zero(nxtgausdensity, g_lengde);
        gausdensity[(g_lengde)-1]=0;
        gausdensity[0]=max_trials;
     //cout<<"nr of walks bunn "<< walks<<endl;

        //END: Gaussian method.
       // cout<< "Should start writing out nxtdensity:--------------- "<< endl;

        //for (int i=0; i< lengde; i++)
        //{
         //   cout<< "zerodensity:++++++++++++++++++++++++++++++++++++++++++++ "<< nxtdensity[i]<< endl;
        //}
        //position=0;
    //  cout<<"nr of walks topp "<< walks<<endl;
        for(int position=0;position<lengde; position++)
        {
            for(int i=0; i<density[position]; i++)
                //applying this model give me a high wave at one step. and a "kick back" at the end.
                //the boundary state is to be reevaluated.
                //first this is just approximate model. secondly we do not have absolutt walles at the ends. But in the start we
                //have chambers so particles my difuse back and at the end they get absurde so a beeter modell will
                //be the following. Wi will allow o move beyond the border as we keep the probabiltiy % all the time and make sure
                //to just neglect all of them and keep suplying the % at the starting position.
//          {if (position<=0)
//            {             move_probability=1;
//            }
//            if (position >=lengde)
//            {                move_probability=0;
//            }
//            if(position > 0 and position < lengde )
//            {              move_probability=move_prob;
//            }
//            if (ran0(&idum) <= move_probability)
//            {  if(position<lengde)
//                {nxtdensity[position+1] += 1;
//                }
//            }
//            else
//            { if(position>0)
//                {  nxtdensity[position-1] += 1;
//                }
//            }
//          }
                // this model was better at we don't get the the highest values at the two ends, also we have solved
                //the previous problem
                //but we got and other problem resulting in to law value at the first step. I think we can solve that by making
                //the model even more realistic, by considering 2 movment of freedom. and equal probabiliy to move up down left
                // and right.  up and down will result in staying in it place.
            {
                                   move_probability=move_prob;

                       if (ran0(&idum) <= move_probability)
                       {  if(position<lengde)
                           {
                               nxtdensity[position+1] += 1;
                           }
                       }
                       else
                       { if(position>0)
                           {  nxtdensity[position-1] += 1;
                           }
                       }
               }

//Code where their is 25% for moving forward backward up and down.

//            {
 // move_probability=move_prob/2;
//                                 //  cout<<"move_probability: "<< move_probability;

//                       if (ran0(&idum) >= 3* move_probability)
//                       {  if(position<lengde)
//                           {nxtdensity[position+1] += 1;
//                           }
//                       }

//                       if ( move_probability < ran0(&idum) && ran0(&idum) <3* move_probability)
//                       {  if(position<lengde)
//                           {nxtdensity[position] += 1;
//                              // cout<<"nr in place= "<< nxtdensity[position] <<endl;
//                           }
//                       }

//                       if(ran0(&idum) <=  move_probability)
//                       { if(position>0)
//                           {  nxtdensity[position-1] += 1;
//                           }
//                       }
//               }
        }
        for (int i=0; i< lengde; i++)
        {density[i]=nxtdensity[i];
           // cout<<"density: "<<nxtdensity[i]<<endl;
        }
      //  cout<<"------------------------------------------"<<endl;
        zero(nxtdensity, lengde);
        density[(lengde)-1]=0;
        density[0]=max_trials;
     //cout<<"nr of walks bunn "<< walks<<endl;
    }  // end of loop over walks



 //END: ---------New Code---------
//  for (int trial=1; trial <= max_trials; trial++)
//  {
//       //cout<<"I m in mc_sampling 1: "<<endl;
//    int position = 0;

//    for (int walks = 1; walks <= number_walks; walks++)
//    {
//    //  cout<<"nr of walks topp "<< walks<<endl;
//        if (position<=0)
//         {
//            move_probability=1;
//          //  cout<<"I changed to move_probability=1 "<<endl;
//         }
//        if (position >=1/D_x)
//         {
//            move_probability=0;
//           // cout<<"I changed to move_probability=0 "<<endl;
//         }


//        if(position > 0 and position < 1/D_x )
//        {
//         move_probability=move_prob;
//        // cout<<"I m position > 0 ; I changed to move_probability=0.5 "<<endl;
//        }

//cout<<"---------------------------------------------------------------------"<<endl;
//      if (ran0(&idum) <= move_probability)
//      {
//         //  cout<<"pos-1: "<<position<< " probability= "<< move_probability<<endl;

//        position += 1;
//      }
//      else
//      {
//       //   cout<<"pos-1: "<<position<< " probability= "<< move_probability<<endl;
//        position -= 1;
//      }

//      walk_cumulative[walks] += position;
//      walk2_cumulative[walks] += position*position;
//      probability[position+number_walks] += 1;
//   //   cout<<"nr of walks bunn "<< walks<<endl;
//    }  // end of loop over walks
    //cout<<"Before density upgrade"<<endl;
    //cout<< "position: "<< position<< endl;

  //  density[position]+=1;
   //cout<<"density at ("<< position <<")= "<< density[position]<<endl;
    //cout<<"nr of trials "<< trial<<endl;
 // } // end of loop over trials

  //cout<<"Ferdig i mc_sampling "<<endl;
}   // end mc_sampling function  
//-------------------------------------------------------------------------------
// Writes the results to two different files
void output(int max_trials, int number_walks, 
            int *walk_cumulative, int *walk2_cumulative, int * probability, int length, int* density, int g_lengde, int* gausdensity)
{
     //cout<<"I m  in output... output output output output output output output!"<<endl;
     //cout<<"length of array density= "<< length<< endl;
    for  (int i=0; i< length; i++)
    {
        cout<<"density: "<< density[i]<<endl;
    }
    for  (int i=0; i< g_lengde; i++)
    {
        cout<<"Gaussian density: "<< "i: "<< i<< " --> "<< gausdensity[i]<<endl;
    }
   //output(max_trials, number_walks, walk_cumulative, walk2_cumulative, probability, length, density);
   // void  output(int, int, int *, int *, int *, int , int *);
  ofstream ofile("testwalkers.dat", fstream::in | fstream::out | fstream::trunc);
  //cout<<"Do u see me, Munir HELOOOOOOOOOOOOOOOOOOOOOOOOOO!"<<endl;
  ofstream probfile("probability.dat");
  //cout<<"Do u see me, Munir HELOOOOOOOOOOOOOOOOOOOOOOOOOO!"<<endl;
  ofstream denfile("density001_timesteps5000.dat");
  ofstream gausfile("gausdensity001_timesteps5000.dat");
  //cout<<"length of density@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2"<< endl;
  for(int i=0; i<length;i++)
  {   denfile<<i<<"  "<< density[i]<< endl;
  }
  denfile.close();
  for(int i=0; i<g_lengde;i++)
  {   gausfile<<i<<"  "<< gausdensity[i]<< endl;
  }
  gausfile.close();
  for( int  i = 1; i <=  number_walks; i++)
  {
    double xaverage = walk_cumulative[i]/((double) max_trials);
    double x2average = walk2_cumulative[i]/((double) max_trials);
    double variance = x2average - xaverage*xaverage;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(6) << i;
    ofile << setw(15) << setprecision(8) << xaverage;
    ofile << setw(15) << setprecision(8) << variance <<endl ;
     //   ofile << setw(15) << setprecision(8) << variance << endl;
  }
  ofile.close();
  // find norm of probability
  double norm = 0.;
  for( int  i = -number_walks; i <=  number_walks; i++)
  {
    norm += (double) probability[i+number_walks];
  }
  // write probability
  for( int  i = -number_walks; i <=  number_walks; i++)
  {
    double histogram = probability[i+number_walks]/norm;
    probfile << setiosflags(ios::showpoint | ios::uppercase);
    probfile << setw(6) << i;
    probfile << setw(15) << setprecision(8) << histogram << endl;
  }
  probfile.close();
 // cout<<"length of density@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2"<< endl;
}  // end of function output 

// Reads in data from screen
void initialise(int& max_trials, int& number_walks, double& move_probability, double& D_x )
{
  cout << "Chose spasial unit length(for project5 0.01 or 0.1)=  ";
  cin >> D_x;
  cout << "Number of Monte Carlo trials (particles in  the system) =";
  cin >> max_trials;
  cout << "Number of attempted walks=";
  cin >> number_walks;
  cout << "Move probability=";
  cin >> move_probability;
}  // end of function initialise   





