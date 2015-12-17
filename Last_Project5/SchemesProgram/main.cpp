#include <iostream>
#include <armadillo>
#include <fstream>
//#include "functions.h"

using namespace arma;
using namespace std;

// Function to read in data from screen
void read_input(double&, double&);

//initialization of v giving as v = u + us
vec initial(vec , double , int );
//---------------------
void ForwardEuler(vec v, double Delta_t, double Delta_x, int n_x, int n_t){
    double a, ui; vec vnew = zeros<vec>(n_x), V=zeros<vec>(n_x);
    a = Delta_t/(Delta_x*Delta_x) ;
    V = v;
    ofstream outfile;

    outfile.open("forwardEuler.txt"); // Writing results to file. Each line in file is
                                     // the spacial solution for time t.
    // algorithm
    for (int j=1; j<=n_t; j++)
    {

        outfile << V(0) + 1 << " "; // Writing in boundary conditions to file.
        for (int i=1; i<n_x-1; i++)
        {

            vnew(i) = a * V(i-1) + (1-2*a) * V(i) + a * V(i+1);
           // cout<<"vnew(i) "<<vnew(i) <<endl;
            //cout<<"V(i) "<<V(i) <<endl;
            //cout<<"v(i) "<<v(i) <<endl;
            // Stands for the sationary part analyticlay expressed as u_s(x)=1-x
            ui = vnew(i) + 1 - i*Delta_x;
            outfile << ui << " ";
        } // Writing in results to file.

        //outfile << V(n_x-1) << endl; // Writing in boundary conditions;
        outfile << V(n_x-1) << endl; // Writing in boundary conditions;
        V = vnew;
    }
    outfile.close();
}

//--------------------
vec tridiagonal(vec a1, vec a2, vec a3, vec b,int n_x){
    // Forward Substitution. Row reducing the matrix equation
    vec v = zeros<vec>(n_x);
    for (int i = 1; i < n_x-1; i++){
        float factor = a1(i)/a2(i-1);
        a2(i) = a2(i) - a3(i-1)*factor;
        b(i) = b(i) - b(i-1)*factor;
    }

    // Backward Substitution. Solving the equation for vector v.
    v(n_x-1) = b(n_x-1)/a2(n_x-1);
    for (int k = n_x-2; k >= 0; k--){
        v(k) = (b(k) - a3(k) * v(k+1))/a2(k);
    }
    return v;
}

void BackwardEuler(vec v, double Delta_t, double Delta_x, int n_x, int n_t){
    double a; vec A1 = zeros<vec>(n_x), A2 = zeros<vec>(n_x), A3 = zeros<vec>(n_x), V;

    V = v;
    cout<<"V(0)= "<< V(0)<< endl;
    a = Delta_t/(Delta_x*Delta_x) ;
    ofstream outfile;
    outfile.open("backwardEuler.txt");

    for (int i=0; i<n_x; i++){
        A1(i) = -a; A2(i) = 1 + 2*a; A3(i) = -a;
    }

    // Time loop
    for (int j=1; j<n_t; j++)
    {
        outfile <<  1 << " ";
        for (int i=1; i<n_x-1; i++)
        {
            outfile << V(i) + 1 - i*Delta_x << " ";
        }
         outfile << 0 << endl; // Writing in boundary conditions;
V(0)=0;
V(n_x-1)=0;
        V = tridiagonal(A1,A2,A3,V,n_x);
    }
    outfile.close();
}

void CrankNicolson(vec v, double Delta_t, double Delta_x, int n_x, int n_t){
    double a,a2,a3; vec A1 = zeros<vec>(n_x), A2 = zeros<vec>(n_x), A3 = zeros<vec>(n_x), vtilde, vnew;
    a = Delta_t / (Delta_x * Delta_x);
    a2 = 2 - 2*a;
    a3 = 2 + 2*a;
    ofstream outfile;

    outfile.open("crankNicolson.txt");
    // Setting the diagonals on of the LHS-matrix.
    for (int i=0; i<n_x; i++){
        A1(i) = -a; A2(i) = a3; A3(i) = -a;
    }

    vnew = v;
    cout<< "vnew(0)= "<< vnew(0)<< endl;
    vtilde = vnew;
    for (int j=1; j<n_t; j++){
        // Chaning the RHS vector v_old into vtilde.
        outfile <<  1 << " ";
        for (int i=1; i<n_x-1; i++)
        {
        outfile << vnew(i) + 1 - i*Delta_x << " ";
        }
        outfile << 0<< endl; // Writing in boundary conditions;

        for (int i=1; i<n_x-1; i++)
        {
            vtilde(i) = a*vnew(i-1) + a2*vnew(i) + a*vnew(i+1);
        }
    vtilde(0)=0;
    vtilde(n_x-1)=0;
        vnew = tridiagonal(A1,A2,A3,vtilde,n_x);
    }
    outfile.close();
}


vec Analytical(double Delta_t, double Delta_x, int n_x, int n_t)
{
    double pi = 3.14159265359;
    vec u;
    ofstream outfile;
    outfile.open("analytical.txt");

    for (int j=0; j <= n_t; j++)
    {
        u = zeros<vec>(n_x);
        for (int i=0; i<n_x; i++)
        {
            for (int n=1; n < 50; n++)
            {
            u(i) += -2/(n*pi)*sin(n*pi*i*Delta_x) * exp(-(n*n*pi*pi*j*Delta_t));
            }
            outfile << u(i) + 1 - i*Delta_x << " ";
        }
            outfile << endl;
    }
    outfile.close();
    return u;
}


//---------------------
int main()
{
        // Declaring variables and conditions on diffusal equation
    double L = 1, final_time, Delta_x, Delta_t;
    // Setting integration points and length of dimension
    int n_x , n_t, mcs;
    //double dx = L / Nx, dt = T / Nt;
    read_input(final_time, Delta_x );
    n_x=(L/Delta_x)+1;
    Delta_t=(Delta_x*Delta_x)/2;
    n_t=(final_time/Delta_t)+1;
    // Declaring vector v. Given as v = u + us, where u is true solution.
    vec v = zeros<vec>(n_x); // zeros call set boundary conditions.
    v=initial(v,Delta_x,n_x);

    //for (int i=1;i<n_x-1;i++){

     //cout<<"initialize vec v element: "<< v(i)<<  endl;
    //}
    ForwardEuler (v,Delta_t,Delta_x,n_x,n_t); // Changing v into spacial solution for t=T
    v=initial(v,Delta_x,n_x);
    BackwardEuler(v,Delta_t,Delta_x,n_x,n_t);
    v=initial(v,Delta_x,n_x);
    CrankNicolson(v,Delta_t,Delta_x,n_x,n_t);
    v = Analytical(Delta_t,Delta_x,n_x,n_t);

    return 0;
}

vec initial(vec v, double Delta_x, int n_x)
{
    v(0) = 0; v(n_x-1) = 0;
    for (int i=1;i<n_x-1;i++){
        v(i) = i*Delta_x - 1;
      //  cout<<"initialize vec v element: "<< v(i)<<  endl;
    }


return v;
}

// read in input data
void read_input( double& final_time, double& Delta_x)
{

  cout << "Give the spacial step length."<< endl;
  cout << "Note that the time integrration step is calculated according to D_t=(D_x)^2/2: ";
  cin >> Delta_x;
  cout << "Give the final time"<< endl;
  cout << "Note the program is on dimensionless basis, the final time should not excide 2: ";
  cin >> final_time;

} // end of function read_input

