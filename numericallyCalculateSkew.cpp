#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <omp.h>
/*
Benjamin M. Roberts
2018-02-05, 16:23

Calculates the expected skew numerically, for a range of x0 and r0t0 values.
Uses Gaussian monopole, and full standard halo model distributions.

Compile: g++ -fopenmp numericallyCalculateSkew.cpp -lm &&

Values are saved to a text file, in form: x0 r0t0 skew
Default filename: numericalSkew.txt
Note: will just over-ride the file if it already exists.

p_s(s) = r0t0*Integrate[ p_c(r)*p_x(s-r) , {r,-infty,infty} ]
        +(1-r0t0)*p_c(s)
where p_c is the 'intrinsic' clock noise (Gaussian).
p_chi is the DM signal distribution (in absence of noise). See paper.

*/

//------------------------------------------------------------------------------
double fv(double v)
/*
Standard halo model for velocity distribution, in laboratory frame.
 f ~ v^2 exp(-v^2)
Note: distribution for DM particles that cross paths with Earth.
We should have: <v> = 370
*/
{

  double ve = 544.; // galactic escape velocity
  double vc = 220.; // local frame velocity

  double Kn = 1.45303*sqrt(M_PI)*pow(vc,3);
  double A = pow(v,2)/Kn;

  double arg1 = -pow((v-vc)/vc,2);

  if(v<=0) return 0; //just for safety - should never be called with v<0

  if(v<=ve-vc){
    double arg2 = -pow((v+vc)/vc,2);
    return A*(exp(arg1)-exp(arg2));
  }else if(v<=ve+vc){
    double arg2 = -pow(ve/vc,2);
    return A*(exp(arg1)-exp(arg2));
  }else{
    return 0;
  }

}

//------------------------------------------------------------------------------
double g(double x, double x0, double s)
/*
Just a Gaussian distribution
*/
{
  double arg = -0.5*pow((x-x0)/s,2);
  double A = sqrt(2*M_PI)*s;
  return exp(arg)/A;
}

//------------------------------------------------------------------------------
double px(double x, double x0)
/*
p_chi for Gaussian monopole model:
*/
{
  if(x<0.001*x0) return 0;
  if(x>9.*x0) return 0;

  double e=2.71828;
  double b = 300.*x0*sqrt(M_PI)/(x);
  double a = b/e;

  double dx=0.8;
  double A=0;
  for(double t=a; t<=b; t+=dx){
    A+=fv(t);
  }
  return A*dx/x;
}

//------------------------------------------------------------------------------
double pxHS(double x, double x0)
/*
p_chi for uniform (hard) sphere model:
*/
{
  if(x<=0||x>x0) return 0.;
  return (2*x/pow(x0,2));
}


//******************************************************************************
int main(void){

  //Check SHM (v) normalisation:
  double dx=1.;
  double A=0;
  double A2=0;
  for(double x=0; x<800; x+=dx){
    A+=fv(x);
    A2+=x*fv(x);
  }
  std::cout<<"Check fv normalisation: "<<A*dx<<";  <v>="<<A2*dx<<"\n";

  //Check P_chi (and Gaussian) normalisation:
  A=0;
  A2=0;
  dx=0.01;
  for(double x=-10.; x<20.; x+=dx){
    A+=px(x,1.);
    A2+=g(x,0.,1.643);
  }
  std::cout<<"Check Px normalisation: "<<A*dx<<"; Gaus Norm: "<<A2*dx<<"\n";


  //define grid for the integral: pc(eta)*px(s-eta) d\eta
  double r0=8.;
  double dr=0.020;

  //Grid for the s integral [find the mean, s.d. and skew]
  double s0=10.;
  double ds=0.01;
  int ngp = (int) 2.*s0/ds;

  //Arrays to hold p_s(s) [Keep the r0t0 and (1-r0t0) parts seperate]
  std::vector<double> ps1(ngp),ps2(ngp);

  // Standard deviation of "intrinsic" clock noise.
  double sigma=1.; //always 1 - just sets the 'units'

  //Output file:
  std::ofstream ofile;
  std::string file_name="numericalSkew.txt";
  ofile.open(file_name.c_str());
  ofile<<"# x0     r0t0     skew\n";

  //Which values of x0 and r0t0 to calculate the skew for:
  int num_steps=10; //
  double x0min = 0.1, x0max = 3.;
  double r0min = 0.0001, r0max = 0.1;
  double dx0 = pow(x0max/x0min,1./(num_steps-1));
  double dr0 = pow(r0max/r0min,1./(num_steps-1));

  //Output info to screen:
  std::cout<<"\nCalculating skew numerically for:\n"
    <<" x0   = "<<x0min<<" -> "<<x0max<<"\n"
    <<" r0t0 = "<<r0min<<" -> "<<r0max<<"  (in "<<num_steps<<" steps)\n";
  std::cout<<"\nOutput filename: "<<file_name<<"\n\n";

  //Loop over chi_0:
  for(int ix=0; ix<num_steps; ix++){
    double x0 = x0min*pow(dx0,ix);
    std::cout<<"x0 = "<<x0<<".."<<std::flush;

    // First, calculate p_s(s)
    // (this is the loop that is slow)
    #pragma omp parallel for
    for(int i=0; i<ngp; i++){
      double s = i*ds - s0;
      double B=0;
      for(double r=-r0; r<r0; r+=dr){
        B+=px(s-r,x0)*g(r,0,sigma);
        //B+=pxHS(s-r,x0)*g(r,0,sigma); //Uniform sphere model
      }
      // Keep the r0t0 and (1-r0t0) parts seperate
      ps1[i]=B*dr;
      ps2[i]=g(s,0,sigma);
    }
    std::cout<<" .. ";

    //Loop over
    for(int ir=0; ir<num_steps; ir++){
      double r0t0 = r0min*pow(dr0,ir);
      double k0=0;  //Normalisation (not used)
      double k1=0;  //Mean
      double k2=0;  //S.deviation [actually, variance]
      double k3=0;  //Skewness* (before /s^3)
      for(int i=0; i<ngp; i++){
        double s = i*ds - s0;
        //k0 += r0t0*ps1[i] + (1.-r0t0)*ps2[i];   //just to check normalisation
        k1 += s*(r0t0*ps1[i] + (1.-r0t0)*ps2[i]);
      }
      k0*=ds; //should be 1 (not used, just to check!)
      k1*=ds;
      //calculate variance and skewness:
      for(int i=0; i<ngp; i++){
        double s = i*ds - s0;
        k2 += pow(s-k1,2)*(r0t0*ps1[i] + (1.-r0t0)*ps2[i]);
        k3 += pow(s-k1,3)*(r0t0*ps1[i] + (1.-r0t0)*ps2[i]);
      }
      k2*=ds;
      k3*=ds;
      double kappa3 = k3/pow(k2,1.5);
      //write results to file:
      ofile<<x0<<" "<<r0t0<<" "<<kappa3<<"\n";
    }//END loop over r0t0

    std::cout<<"..done\n";
    ofile<<"\n"; //seperate each 'block'
  }//END loop over x0

  ofile.close();

  return 0;
}
