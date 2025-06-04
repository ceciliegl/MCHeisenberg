#ifndef BASICFUNCTIONS
#define BASICFUNCTIONS

#include <limits>
#include <numeric>
#include <algorithm>

#include "vec3.hpp"
//#include "myrandom.hpp"

const double PI=3.14159265358979323846264338;
const double TWOPI=2.*PI;
const double SQRTTHREEOVERTWO = 0.86602540378443864676;
const complex<double> zero(0.0,0.0);

using namespace std;

complex<double> exponential(double phi)
{
  //returns exp(i\phi)
  complex<double> a(cos(phi), sin(phi));
  return a;
}

vector<double> mean_stdev(vector<double> v)
{
  double sum = accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();

  double stdev = 0.;
  for (int i = 0; i < v.size(); i++) stdev += (v[i] - mean)*(v[i] - mean);

  stdev /= double(v.size()-1);

  double usikkerhet = sqrt(stdev/double(v.size())); //Usikkerhet i gjennomsnittet


  /*vector<double> diff(v.size());
  transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = sqrt(sq_sum / v.size());*/

  return {mean, usikkerhet};
}

vector<complex<double>> mean_stdev(vector<complex<double>> vcompl)
{
  vector<double> ansRe, ansIm;

  int n = vcompl.size();
  vector<double> v(n);

  //REAL PART:
  for(int i = 0; i < n; i++) v[i] = vcompl[i].real();
  ansRe = mean_stdev(v);

  //IMAGINARY PART:
  for(int i = 0; i < n; i++) v[i] = vcompl[i].imag();
  ansIm = mean_stdev(v);

  complex<double> totmean = {ansRe[0],ansIm[0]};
  complex<double> totstdev = {ansRe[1],ansIm[1]};

  return {totmean, totstdev};
}

vector<cvec3> mean_stdev(vector<cvec3> cv)
{
  //first convert all 6 elements of the cvec3 to separate vectors
  vector<vector<double>> v(6, vector<double>(cv.size(), 0.0));
  for(int i = 0; i < cv.size(); i++)
  {
    v[0][i] = cv[i].x.real();
    v[1][i] = cv[i].x.imag();
    v[2][i] = cv[i].y.real();
    v[3][i] = cv[i].y.imag();
    v[4][i] = cv[i].z.real();
    v[5][i] = cv[i].z.imag();
  }

  //convert back to cvec3
  cvec3 mean;
  cvec3 stdev;

  vector<double> ansreal;
  vector<double> ansimag;
  ansreal = mean_stdev(v[0]);
  ansimag = mean_stdev(v[1]);
  mean.x  = {ansreal[0], ansimag[0]};
  stdev.x = {ansreal[1], ansimag[1]};
  ansreal = mean_stdev(v[2]);
  ansimag = mean_stdev(v[3]);
  mean.y  = {ansreal[0], ansimag[0]};
  stdev.y = {ansreal[1], ansimag[1]};
  ansreal = mean_stdev(v[4]);
  ansimag = mean_stdev(v[5]);
  mean.z  = {ansreal[0], ansimag[0]};
  stdev.z = {ansreal[1], ansimag[1]};

  return {mean, stdev};
}

vec3 randomunitvector(RandomNumberGenerator& RANGEN)
{
  /*
  double phi      = TWOPI*RANGEN.Get();
  double costheta = 1-2*RANGEN.Get();
  double sintheta = sqrt(1-costheta*costheta);

  return vec3(cos(phi)*sintheta,sin(phi)*sintheta,costheta);
  */

  //Eventuelt:
  while(true)
  {
    double x = RANGEN.Get0inc1()-0.5;
    double y = RANGEN.Get0inc1()-0.5;
    double z = RANGEN.Get0inc1()-0.5;
    double r = sqrt(x*x+y*y+z*z);

    if(r > 0.5) continue;

    x /= r;
    y /= r;
    z /= r;

    return vec3(x,y,z);
  }
}

vec3 randomGaussianvector(RandomNumberGenerator& RANGEN)
{
  double x = RANGEN.RanNorm(0,1); //Draw a random number from a Gaussian distribution.
  double y = RANGEN.RanNorm(0,1); //Draw a random number from a Gaussian distribution.
  double z = RANGEN.RanNorm(0,1); //Draw a random number from a Gaussian distribution.
  
  return vec3(x, y, z);
}

vec3 RotPI(vec3 inspin, vec3 rotax, RandomNumberGenerator RANGEN)
{
  if(dot(rotax,rotax) == 0){rotax = randomunitvector(RANGEN); cout << "rotax = 0, set to random." << endl;}
  vec3 para = dot(inspin,rotax)/dot(rotax,rotax)*rotax;
  vec3 orto = inspin - para;

  vec3 rotspin = para - orto;

  return rotspin;
}

double findminimum(vector<double> A)
{
  double ans = numeric_limits<double>::max();
  for(int i = 0; i < A.size(); i++)
  {
    if(A[i] < ans) ans = A[i];
  }

  return ans;
}

#endif
