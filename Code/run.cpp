// Write a code that does MC on the classical Heisenberg model. Add the possibility of a static hole.
#include<vector>
#include<complex>
#include<iostream>
#include <iomanip>
#include<string>
#include <time.h>

using namespace std;

#include "myrandom.hpp"

#include "ReadInputFiles.hpp"
#include "BasicFunctions.hpp"
#include "Solver.hpp"


int main(int argc, char const *argv[])
{
  double starttime = clock();

  string project = argv[1];
  string run_number;

  //double beta = 1;

  if(atoi(argv[2]) < 10)
  {
    run_number = std::string("00") + std::string(argv[2]);
  }
  else if(atoi(argv[2]) > 9 && atoi(argv[2]) < 100)
  {
    run_number = std::string("0") + std::string(argv[2]);
  }
  else if(atoi(argv[2]) > 99)
  {
    run_number = std::string(argv[2]);
  }

  string folder = std::string(getenv("HOME")) + std::string("/Documents/MCHeisenberg/Data/") + std::string(project) + std::string("/Run") + std::string(run_number) + std::string("/");

  cout << folder << endl;

  ReadInputFiles params(folder + "parameters.txt");
  params.generate();

  Lattice mylattice(params);

  #ifdef ENERGYMINIMISATION
  vector<double> betas = {0};
  #else
  //vector<double> betas = {20};
  vector<double> betas = {1001}; //, 2, 10, 20, 100, 200, 1000};
  #endif

  Solver mysolver(folder, mylattice, params);

  if(params.RESETFILES){mysolver.resetdatafiles(betas);}

  if(params.configprint > 0){mysolver.print_Rvecs();}
  
  for(int b = 0; b < betas.size(); b++) mysolver.solve(betas[b]);

  double endtime = clock();

  double tottime = (endtime-starttime)/CLOCKS_PER_SEC;

  mysolver.Runlog << "Total time    = " << tottime << endl;

  return 0;
}
