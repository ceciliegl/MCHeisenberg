#ifndef SOLVER
#define SOLVER

#include <fftw3.h>
#include "vec3.hpp"
#include "cvec3.hpp"
//#include "myrandom.hpp"
#include "Lattice.hpp"
#include "BasicFunctions.hpp"

class Solver
{
public:
  string dir;
  int Nh, Ns;
  int N;
  int NN, NsNs;
  int Nq;

  int holepos;

  int ncyc, nequ, nbins, nconew;
  int dnconew;
  double EQUENERGY;
  int print;
  int configprint;

  double conew;

  Lattice mylattice;

  double beta;

  vector<vec3> spins; //Save the spins in a vector of length 3N? doubles? between -1 and 1.

  double E;
  vec3 M;

  int naccepts;
  int currentaccepts;

  const int nexpvals = 8;
  vector<vector<double>> expectationvalues;
  vector<vector<double>> expOP;
  #ifdef COMPUTESQ
    vector<vector<cvec3>>  expSq;
  #endif

  vector<int> nnbrs;
  vector<vector<int>>    neighbours;
  vector<vector<double>> couplings;   //Currently saving the couplings as doubles. Could optionally store index of a J vector J={J1,J2,...} to save space.
  double DXY;

  RandomNumberGenerator MYRAN;

  ofstream Runlog;
  ofstream OutfileEnergy;
  ofstream OutfileExpect;
  ofstream OutfileOP;
  #ifdef COMPUTESQ
    ofstream OutfileExpectSq;
  #endif
  ofstream OutfileExpectLog;
  #ifdef COMPUTESQ
    ofstream OutfileExpectSqLog;
  #endif
  ofstream OutfileAcceptanceLog;
  ofstream OutfileRvecs;
  ofstream OutfileConfigs;
  ofstream OutfileEquConfigs;

  Solver();
  Solver(string dir0, Lattice mylattice, ReadInputFiles params);

  void solve(double beta0);

  void initialize_random();

  void calc_energy();
  void calc_magnetisation();

  vec3 localfield(int site);
  void Equilibrate();

  void MonteCarloUpdate();
  void Metropolis(); //Choose random spin, choose random direction, rotate pi around this direction. Compute deltaE and do Metropolis.
  void conew_update(double R); //Update the cone width for Gaussian Metropolis.
  void OverRelaxation(); //Rotate by pi about local field.

  void print_Rvecs();
  void print_config();
  void print_equconfig();

  //I want to compute the spin spin correlation function.
  //How to do that? FFT of the spins? And then compute expectation value of SqSq? I guess.
  //Or do I want the spatial ones. Do q-space for now.

  vector<cvec3> Sq();
  complex<double> OP(int subi);

  void resetdatafiles(vector<double> betas);
};

Solver::Solver(){}

Solver::Solver(string dir0, Lattice mylattice0, ReadInputFiles params)
{
  MYRAN = RandomNumberGenerator(params.RANSEED); //USE THIS INSTEAD OF THE GLOBAL RAN
  //MYRAN = RAN;
  dir = dir0;
  Nh = params.Nh;

  ncyc  = params.ncyc;
  nequ  = params.nequ;
  nbins = params.nbins;
  nconew = params.nconew;
  EQUENERGY = params.EQUENERGY;
  print = params.print;
  configprint = params.configprint;

  #ifdef ADAPTIVEGAUSSIAN
    dnconew = nequ/nconew;
  #endif

  cout << "Nbins = " << nbins << endl;

  mylattice = mylattice0;

  nnbrs = mylattice.nnbrs;

  neighbours = mylattice.neighbours;
  couplings  = mylattice.couplings;
  DXY = mylattice.DXY; //Single ion anisotropy, let it point in the z-direction, if D is positive, XY is favoured.

  N = mylattice.Ntot; //Number of sites.
  NN = N*N;
  Nq = mylattice.Nq;

  holepos = params.holepos;

  Ns = N - Nh; //Number of spins;
  NsNs = Ns*Ns;

  spins = vector<vec3>(N, 0.0); //The spin vector should also have the hole sites, where the spins are set to zero. You could also have set the bonds to zero for these sites.

  expectationvalues = vector<vector<double>>(nexpvals, vector<double>(nbins, 0.0));
  expOP = vector<vector<double>>(mylattice.subn, vector<double>(nbins, 0.0));
  #ifdef COMPUTESQ
    expSq = vector<vector<cvec3>>(Nq, vector<cvec3>(nbins, cvec3(zero, zero, zero)));
  #endif

  conew = 60.0; //Initial cone width for Gaussian Metropolis.

  Runlog = ofstream(dir + "runlog.txt", std::ios_base::app);
  if (!Runlog.is_open())
   cout<< "Could not open file." << endl;

  #ifdef ENERGYMINIMISATION
    Runlog << "RUNNING ENERGY MINIMISATION" << endl;
    Runlog << endl;
  #else
    #ifdef ADAPTIVEGAUSSIAN
      Runlog << "RUNNING ADAPTIVE GAUSSIAN" << endl;
      Runlog << endl;
    #else
      Runlog << "RUNNING METROPOLIS" << endl;
      Runlog << endl;
    #endif
  #endif

  OutfileEnergy = ofstream(dir + "energy.txt", std::ios_base::app);
  if (!OutfileEnergy.is_open())
    cout<< "Could not open file." << endl;

  OutfileExpect = ofstream(dir + "expvals.txt", std::ios_base::app);
  if (!OutfileExpect.is_open())
    cout<< "Could not open file." << endl;

  OutfileOP = ofstream(dir + "OP.txt", std::ios_base::app);
  if (!OutfileOP.is_open())
    cout<< "Could not open file." << endl;

  #ifdef COMPUTESQ
    OutfileExpectSq = ofstream(dir + "expSq.txt", std::ios_base::app);
    if (!OutfileExpectSq.is_open())
      cout<< "Could not open file." << endl;
  #endif

  OutfileExpectLog = ofstream(dir + "expvalslog.txt", std::ios_base::app);
  if (!OutfileExpectLog.is_open())
    cout<< "Could not open file." << endl;

  #ifdef COMPUTESQ
    OutfileExpectSqLog = ofstream(dir + "expSqlog.txt", std::ios_base::app);
    if (!OutfileExpectSqLog.is_open())
      cout<< "Could not open file." << endl;
  #endif

  OutfileAcceptanceLog = ofstream(dir + "acceptancelog.txt", std::ios_base::app);
  if (!OutfileAcceptanceLog.is_open())
    cout<< "Could not open file." << endl;

  OutfileRvecs = ofstream(dir + "Rvecs.txt", std::ios_base::app);
  if (!OutfileRvecs.is_open())
    cout<< "Could not open file." << endl;

  #ifdef EQUPRINT
    OutfileEquConfigs = ofstream(dir + "equconfigs.txt", std::ios_base::app);
    if (!OutfileEquConfigs.is_open())
      cout<< "Could not open file." << endl;
  #endif

  /*for(int i = 0; i < N; i++)
  {
    cout << "Site " << i << ":" << endl;

    for(int j = 0; j < nnbrs[i]; j++) cout << "coupling " << couplings[i][j] << " to site " << neighbours[i][j] << endl;
  }*/
}

void Solver::solve(double beta0)
{
  beta = beta0;

  OutfileConfigs = ofstream(dir + "configs_" + to_string(beta) + ".txt", std::ios_base::app);
  if (!OutfileConfigs.is_open())
    cout<< "Could not open file." << endl;

  expectationvalues = vector<vector<double>>(nexpvals, vector<double>(nbins, 0.0));
  expOP = vector<vector<double>>(mylattice.subn, vector<double>(nbins, 0.0));
  #ifdef COMPUTESQ
    expSq = vector<vector<cvec3>>(Nq, vector<cvec3>(nbins, cvec3(zero, zero, zero)));
  #endif

  for(int ibin = 0; ibin < nbins; ibin++)
  {
    try 
    {
      string filename = "initspins.txt";
      // Try to read in init spins from file.
      ifstream initfile(dir + filename);
      if (!initfile.is_open())
      {
        throw filename; // Throw an exception when a problem arise
      }
      for (int i = 0; i < N; i++)
      {
        initfile >> spins[i].x;
        initfile >> spins[i].y;
        initfile >> spins[i].z;
      }
      initfile.close();
      cout<< "Initialized from file." << endl;
    }
    catch (string filename) 
    {
      // If init spins are not read in, set them to random values.
      initialize_random();
      cout<< "Initialized random." << endl;
    }

    beta = beta0;

    Equilibrate(); //Run nequ MC cycles before starting to store data.

    calc_energy();
    calc_magnetisation();

    //cout << "M = " << M << "   M2 = " << dot(M,M) << endl;

    naccepts = 0;

    for(int cycle = 0; cycle < ncyc; cycle++)
    {
      //cout << cycle << endl;
      MonteCarloUpdate();

      expectationvalues[0][ibin] += E/N;
      expectationvalues[1][ibin] += E*E/NN;
      expectationvalues[2][ibin] += M.x/N;
      //cout << M.x/N << "   " << M.x << "   " << N << endl;
      expectationvalues[3][ibin] += M.y/N;
      expectationvalues[4][ibin] += M.z/N;
      //cout << "M = " << M << "   M2 = " << dot(M,M) << endl;
      expectationvalues[5][ibin] += dot(M,M)/NN;
      expectationvalues[6][ibin] += sqrt(dot(M,M))/N;
      expectationvalues[7][ibin] += M.z*M.z/NN;

      /*if(mylattice.latticetype == "PYROCHLORE")
      {
        for(int isub = 0; isub < mylattice.subn; isub++)
        {
          expOP[isub][ibin] += abs(OP(isub));
          //cout << expOP[isub][ibin] << "     ";
        }
        //cout << endl;
      }*/

      #ifdef COMPUTESQ
        vector<cvec3> mySqs = Sq();
        for(int iq = 0; iq < Nq; iq++)
        {
          expSq[iq][ibin].x += conj(mySqs[iq].x)*mySqs[iq].x; //Divide by NN?
          expSq[iq][ibin].y += conj(mySqs[iq].y)*mySqs[iq].y;
          expSq[iq][ibin].z += conj(mySqs[iq].z)*mySqs[iq].z;
        }
      #endif

      if ((cycle+1)%print == 0) //Write data to file
      {
        //cout << cycle << endl;
        OutfileEnergy << E << "\n";
        OutfileExpectLog << beta << endl;
        OutfileExpectLog << cycle+1;
        for(int i = 0; i < nexpvals; i++) OutfileExpectLog << " " << expectationvalues[i][ibin]/(cycle+1);
        OutfileExpectLog << "\n";
        #ifdef COMPUTESQ
          for(int iq = 0; iq < Nq; iq++)
          {
            OutfileExpectSqLog << beta << endl;
            OutfileExpectSqLog << cycle+1 << endl;
            OutfileExpectSqLog << setprecision(10) << mylattice.Qvecs[iq] << "   " << expSq[iq][ibin].x.real()/(cycle+1) << "   " << expSq[iq][ibin].x.imag()/(cycle+1)
                                                                          << "   " << expSq[iq][ibin].y.real()/(cycle+1) << "   " << expSq[iq][ibin].y.imag()/(cycle+1)
                                                                          << "   " << expSq[iq][ibin].z.real()/(cycle+1) << "   " << expSq[iq][ibin].z.imag()/(cycle+1)
                                                                          << "\n";
          }
        #endif
        //outfile_accepts << count_accepts << "  " << (cycle+1)*NN<< "\n";

        OutfileAcceptanceLog << "beta: " << beta << "    bin: " << ibin << "    conew: " << conew << "    cycle: " << cycle+1 << "    %accept: " << double(currentaccepts)/(N)*100 << endl;
      }
    }
    Runlog << "beta: " << beta << "    bin: " << ibin << "    conew: " << conew << "    %accept: " << double(naccepts)/(N*ncyc)*100 << endl;
  }

  vector<double> expectationvalues_mean(nexpvals, 0.0);
  vector<double> expectationvalues_std(nexpvals, 0.0);

  vector<double> expOP_mean(mylattice.subn, 0.0);
  vector<double> expOP_std(mylattice.subn, 0.0);

  #ifdef COMPUTESQ
    vector<cvec3> expSq_mean(Nq, cvec3(0,0,0));
    vector<cvec3> expSq_std(Nq, cvec3(0,0,0));
  #endif

  cout << "Ncyc = " << ncyc << endl;

  for(int j = 0; j < nexpvals; j++)
  {
    vector<double> ans = mean_stdev(expectationvalues[j]);
    expectationvalues_mean[j] = ans[0]/ncyc;
    expectationvalues_std[j]  = ans[1]/ncyc;
  }

  for(int j = 0; j < mylattice.subn; j++)
  {
    vector<double> ans = mean_stdev(expOP[j]);
    expOP_mean[j] = ans[0]/ncyc;
    expOP_std[j]  = ans[1]/ncyc;
  }

  #ifdef COMPUTESQ
    for(int j = 0; j < Nq; j++)
    {
      vector<cvec3> ans = mean_stdev(expSq[j]);
      expSq_mean[j] = ans[0]/ncyc;
      expSq_std[j]  = ans[1]/ncyc;

    }
  #endif

  OutfileExpect << beta << " " << expectationvalues_mean[0] << " " << expectationvalues_std[0]
                        << " " << expectationvalues_mean[1] << " " << expectationvalues_std[1]
                        << " " << expectationvalues_mean[2] << " " << expectationvalues_std[2]
                        << " " << expectationvalues_mean[3] << " " << expectationvalues_std[3]
                        << " " << expectationvalues_mean[4] << " " << expectationvalues_std[4]
                        << " " << expectationvalues_mean[5] << " " << expectationvalues_std[5]
                        << " " << expectationvalues_mean[6] << " " << expectationvalues_std[6]
                        << " " << expectationvalues_mean[7] << " " << expectationvalues_std[7]
                        << endl;

  OutfileOP << beta;
  for(int subi = 0; subi < mylattice.subn; subi++) OutfileOP << " " << expOP_mean[subi] << " " << expOP_std[subi] << " ";
  OutfileOP << endl;

  #ifdef COMPUTESQ
    OutfileExpectSq << beta << endl;
    for(int iq = 0; iq < Nq; iq++)
    {
      OutfileExpectSq << setprecision(10) << mylattice.Qvecs[iq] << "   " << expSq_mean[iq].x.real() << "   " << expSq_std[iq].x.real() << "   " << expSq_mean[iq].x.imag() << "   " << expSq_std[iq].x.imag()
                                            << "   " << expSq_mean[iq].y.real() << "   " << expSq_std[iq].y.real() << "   " << expSq_mean[iq].y.imag() << "   " << expSq_std[iq].y.imag()
                                            << "   " << expSq_mean[iq].z.real() << "   " << expSq_std[iq].z.real() << "   " << expSq_mean[iq].z.imag() << "   " << expSq_std[iq].z.imag()
                                            << endl;
    }
  #endif

  if (configprint == 1) print_config();

  OutfileConfigs.close();

  cout << "Done with beta = " << beta << "." << endl;
}

void Solver::initialize_random()
{
  for(int i = 0; i < N; i++)
  {
    spins[i] = randomunitvector(RAN);
  }

  if(Nh == 1)//Let holepos be zero if Nh == 1
  {
    spins[holepos] = vec3(0.0,0.0,0.0);
  }
}

void Solver::calc_energy()
{
  E = 0;
  //Heisenberg interactions:
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < nnbrs[i]; j++)
    {
      E += couplings[i][j]*dot(spins[i],spins[neighbours[i][j]]);
    }
    E += 2*DXY*spins[i].z*spins[i].z; //Multiply by 2 as we divide by 2 later.
  }
  E *= 0.5;
}

void Solver::calc_magnetisation()
{
  M = vec3(0,0,0);

  for(int i = 0; i < N; i++) M = M + spins[i];
}

vec3 Solver::localfield(int site)
{
  vec3 locfield(0,0,0);
  for(int j = 0; j < nnbrs[site]; j++) locfield = locfield + couplings[site][j]*spins[neighbours[site][j]];
  return locfield;
}

void Solver::Equilibrate()
{

  calc_energy();
  cout << "E = " << E << endl;
  calc_magnetisation();

  #ifdef ENERGYMINIMISATION
  int i = 0;
    while(i < nequ) // !& E > EQUENERGY)
    {
      MonteCarloUpdate();
      i++;
    }
    //This part could be improved to restart in a new random state if equ not reached:
    if (E > EQUENERGY) {Runlog << "DID NOT REACH ENERGY BELOW " << setprecision(17) << EQUENERGY << " AFTER " << nequ << " EQUILIBRATION STEPS! E = " << setprecision(17) << E << endl; Runlog << "STARTING EQUILIBRATION FROM NEW RANDOM CONFIGURATION." << endl; initialize_random(); Equilibrate();}
    else {Runlog << "EQUILIBRATION REACHED AFTER " << i << " STEPS, WITH ENERGY " << setprecision(17) << E << " EQUENERGY = " << setprecision(17) << EQUENERGY << endl;}
  #else
    double bstart;
    #ifdef ADIABATIC
      bstart = 0.1;
      Runlog << "Equilibrating adiabatically, starting from beta = " << bstart << "." << endl;
    #else
      bstart = beta;
      Runlog << "Equilibrating at constant beta." << endl;
    #endif

    #ifdef ADIABATICFOCUS
      double bstart2 = 20;
      double bend2 = 100;
      if (beta < bstart2) 
      {
        bstart2 = 0.1;
        bend2 = beta;
      }

      bool betaless = false;
      if (beta < bend2) {bend2 = beta; betaless = true;}

      int ninit = 1000*(nequ/1000 > 2) + int(0.01*nequ)*(nequ/1000 <= 2);

      double db1, db2, db3;
      db1 = (bstart2-bstart)/(ninit);
      if (betaless) db2 = (bend2-bstart2)/(int(0.9*nequ)-1-ninit);
      else db2 = (bend2-bstart2)/(int(0.8*(0.9*nequ-ninit)));
      db3 = (beta-bend2)/(int(0.9*nequ) - int(0.8*(0.9*nequ-ninit))-1);
    #else
      double db = (beta-bstart)/(int(0.9*nequ)-1);
    #endif

    #ifdef EQUPRINT
      vector<double> printbetas = {1, 2, 5, 10, 20, 25, 30, 35, 40, 45, 50, 55, 60, 70, 80, 90, 100, 200, 500, 1000, 5000, 10000, beta};
      sort(printbetas.begin(), printbetas.end());
      int ip = 0;
    #endif

    int nconewupdates = 0;
    currentaccepts = 0;
    for (int i = 0; i < int(0.9*nequ); i++)
    {
      #ifdef ADIABATICFOCUS
        if (i < ninit) beta = bstart + db1*i;
        else if (betaless) beta = bstart2 + db2*(i-ninit);
        else if (i < int(0.8*(0.9*nequ-ninit))) beta = bstart2 + db2*(i-ninit);
        else beta = bend2 + db3*(i-int(0.8*(0.9*nequ-ninit)));
      #else
        beta = bstart + db*i;
      #endif

      //cout << beta << endl;

      MonteCarloUpdate();

      #ifdef ADAPTIVEGAUSSIAN
        if ((i+1)%dnconew == 0)
        {
          conew_update(double(currentaccepts)/(N*dnconew));
          currentaccepts = 0;
          //cout << conew << endl;
          nconewupdates++;
        }
      #endif

      #ifdef EQUPRINT
        if (beta >= printbetas[ip]) 
        {
          print_equconfig();
          Runlog << "Printing equilibration configuration at beta = " << beta << " and conew = " << conew << "." << endl;
          ip++;
        }
      #endif
    }
    for (int i = int(0.9*nequ); i < nequ; i++)
    {
      MonteCarloUpdate();
      #ifdef ADAPTIVEGAUSSIAN
        if ((i+1)%dnconew == 0)
        {
          conew_update(double(currentaccepts)/double(N*dnconew));
          currentaccepts = 0;
          //cout << conew << endl;
          nconewupdates++;
        }
      #endif
    }
    #ifdef ADAPTIVEGAUSSIAN
      Runlog << "Equilibration: " << nequ << " steps, " << nconewupdates << " conew updates." << endl;
    #else
      Runlog << "Equilibration: " << nequ << " steps." << endl;
    #endif
  #endif
}

void Solver::MonteCarloUpdate()
{
  //Do N MC updates.

  for(int iMC = 0; iMC < N; iMC++)
  {
    Metropolis();
    OverRelaxation();
  }
}

void Solver::Metropolis()
{
  int site = int(MYRAN.RanInt(N-1-Nh)+Nh+holepos)%N; //Is this correct? //Draw a random spin, but not spin at holeposition if Nh = 1.

  if(site >= N || site < 0) cout << "Site " << site << endl;
  if(Nh ==1 && site == holepos) cout << "Site " << site << endl;

  vec3 newspin;

  #ifdef ADAPTIVEGAUSSIAN
    newspin = spins[site] + conew*randomGaussianvector(MYRAN); //Update spin according to Gaussian moves.
    newspin = newspin/sqrt(dot(newspin,newspin)); //Normalize the vector.
  #else
  #ifdef ENERGYMINIMISATION
    vec3 locf = localfield(site); //Compute the local field at the site.
    newspin = (-1)*locf/sqrt(dot(locf,locf)); //Set the spin antiparallel to the local field.
  #else
    newspin = randomunitvector(MYRAN); //Draw a random unit vector.
  #endif
  #endif

  double dE = dot(localfield(site),(newspin-spins[site])); //Is this correct? Or wrong by a factor of 2? //Compute the energy change.
  dE += DXY*(newspin.z*newspin.z - spins[site].z*spins[site].z); //Single ion anisotropy.

  #ifdef ENERGYMINIMISATION
  bool crit = dE < 0;
  #else
  bool crit = (dE < 0) || (MYRAN.Get0inc1() <= exp(-beta*dE));
  #endif

  if(crit)
  {
    //cout << "Accepting." << endl;
    E += dE;
    M = M + (newspin - spins[site]);
    spins[site] = newspin;
    currentaccepts++;
    naccepts++;
  }
}

void Solver::conew_update(double R) //R is the acceptance ratio in the previous MC step.
{
  double f = 0.5/(1.-R);
  conew *= f;
  //cout << f << endl;
  if (conew > 60) conew = 60.0;
  cout << "Acceptance ratio: " << R << "   sigma: " << conew << endl;
}

void Solver::OverRelaxation()
{
  int site = int(MYRAN.RanInt(N-1-Nh)+Nh+holepos)%N; //Is this correct? //Draw a random spin, but not spin at holeposition if Nh = 1.
  vec3 locfield = localfield(site);

  //cout << "Loc: " << locfield << endl;

  M = M - spins[site];

  spins[site] = RotPI(spins[site], locfield, MYRAN);

  M = M + spins[site];
}

vector<cvec3> Solver::Sq() //Returns SqSq for each q in mylattice.Qvecs. Maybe do {SqxSqx, SqySqy, SqzSqz}? Not for now.
{
  vector<vec3> SqSqvec(Nq, vec3(0,0,0)); //Should have real elements.
  vector<cvec3> Sqvec(Nq, cvec3(zero,zero,zero));  //Elements are in general complex.

  if(mylattice.latticetype == "PYRO16" || mylattice.latticetype == "TRIANGULARLADDER" || mylattice.latticetype == "TRIANGULAR" || mylattice.latticetype == "PYROCHLORE" || mylattice.latticetype == "KAGOME")
  {
    for(int iq = 0; iq < Nq; iq++)
    {
      vec3 q = mylattice.Qvecs[iq];
      for(int i = 0; i < N; i++) Sqvec[iq] = Sqvec[iq] + spins[i]*exponential(dot(q, mylattice.Rvecs[i]));
      Sqvec[iq] = Sqvec[iq]/complex<double>(sqrt(N));
    }
  }
  else
  {
    //For the cases where we can use FFTW:

    vector<complex<double>> Svec(N*3, zero);
    int dim = 3;
    vector<int> n = {mylattice.Nx, mylattice.Ny, mylattice.Nz};

    //MAYBE MORE EFFICIENT TO ADD THIS PLAN ONCE IN THE BEGINNING.
    fftw_plan SPINS_RTOQ = fftw_plan_many_dft(dim, &n[0], 3*mylattice.subn, reinterpret_cast<fftw_complex*>(&Svec[0]),
                                              NULL, 3*mylattice.subn, 1, reinterpret_cast<fftw_complex*>(&Svec[0]),
                                              NULL, 3*mylattice.subn, 1, 1, FFTW_MEASURE);

    //Start by doing each sublattice and x y z separately? Or make a long vector and use stride.

    for(int i = 0; i < N; i++)
    {
      Svec[3*i+0] = spins[i].x;
      Svec[3*i+1] = spins[i].y;
      Svec[3*i+2] = spins[i].z;
    }

    //Then do FFT using stride.

    fftw_execute(SPINS_RTOQ);

    //Then you need to sum correctly over the alphas? Or multiply them in? Or do this outside? Do it here for now? No, because it will depend on the number of zones etc.
    //So if all qs in all zones are stored in Qvecs, we will treat it here. The mylattice.N first elements of Qvecs are the ones in the first BZ?

    //Also, this function is returning a cvec3. You want to return a vector<cvec3> instead with all the qs.

    vector<cvec3> Sqvec0(N, zero); //Number of q's corresponds to the number of unit cells.

    for(int i = 0; i < N; i++)
    {
      Sqvec0[i].x = Svec[3*i+0]/complex<double>(sqrt(N));
      Sqvec0[i].y = Svec[3*i+1]/complex<double>(sqrt(N));
      Sqvec0[i].z = Svec[3*i+2]/complex<double>(sqrt(N));
    }

    //Sqvec now has Sq for each sublattice and each q in the first BZ.
    //Lastly, we should convert to SqSq for each q in Qvecs:

    //How to do this correctly? Need to know which qs that belong to eachother in the different zones.
    //The easiest is to assume that the qs are ordered such that the mylattice.N first qs corresond to one zone,
    //the mylattice.N next to the next zone, and so on.

    for(int iq = 0; iq < Nq; iq++)
    {
      for(int j = 0; j < mylattice.subn; j++) Sqvec[iq] = Sqvec[iq] + Sqvec0[mylattice.subn*(iq%mylattice.N)+j]*exponential(dot(mylattice.Qvecs[iq],mylattice.subvecs[j]));
    }

    fftw_destroy_plan(SPINS_RTOQ);
  }

  return Sqvec;
}

complex<double> Solver::OP(int subi)
{
  //for(int j = subi; j < N; j++) cout << spins[j] << endl;
  complex<double> ans(zero);
  if(mylattice.latticetype == "PYROCHLORE")
  {
    vector<int> trispins;
    vector<complex<double>> w = {{1,0}, exponential(TWOPI/3), exponential(-TWOPI/3)};

    for(int j = subi; j < N; j+=mylattice.subn)
    {
      //cout << j << "   " << ans << endl;
      trispins = mylattice.GetTriangle(j);
      //cout << trispins[0] << "   " << trispins[1] << "   " << trispins[2] << endl;
      ans += w[0]*dot(spins[trispins[0]],spins[trispins[1]]) + w[1]*dot(spins[trispins[1]],spins[trispins[2]]) + w[2]*dot(spins[trispins[2]],spins[trispins[0]]);
    }
  }
  //cout << ans << endl;
  return ans*complex<double>(1./mylattice.N);
}

void Solver::print_Rvecs()
{
  for(int i = 0; i < N; i++)
  {
    OutfileRvecs << mylattice.Rvecs[i] << endl;
  }
}

void Solver::print_config()
{
  for(int i = 0; i < N; i++)
  {
    OutfileConfigs << spins[i] << "    ";
  }
  OutfileConfigs << endl;
}

void Solver::print_equconfig()
{
  OutfileEquConfigs << beta << endl;
  for(int i = 0; i < N; i++)
  {
    OutfileEquConfigs << spins[i] << "    ";
  }
  OutfileEquConfigs << endl;
}

void Solver::resetdatafiles(vector<double> betas)
{
  ofstream Runlogfile(dir + "runlog.txt");
  if (!Runlogfile.is_open())
     cout<<"Could not open file" << endl;

  ofstream Energyfile(dir + "energy.txt");
  if (!Energyfile.is_open())
    cout<<"Could not open file" << endl;

  ofstream Expectfile(dir + "expvals.txt");
  if (!Expectfile.is_open())
    cout<<"Could not open file" << endl;

  ofstream ExpectOPfile(dir + "OP.txt");
  if (!ExpectOPfile.is_open())
    cout<<"Could not open file" << endl;

  #ifdef COMPUTESQ
    ofstream ExpectSqfile(dir + "expSq.txt");
    if (!ExpectSqfile.is_open())
      cout<<"Could not open file" << endl;
  #endif

  ofstream Expectlogfile(dir + "expvalslog.txt");
  if (!Expectlogfile.is_open())
    cout<<"Could not open file" << endl;

  #ifdef COMPUTESQ
    ofstream ExpectSqlogfile(dir + "expSqlog.txt");
    if (!ExpectSqlogfile.is_open())
      cout<<"Could not open file" << endl;
  #endif

  ofstream OutfileAcceptanceLog(dir + "acceptancelog.txt");
  if (!OutfileAcceptanceLog.is_open())
    cout<< "Could not open file." << endl;

  ofstream OutfileRvecs(dir + "Rvecs.txt");
  if (!OutfileRvecs.is_open())
    cout<< "Could not open file." << endl;

  for (int i = 0; i < betas.size(); i++)
  {
    ofstream OutfileConfigs(dir + "configs_" + to_string(betas[i]) + ".txt");
    if (!OutfileConfigs.is_open())
      cout<< "Could not open file." << endl;
  }

  #ifdef EQUPRINT
    ofstream OutfileEquConfigs(dir + "equconfigs.txt");
    if (!OutfileEquConfigs.is_open())
      cout<< "Could not open file." << endl;
  #endif
}

#endif
