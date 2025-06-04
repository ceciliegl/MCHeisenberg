#include <fstream>
#include<string>

using namespace std;

class ReadInputFiles
{
public:
  string filename_parameters;

  string latticetype;

  int Nx, Ny, Nz;
  int Nh;

  int holepos;

  double J1, J2, J3, J3a, J3b;
  double DXY;

  int ncyc, nequ, nbins, nconew, print, configprint;

  unsigned long int RANSEED;

  double EQUENERGY;

  bool OBC;

  bool RESETFILES;

  ReadInputFiles(string filename_param);
  void generate();
};

ReadInputFiles::ReadInputFiles(string filename_param)
{
  filename_parameters = filename_param;
}



void ReadInputFiles::generate()
{
  ifstream parameters(filename_parameters);
  if (!parameters)
  {
    cerr << "Unable to open file " << filename_parameters << "." << endl;
    exit(1);   // call system to stop
  }

  string nothing;
  parameters >> nothing >> nothing >> latticetype;
  parameters >> nothing >> nothing >> Nx;
  parameters >> nothing >> nothing >> Ny;
  parameters >> nothing >> nothing >> Nz;
  parameters >> nothing >> nothing >> Nh;
  parameters >> nothing >> nothing >> holepos;
  parameters >> nothing >> nothing >> J1; parameters >> nothing >> nothing >> J2; parameters >> nothing >> nothing >> J3; parameters >> nothing >> nothing >> J3a; parameters >> nothing >> nothing >> J3b;
  parameters >> nothing >> nothing >> DXY;
  parameters >> nothing >> nothing >> ncyc; parameters >> nothing >> nothing >> nequ; parameters >> nothing >> nothing >> nbins; parameters >> nothing >> nothing >> nconew; 
  parameters >> nothing >> nothing >> print; cout << print << endl;
  parameters >> nothing >> nothing >> configprint;
  parameters >> nothing >> nothing >> RANSEED; cout << RANSEED << endl;
  parameters >> nothing >> nothing >> EQUENERGY; cout << EQUENERGY << endl;
  parameters >> nothing >> nothing >> OBC;
  parameters >> nothing >> nothing >> RESETFILES;

  parameters.close();
}
