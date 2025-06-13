#ifndef LATTICE
#define LATTICE

#include "vec3.hpp"
#include "BasicFunctions.hpp"



class Lattice
{
public:
  string latticetype;

  int Nx, Ny, Nz, N, Ntot, Nq;
  int subn, subsize;

  int Ncouplings;

  double J1, J2, J3, J3a, J3b;
  double DXY;

  vec3 a1, a2, a3, b1, b2, b3;
  vector<vec3> subvecs; //Sublattice vectors
  vector<vec3> Rvecs;
  vector<vec3> Qvecs;

  vector<int> nnbrs;
  vector<vector<int>>    neighbours;
  vector<vector<double>> couplings;
  //vector<double> SIanisotropy; Make it simple for now.

  Lattice();
  Lattice(ReadInputFiles params);

  vector<int> site_to_coord(int site);
  int coord_to_site(vector<int> coord);

  vector<int> GetTriangle(int site);
};

Lattice::Lattice(){}

Lattice::Lattice(ReadInputFiles params)
{
  latticetype = params.latticetype;
  //Må beregne alle q-vektorene?
  //Tenker å lage en vektor med q-koordinater.
  Nx = params.Nx; Ny = params.Ny; Nz = params.Nz;
  N = Nx*Ny*Nz;
  Nq = N;

  if(latticetype == "SQUARE" || latticetype == "SQUAREOBC" || latticetype == "CUBIC" || latticetype == "TRIANGULAR"
      || latticetype == "CHAIN" || latticetype == "TRIANGULARLADDER"){subn = 1;}
  else if(latticetype == "KAGOME"){subn = 3;}
  else if(latticetype == "PYROCHLORE"){subn = 4;}
  else if(latticetype == "PYRO16"){subn = 16;}
  else
  {
    cout << "Number of sublattices is not defined for chosen lattice type." << endl;
    exit(1);
  }

  Ntot = N*subn;

  int nZones = 1;

  if(latticetype == "PYRO16"){Ntot = N*subn; nZones = 4*4*4;}
  else if(latticetype == "PYROCHLORE"){Ntot = N*subn; nZones = 8;}
  else if(latticetype == "TRIANGULARLADDER"){Ny = 2; Nz = 1; N = Nx*Ny*Nz; Ntot = N; nZones = 1;}
  else if(latticetype == "KAGOME"){Ntot = N*subn; nZones = 4;}

  Nq = N*nZones;

  vector<vec3> extended(nZones);

  J1 = params.J1; J2 = params.J2; J3 = params.J3; J3a = params.J3a; J3b = params.J3b;
  DXY = params.DXY;

  int nnbrmax;

  if(latticetype == "CHAIN" || latticetype == "TRIANGULARLADDER") nnbrmax = 4;
  else if(latticetype == "SQUARE") nnbrmax = 12;
  else if(latticetype == "SQUAREOBC") nnbrmax = 12;
  else if(latticetype == "TRIANGULAR") nnbrmax = 6;
  else if(latticetype == "KAGOME") nnbrmax = 4;
  else if(latticetype == "PYRO16") nnbrmax = 6;
  else if(latticetype == "PYROCHLORE") nnbrmax = 30;
  else
  {
    cout << "NNBRMAX NOT DEFINED FOR " << latticetype << endl;
    exit(1);
  }

  nnbrs = vector<int>(Ntot, 0);
  neighbours = vector<vector<int>>(Ntot, vector<int>(nnbrmax, Ntot));
  couplings  = vector<vector<double>>(Ntot, vector<double>(nnbrmax, 0.0));

  //Må definere gittervektorene utifra hvilket gitter jeg vil ha.
  Rvecs = vector<vec3>(Ntot);
  Qvecs = vector<vec3>(Nq);

  subsize = subn*subn;
  subvecs = vector<vec3>(subn);
  subvecs[0] = vec3(0,0,0);

  if(latticetype == "CHAIN" || latticetype == "SQUARE" || latticetype == "SQUAREOBC" || latticetype == "CUBIC" || latticetype == "PYRO16")
  {
    a1 = vec3(1.,0,0);
    a2 = vec3(0,1.,0);
    a3 = vec3(0,0,1.);
    if(latticetype == "PYRO16")
    {
      subvecs[0]  = vec3(0,0,0);
      subvecs[1]  = vec3(0,0.25,0.25);
      subvecs[2]  = vec3(0.25,0,0.25);
      subvecs[3]  = vec3(0.25,0.25,0);
      subvecs[4]  = vec3(0,0.5,0.5);
      subvecs[5]  = subvecs[4] + subvecs[1];
      subvecs[6]  = subvecs[4] + subvecs[2];
      subvecs[7]  = subvecs[4] + subvecs[3];
      subvecs[8]  = vec3(0.5,0,0.5);
      subvecs[9]  = subvecs[8] + subvecs[1];
      subvecs[10] = subvecs[8] + subvecs[2];
      subvecs[11] = subvecs[8] + subvecs[3];
      subvecs[12] = vec3(0.5,0.5,0);
      subvecs[13] = subvecs[12] + subvecs[1];
      subvecs[14] = subvecs[12] + subvecs[2];
      subvecs[15] = subvecs[12] + subvecs[3];
    }
  }
  else if(latticetype == "TRIANGULAR" || latticetype == "TRIANGULARLADDER" || latticetype == "KAGOME")
  {
    a1 = vec3(1.,0,0);
    a2 = vec3(0.5,SQRTTHREEOVERTWO,0);
    a3 = vec3(0,0,1.);

    if(latticetype == "KAGOME")
    {
      subvecs[1] = 0.5*a1;
      subvecs[2] = 0.5*a2;
    }
  }
  else if(latticetype == "PYROCHLORE")
  {
    a1 = vec3(0,0.5,0.5);
    a2 = vec3(0.5,0,0.5);
    a3 = vec3(0.5,0.5,0);
    if(latticetype == "PYROCHLORE")
    {
      subvecs[1] = vec3(0,0.25,0.25);
      subvecs[2] = vec3(0.25,0,0.25);
      subvecs[3] = vec3(0.25,0.25,0);
    }
  }
  else
  {
    cout << "Lattice vectors are not defined for chosen lattice type." << endl;
    exit(1);
  }

  if(latticetype == "TRIANGULARLADDER")
  {
    int j = 0;

    for(int iz = 0; iz < Nz; iz++)
      for(int ix = 0; ix < Nx; ix++)
        for(int iy = 0; iy < Ny; iy++)
        {
          Rvecs[j] = ix*a1 + iy*a2; //Genererer Ntot r-vektorer.
          j++;
        }
  }
  else if(latticetype == "CHAIN")
  {
    int j = 0;

    for(int iz = 0; iz < Nz; iz++)
      for(int iy = 0; iy < Ny; iy++)
        for(int ix = 0; ix < Nx; ix++)
        {
          Rvecs[j] = ix*a1 + iy*a2 + iz*a3; //Genererer Ntot r-vektorer.
          j++;
        }
  }
  else
  {
    int j = 0;

    for(int iz = 0; iz < Nz; iz++)
      for(int iy = 0; iy < Ny; iy++)
        for(int ix = 0; ix < Nx; ix++)
          for(int isubn = 0; isubn < subn; isubn++)
          {
            Rvecs[j] = ix*a1 + iy*a2 + iz*a3 + subvecs[isubn]; //Genererer Ntot r-vektorer.
            j++;
          }
    cout << "WARNING: CHECK THAT R-VECS ARE CORRECT FOR " << latticetype << " but only needed if you're not using FFT." << endl;
  }

  /*
  for(int i = 0; i < Ntot; i++)
  {
    cout << Rvecs[i] << endl;
  }
  */


  //Calculate the reciprocal lattice vectors
  b1 = TWOPI*cross(a2,a3)/dot(a1, cross(a2, a3));
  b2 = TWOPI*cross(a3,a1)/dot(a2, cross(a3, a1));
  b3 = TWOPI*cross(a1,a2)/dot(a3, cross(a1, a2));

  //cout << b1.x << "    " << b1.y << "    " << b1.z << endl;
  //cout << b2.x << "    " << b2.y << "    " << b2.z << endl;
  //cout << b3.x << "    " << b3.y << "    " << b3.z << endl;

  if(latticetype == "PYRO16")
  {
    int count = 0;
    for(int i1 = 0; i1 < 4; i1++)
    for(int i2 = 0; i2 < 4; i2++)
    for(int i3 = 0; i3 < 4; i3++)
    {
      extended[count] = i1*b1+i2*b2+i3*b3;
      cout << setprecision(10) << extended[count] << endl;
      count++;
    }
  }
  else if(latticetype == "PYROCHLORE")
  {
    extended = {vec3(0,0,0), b1, b2, b3, b1+b2, b2+b3, b3+b1, b1+b2+b3};
  }
  else if(latticetype == "KAGOME")
  {
    extended = {vec3(0,0,0), b1, b2, b1+b2};
  }
  else if (nZones == 1)
  {
    extended = {vec3(0,0,0)};
  }
  else
  {
    cout << "EXTENDED ZONES NOT DEFINED FOR " << latticetype << endl;
    exit(1);
  }

  if(latticetype == "TRIANGULARLADDER")
  {
    int j = 0;

    double Nxinv = 1./Nx;
    for(int iz = 0; iz < Nz; iz++)
      for(int iy = 0; iy < Ny; iy++)
        for(int ix = 0; ix < Nx; ix++)
        {
          Qvecs[j] = vec3(ix*TWOPI*Nxinv, iy*PI/SQRTTHREEOVERTWO, 0);
          j++;
        }
  }
  else
  {
    int j = 0;

    double Nxinv = 1./Nx; double Nyinv = 1./Ny; double Nzinv = 1./Nz;
    for (int izone = 0; izone < nZones; izone++)
      for(int iz = 0; iz < Nz; iz++)
        for(int iy = 0; iy < Ny; iy++)
          for(int ix = 0; ix < Nx; ix++)
          {
            Qvecs[j] = ix*Nxinv*b1 + iy*Nyinv*b2 + iz*Nzinv*b3 + extended[izone]; //Genererer N q-vektorer.
            j++;
          }
  }

  //CURRENTLY DOUBLE COUNTING ALL BONDS, BUT DIVIDING BY TWO IN THE ENERGY.
  if(latticetype == "CHAIN" || latticetype == "TRIANGULARLADDER")
  {
    for(int i = 0; i < Ntot; i++)
    {
      neighbours[i][0] = (i-1+Ntot)%Ntot; couplings[i][0]  = J1;
      neighbours[i][1] = (i+1+Ntot)%Ntot; couplings[i][1]  = J1;

      nnbrs[i] += 2;

      neighbours[i][2] = (i-2+Ntot)%Ntot; couplings[i][2]  = J2;
      neighbours[i][3] = (i+2+Ntot)%Ntot; couplings[i][3]  = J2;

      nnbrs[i] += 2;
    }
  }
  else if(latticetype == "PYRO16")
  {
    /*neighbours[0]  = {1,2,3,5,10,15};
    neighbours[1]  = {0,2,3,4,11,14};
    neighbours[2]  = {0,1,3,7,8,13};
    neighbours[3]  = {0,1,2,6,9,12};
    neighbours[4]  = {1,5,6,7,11,14};
    neighbours[5]  = {0,4,6,7,10,15};
    neighbours[6]  = {3,4,5,7,9,12};
    neighbours[7]  = {2,4,5,6,8,13};
    neighbours[8]  = {2,7,9,10,11,13};
    neighbours[9]  = {3,6,8,10,11,12};
    neighbours[10] = {0,5,8,9,11,15};
    neighbours[11] = {1,4,8,9,10,14};
    neighbours[12] = {3,6,9,13,14,15};
    neighbours[13] = {2,7,8,12,14,15};
    neighbours[14] = {1,4,11,12,13,15};
    neighbours[15] = {0,5,10,12,13,14};*/

    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)

      //Nearest neighbours:
      if(coord[3] == 0)
      {
        neighbours[i][0] = coord_to_site({coord[0],     coord[1],   coord[2],     1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],     coord[1],   coord[2],     2}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],     coord[1],   coord[2],     3}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],     coord[1]-1, coord[2]-1,   5}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0]-1,   coord[1],   coord[2]-1,  10}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0]-1,   coord[1]-1, coord[2],    15}); couplings[i][5] = J1;
      }
      if(coord[3] == 1)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1], coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1], coord[2],   2}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1], coord[2],   3}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1], coord[2],   4}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0]-1, coord[1], coord[2],  11}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0]-1, coord[1], coord[2],  14}); couplings[i][5] = J1;
      }
      if(coord[3] == 2)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1],   coord[2],   1}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1],   coord[2],   3}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1]-1, coord[2],   7}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1],   coord[2],   8}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1]-1, coord[2],  13}); couplings[i][5] = J1;
      }
      if(coord[3] == 3)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1], coord[2],     0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1], coord[2],     1}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1], coord[2],     2}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1], coord[2]-1,   6}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1], coord[2]-1,   9}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1], coord[2],    12}); couplings[i][5] = J1;
      }
      if(coord[3] == 4)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1], coord[2],   1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1], coord[2],   5}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1], coord[2],   6}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1], coord[2],   7}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0]-1, coord[1], coord[2],  11}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0]-1, coord[1], coord[2],  14}); couplings[i][5] = J1;
      }
      if(coord[3] == 5)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1]+1, coord[2]+1,   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1],   coord[2],     4}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],     6}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1],   coord[2],     7}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],    10}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1,  15}); couplings[i][5] = J1;
      }
      if(coord[3] == 6)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1], coord[2]+1,   3}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1], coord[2],     4}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1], coord[2],     5}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1], coord[2],     7}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1], coord[2],     9}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1], coord[2]+1,  12}); couplings[i][5] = J1;
      }
      if(coord[3] == 7)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1]+1, coord[2],   2}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1],   coord[2],   4}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1],   coord[2],   5}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1],   coord[2],   6}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1]+1, coord[2],   8}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1],   coord[2],  13}); couplings[i][5] = J1;
      }
      if(coord[3] == 8)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1],   coord[2],   2}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1]-1, coord[2],   7}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1],   coord[2],   9}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1],   coord[2],  10}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1],   coord[2],  11}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1]-1, coord[2],  13}); couplings[i][5] = J1;
      }
      if(coord[3] == 9)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1], coord[2]+1,   3}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1], coord[2],     6}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1], coord[2],     8}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1], coord[2],    10}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1], coord[2],    11}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1], coord[2]+1,  12}); couplings[i][5] = J1;
      }
      if(coord[3] == 10)
      {
        neighbours[i][0] = coord_to_site({coord[0]+1, coord[1],   coord[2]+1,   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],     5}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],     8}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1],   coord[2],     9}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],    11}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1,  15}); couplings[i][5] = J1;
      }
      if(coord[3] == 11)
      {
        neighbours[i][0] = coord_to_site({coord[0]+1, coord[1], coord[2],   1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1], coord[2],   4}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1], coord[2],   8}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1], coord[2],   9}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1], coord[2],  10}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1], coord[2],  14}); couplings[i][5] = J1;
      }
      if(coord[3] == 12)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1], coord[2],     3}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1], coord[2]-1,   6}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1], coord[2]-1,   9}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1], coord[2],    13}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1], coord[2],    14}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1], coord[2],    15}); couplings[i][5] = J1;
      }
      if(coord[3] == 13)
      {
        neighbours[i][0] = coord_to_site({coord[0], coord[1]+1, coord[2],   2}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0], coord[1],   coord[2],   7}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0], coord[1]+1, coord[2],   8}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0], coord[1],   coord[2],  12}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0], coord[1],   coord[2],  14}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0], coord[1],   coord[2],  15}); couplings[i][5] = J1;
      }
      if(coord[3] == 14)
      {
        neighbours[i][0] = coord_to_site({coord[0]+1, coord[1], coord[2],   1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1], coord[2],   4}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1], coord[2],  11}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1], coord[2],  12}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1], coord[2],  13}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1], coord[2],  15}); couplings[i][5] = J1;
      }
      if(coord[3] == 15)
      {
        neighbours[i][0] = coord_to_site({coord[0]+1, coord[1]+1, coord[2],     0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1,   5}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1,  10}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1],   coord[2],    12}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],    13}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1],   coord[2],    14}); couplings[i][5] = J1;
      }

      nnbrs[i] += 6;
    }
  }
  else if(latticetype == "SQUARE")
  {
    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)

      neighbours[i][0]  = coord_to_site({coord[0]+1, coord[1],   coord[2],   coord[3]}); couplings[i][0]  = J1;
      neighbours[i][1]  = coord_to_site({coord[0]-1, coord[1],   coord[2],   coord[3]}); couplings[i][1]  = J1;
      neighbours[i][2]  = coord_to_site({coord[0],   coord[1]+1, coord[2],   coord[3]}); couplings[i][2]  = J1;
      neighbours[i][3]  = coord_to_site({coord[0],   coord[1]-1, coord[2],   coord[3]}); couplings[i][3]  = J1;
      nnbrs[i] += 4;

      neighbours[i][4]  = coord_to_site({coord[0]+1, coord[1]+1, coord[2],   coord[3]}); couplings[i][4]  = J2;
      neighbours[i][5]  = coord_to_site({coord[0]-1, coord[1]-1, coord[2],   coord[3]}); couplings[i][5]  = J2;
      neighbours[i][6]  = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   coord[3]}); couplings[i][6]  = J2;
      neighbours[i][7]  = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   coord[3]}); couplings[i][7]  = J2;
      nnbrs[i] += 4;

      neighbours[i][8]  = coord_to_site({coord[0]+2, coord[1],   coord[2],   coord[3]}); couplings[i][8]  = J3;
      neighbours[i][9]  = coord_to_site({coord[0]-2, coord[1],   coord[2],   coord[3]}); couplings[i][9]  = J3;
      neighbours[i][10] = coord_to_site({coord[0],   coord[1]+2, coord[2],   coord[3]}); couplings[i][10] = J3;
      neighbours[i][11] = coord_to_site({coord[0],   coord[1]-2, coord[2],   coord[3]}); couplings[i][11] = J3;
      nnbrs[i] += 4;
    }
  }
  else if(latticetype == "SQUAREOBC")
  {
    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)


      if (coord[0]+1 < Nx)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]+1, coord[1],   coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J1; nnbrs[i] += 1;}
      if (coord[0]-1 > -1)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]-1, coord[1],   coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J1; nnbrs[i] += 1;}
      if (coord[1]+1 < Ny)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0],   coord[1]+1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J1; nnbrs[i] += 1;}
      if (coord[1]-1 > -1)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0],   coord[1]-1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J1; nnbrs[i] += 1;}


      if ((coord[0]+1 < Nx) && (coord[1]+1 < Ny)) {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]+1, coord[1]+1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J2; nnbrs[i] += 1;}
      if ((coord[0]-1 > -1) && (coord[1]-1 > -1)) {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]-1, coord[1]-1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J2; nnbrs[i] += 1;}
      if ((coord[0]+1 < Nx) && (coord[1]-1 > -1)) {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J2; nnbrs[i] += 1;}
      if ((coord[0]-1 > -1) && (coord[1]+1 < Ny)) {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J2; nnbrs[i] += 1;}

      if (coord[0]+2 < Nx)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]+2, coord[1],   coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J3; nnbrs[i] += 1;}
      if (coord[0]-2 > -1)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0]-2, coord[1],   coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J3; nnbrs[i] += 1;}
      if (coord[1]+2 < Ny)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0],   coord[1]+2, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J3; nnbrs[i] += 1;}
      if (coord[1]-2 > -1)                        {neighbours[i][nnbrs[i]]  = coord_to_site({coord[0],   coord[1]-2, coord[2],   coord[3]}); couplings[i][nnbrs[i]]  = J3; nnbrs[i] += 1;}
    }
  }
  else if(latticetype == "TRIANGULAR")
  {
    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)
      neighbours[i][0] = coord_to_site({coord[0]+1, coord[1],   coord[2],   coord[3]}); couplings[i][0] = J1;
      neighbours[i][1] = coord_to_site({coord[0]-1, coord[1],   coord[2],   coord[3]}); couplings[i][1] = J1;
      neighbours[i][2] = coord_to_site({coord[0],   coord[1]+1, coord[2],   coord[3]}); couplings[i][2] = J1;
      neighbours[i][3] = coord_to_site({coord[0],   coord[1]-1, coord[2],   coord[3]}); couplings[i][3] = J1;
      neighbours[i][4] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   coord[3]}); couplings[i][4] = J1;
      neighbours[i][5] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   coord[3]}); couplings[i][5] = J1;
      nnbrs[i] += 6;
    }
  }
  else if(latticetype == "KAGOME")
  {
    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)
      if(coord[3] == 0)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]-1, coord[1],   coord[2],   1}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   2}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1]-1, coord[2],   2}); couplings[i][3] = J1;
      }
      if(coord[3] == 1)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1],   coord[2],   0}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   2}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   2}); couplings[i][3] = J1;
      }
      if(coord[3] == 2)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1]+1, coord[2],   0}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   1}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   1}); couplings[i][3] = J1;
      }
      nnbrs[i] += 4;
    }
  }
  else if(latticetype == "PYROCHLORE")
  {
    //Store sites as nsub*(Nx*(Nz*(nz)+ny)+nx)+s

    for(int i = 0; i < Ntot; i++)
    {
      vector<int> coord = site_to_coord(i); //coordinates of a site in (a1,a2,a3,subind)

      //Nearest neighbours:
      if(coord[3] == 0)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   1}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]-1, coord[1],   coord[2],   1}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   2}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0],   coord[1]-1, coord[2],   2}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],   3}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1],   coord[2]-1, 3}); couplings[i][5] = J1;
      }
      if(coord[3] == 1)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0]+1, coord[1],   coord[2],   0}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   2}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   2}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],   3}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, 3}); couplings[i][5] = J1;
      }
      if(coord[3] == 2)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1]+1, coord[2],   0}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   1}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   1}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],   3}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, 3}); couplings[i][5] = J1;
      }
      if(coord[3] == 3)
      {
        neighbours[i][0] = coord_to_site({coord[0],   coord[1],   coord[2],   0}); couplings[i][0] = J1;
        neighbours[i][1] = coord_to_site({coord[0],   coord[1],   coord[2]+1, 0}); couplings[i][1] = J1;
        neighbours[i][2] = coord_to_site({coord[0],   coord[1],   coord[2],   1}); couplings[i][2] = J1;
        neighbours[i][3] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, 1}); couplings[i][3] = J1;
        neighbours[i][4] = coord_to_site({coord[0],   coord[1],   coord[2],   2}); couplings[i][4] = J1;
        neighbours[i][5] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, 2}); couplings[i][5] = J1;
      }

      nnbrs[i] += 6;

      /*
      //Next nearest neighbours:
      if(coord[3] == 0)
      {
        neighbours[i][6]  = coord_to_site({coord[0],   coord[1]-1, coord[2],   1}); couplings[i][6]  = J2;
        neighbours[i][7]  = coord_to_site({coord[0],   coord[1],   coord[2]-1, 1}); couplings[i][7]  = J2;
        neighbours[i][8]  = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   1}); couplings[i][8]  = J2;
        neighbours[i][9]  = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, 1}); couplings[i][9]  = J2;
        neighbours[i][10] = coord_to_site({coord[0]-1, coord[1],   coord[2],   2}); couplings[i][10] = J2;
        neighbours[i][11] = coord_to_site({coord[0],   coord[1],   coord[2]-1, 2}); couplings[i][11] = J2;
        neighbours[i][12] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   2}); couplings[i][12] = J2;
        neighbours[i][13] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, 2}); couplings[i][13] = J2;
        neighbours[i][14] = coord_to_site({coord[0]-1, coord[1],   coord[2],   3}); couplings[i][14] = J2;
        neighbours[i][15] = coord_to_site({coord[0],   coord[1]-1, coord[2],   3}); couplings[i][15] = J2;
        neighbours[i][16] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, 3}); couplings[i][16] = J2;
        neighbours[i][17] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, 3}); couplings[i][17] = J2;
      }
      if(coord[3] == 1)
      {
        neighbours[i][6]  = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   0}); couplings[i][6]  = J2;
        neighbours[i][7]  = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, 0}); couplings[i][7]  = J2;
        neighbours[i][8]  = coord_to_site({coord[0],   coord[1]+1, coord[2],   0}); couplings[i][8]  = J2;
        neighbours[i][9]  = coord_to_site({coord[0],   coord[1],   coord[2]+1, 0}); couplings[i][9]  = J2;
        neighbours[i][10] = coord_to_site({coord[0]+1, coord[1],   coord[2],   2}); couplings[i][10] = J2;
        neighbours[i][11] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, 2}); couplings[i][11] = J2;
        neighbours[i][12] = coord_to_site({coord[0],   coord[1]-1, coord[2],   2}); couplings[i][12] = J2;
        neighbours[i][13] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, 2}); couplings[i][13] = J2;
        neighbours[i][14] = coord_to_site({coord[0]+1, coord[1],   coord[2],   3}); couplings[i][14] = J2;
        neighbours[i][15] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   3}); couplings[i][15] = J2;
        neighbours[i][16] = coord_to_site({coord[0],   coord[1],   coord[2]-1, 3}); couplings[i][16] = J2;
        neighbours[i][17] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, 3}); couplings[i][17] = J2;
      }
      if(coord[3] == 2)
      {
        neighbours[i][6]  = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   0}); couplings[i][6]  = J2;
        neighbours[i][7]  = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, 0}); couplings[i][7]  = J2;
        neighbours[i][8]  = coord_to_site({coord[0]+1, coord[1],   coord[2],   0}); couplings[i][8]  = J2;
        neighbours[i][9]  = coord_to_site({coord[0],   coord[1],   coord[2]+1, 0}); couplings[i][9]  = J2;
        neighbours[i][10] = coord_to_site({coord[0],   coord[1]+1, coord[2],   1}); couplings[i][10] = J2;
        neighbours[i][11] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, 1}); couplings[i][11] = J2;
        neighbours[i][12] = coord_to_site({coord[0]-1, coord[1],   coord[2],   1}); couplings[i][12] = J2;
        neighbours[i][13] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, 1}); couplings[i][13] = J2;
        neighbours[i][14] = coord_to_site({coord[0],   coord[1]+1, coord[2],   3}); couplings[i][14] = J2;
        neighbours[i][15] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   3}); couplings[i][15] = J2;
        neighbours[i][16] = coord_to_site({coord[0],   coord[1],   coord[2]-1, 3}); couplings[i][16] = J2;
        neighbours[i][17] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, 3}); couplings[i][17] = J2;
      }
      if(coord[3] == 3)
      {
        neighbours[i][6]  = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, 0}); couplings[i][6]  = J2;
        neighbours[i][7]  = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, 0}); couplings[i][7]  = J2;
        neighbours[i][8]  = coord_to_site({coord[0]+1, coord[1],   coord[2],   0}); couplings[i][8]  = J2;
        neighbours[i][9]  = coord_to_site({coord[0],   coord[1]+1, coord[2],   0}); couplings[i][9]  = J2;
        neighbours[i][10] = coord_to_site({coord[0],   coord[1],   coord[2]+1, 1}); couplings[i][10] = J2;
        neighbours[i][11] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, 1}); couplings[i][11] = J2;
        neighbours[i][12] = coord_to_site({coord[0]-1, coord[1],   coord[2],   1}); couplings[i][12] = J2;
        neighbours[i][13] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   1}); couplings[i][13] = J2;
        neighbours[i][14] = coord_to_site({coord[0],   coord[1],   coord[2]+1, 2}); couplings[i][14] = J2;
        neighbours[i][15] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, 2}); couplings[i][15] = J2;
        neighbours[i][16] = coord_to_site({coord[0],   coord[1]-1, coord[2],   2}); couplings[i][16] = J2;
        neighbours[i][17] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   2}); couplings[i][17] = J2;
      }

      nnbrs[i] += 12;

      //Next-next nearest neighbours:

      if(coord[3] == 0)
      {
        neighbours[i][18] = coord_to_site({coord[0]+1, coord[1],   coord[2],   coord[3]}); couplings[i][18] = J3a;
        neighbours[i][19] = coord_to_site({coord[0]-1, coord[1],   coord[2],   coord[3]}); couplings[i][19] = J3a;
        neighbours[i][20] = coord_to_site({coord[0],   coord[1]+1, coord[2],   coord[3]}); couplings[i][20] = J3a;
        neighbours[i][21] = coord_to_site({coord[0],   coord[1]-1, coord[2],   coord[3]}); couplings[i][21] = J3a;
        neighbours[i][22] = coord_to_site({coord[0],   coord[1],   coord[2]+1, coord[3]}); couplings[i][22] = J3a;
        neighbours[i][23] = coord_to_site({coord[0],   coord[1],   coord[2]-1, coord[3]}); couplings[i][23] = J3a;

        neighbours[i][24] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   coord[3]}); couplings[i][24] = J3b;
        neighbours[i][25] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   coord[3]}); couplings[i][25] = J3b;
        neighbours[i][26] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, coord[3]}); couplings[i][26] = J3b;
        neighbours[i][27] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, coord[3]}); couplings[i][27] = J3b;
        neighbours[i][28] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, coord[3]}); couplings[i][28] = J3b;
        neighbours[i][29] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, coord[3]}); couplings[i][29] = J3b;
      }
      if(coord[3] == 1)
      {
        neighbours[i][18] = coord_to_site({coord[0]+1, coord[1],   coord[2],   coord[3]}); couplings[i][18] = J3a;
        neighbours[i][19] = coord_to_site({coord[0]-1, coord[1],   coord[2],   coord[3]}); couplings[i][19] = J3a;
        neighbours[i][20] = coord_to_site({coord[0]+1, coord[1]-1, coord[2],   coord[3]}); couplings[i][20] = J3a;
        neighbours[i][21] = coord_to_site({coord[0]-1, coord[1]+1, coord[2],   coord[3]}); couplings[i][21] = J3a;
        neighbours[i][22] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, coord[3]}); couplings[i][22] = J3a;
        neighbours[i][23] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, coord[3]}); couplings[i][23] = J3a;

        neighbours[i][24] = coord_to_site({coord[0], coord[1]-1, coord[2],   coord[3]}); couplings[i][24] = J3b;
        neighbours[i][25] = coord_to_site({coord[0], coord[1]+1, coord[2],   coord[3]}); couplings[i][25] = J3b;
        neighbours[i][26] = coord_to_site({coord[0], coord[1],   coord[2]-1, coord[3]}); couplings[i][26] = J3b;
        neighbours[i][27] = coord_to_site({coord[0], coord[1],   coord[2]+1, coord[3]}); couplings[i][27] = J3b;
        neighbours[i][28] = coord_to_site({coord[0], coord[1]+1, coord[2]-1, coord[3]}); couplings[i][28] = J3b;
        neighbours[i][29] = coord_to_site({coord[0], coord[1]-1, coord[2]+1, coord[3]}); couplings[i][29] = J3b;
      }
      if(coord[3] == 2)
      {
        neighbours[i][18] = coord_to_site({coord[0],   coord[1]+1,   coord[2],   coord[3]}); couplings[i][18] = J3a;
        neighbours[i][19] = coord_to_site({coord[0],   coord[1]-1,   coord[2],   coord[3]}); couplings[i][19] = J3a;
        neighbours[i][20] = coord_to_site({coord[0]-1, coord[1]+1,   coord[2],   coord[3]}); couplings[i][20] = J3a;
        neighbours[i][21] = coord_to_site({coord[0]+1, coord[1]-1,   coord[2],   coord[3]}); couplings[i][21] = J3a;
        neighbours[i][22] = coord_to_site({coord[0],   coord[1]+1,   coord[2]-1, coord[3]}); couplings[i][22] = J3a;
        neighbours[i][23] = coord_to_site({coord[0],   coord[1]-1,   coord[2]+1, coord[3]}); couplings[i][23] = J3a;

        neighbours[i][24] = coord_to_site({coord[0]-1, coord[1], coord[2],   coord[3]}); couplings[i][24] = J3b;
        neighbours[i][25] = coord_to_site({coord[0]+1, coord[1], coord[2],   coord[3]}); couplings[i][25] = J3b;
        neighbours[i][26] = coord_to_site({coord[0]+1, coord[1], coord[2]-1, coord[3]}); couplings[i][26] = J3b;
        neighbours[i][27] = coord_to_site({coord[0]-1, coord[1], coord[2]+1, coord[3]}); couplings[i][27] = J3b;
        neighbours[i][28] = coord_to_site({coord[0],   coord[1], coord[2]-1, coord[3]}); couplings[i][28] = J3b;
        neighbours[i][29] = coord_to_site({coord[0],   coord[1], coord[2]+1, coord[3]}); couplings[i][29] = J3b;
      }
      if(coord[3] == 3)
      {
        neighbours[i][18] = coord_to_site({coord[0],   coord[1],   coord[2]+1, coord[3]}); couplings[i][18] = J3a;
        neighbours[i][19] = coord_to_site({coord[0],   coord[1],   coord[2]-1, coord[3]}); couplings[i][19] = J3a;
        neighbours[i][20] = coord_to_site({coord[0]-1, coord[1],   coord[2]+1, coord[3]}); couplings[i][20] = J3a;
        neighbours[i][21] = coord_to_site({coord[0]+1, coord[1],   coord[2]-1, coord[3]}); couplings[i][21] = J3a;
        neighbours[i][22] = coord_to_site({coord[0],   coord[1]-1, coord[2]+1, coord[3]}); couplings[i][22] = J3a;
        neighbours[i][23] = coord_to_site({coord[0],   coord[1]+1, coord[2]-1, coord[3]}); couplings[i][23] = J3a;

        neighbours[i][24] = coord_to_site({coord[0]-1, coord[1],   coord[2], coord[3]}); couplings[i][24] = J3b;
        neighbours[i][25] = coord_to_site({coord[0]+1, coord[1],   coord[2], coord[3]}); couplings[i][25] = J3b;
        neighbours[i][26] = coord_to_site({coord[0],   coord[1]+1, coord[2], coord[3]}); couplings[i][26] = J3b;
        neighbours[i][27] = coord_to_site({coord[0],   coord[1]-1, coord[2], coord[3]}); couplings[i][27] = J3b;
        neighbours[i][28] = coord_to_site({coord[0]+1, coord[1]-1, coord[2], coord[3]}); couplings[i][28] = J3b;
        neighbours[i][29] = coord_to_site({coord[0]-1, coord[1]+1, coord[2], coord[3]}); couplings[i][29] = J3b;
      }

      nnbrs[i] += 12;

      */

    }
  }
  else
  {
    cout << "ERROR: Couplings not initialized for " << latticetype << "." << endl;
    exit(1);
  }
}

vector<int> Lattice::site_to_coord(int site)
{
  //Store sites as nsub*(Nx*(Nz*(nz)+ny)+nx)+s
  vector<int> coord(4);

  coord[3] = site%subn;
  site /= subn;
  coord[0] = site%Nx;
  site /= Nx;
  coord[1] = site%Ny;
  site /= Ny;
  coord[2] = site;

  return coord;
}

int Lattice::coord_to_site(vector<int> coord)
{
  //Store sites as nsub*(Nx*(Nz*(nz)+ny)+nx)+s
  int nx = (coord[0]+Nx)%Nx;
  int ny = (coord[1]+Ny)%Ny;
  int nz = (coord[2]+Nz)%Nz;
  int s = (coord[3]+subn)%subn;
  return subn*(Nx*(Ny*(nz)+ny)+nx)+s;
}

vector<int> Lattice::GetTriangle(int site)
{
  vector<int> s0, s1, s2;
  s0 = site_to_coord(site);

  if(latticetype == "PYROCHLORE")
  {
    if(s0[3] == 0)
    {
      s1 = {s0[0]+1,s0[1]-1,s0[2],0};
      s2 = {s0[0],s0[1]-1,s0[2]+1,0};
    }
    else if(s0[3] == 1)
    {
      s1 = {s0[0],s0[1]+1,s0[2],1};
      s2 = {s0[0],s0[1],s0[2]+1,1};
    }
    else if(s0[3] == 2)
    {
      s1 = {s0[0],s0[1],s0[2]+1,2};
      s2 = {s0[0]+1,s0[1],s0[2],2};
    }
    else if(s0[3] == 3)
    {
      s1 = {s0[0]+1,s0[1],s0[2],3};
      s2 = {s0[0],s0[1]+1,s0[2],3};
    }
  }

  return {site, coord_to_site(s1), coord_to_site(s2)};
}


#endif
