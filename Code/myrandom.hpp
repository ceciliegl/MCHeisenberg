#include<fstream>
#include<string>
#include "MersenneTwister.hpp"

// the class name is MTRand
#define RANDOMNUMBERGENERATOR MTRand

const string RANDOMSEED="randomseed";

unsigned long int GetRandomSeed(){
  ifstream file(RANDOMSEED.c_str());
  unsigned long int seed;
  file >> seed; return seed;
}

class RandomNumberGenerator{
public:
  RandomNumberGenerator(unsigned long int seed=21312512):
    generator(seed)
  {
  }
  double operator()(){return generator.randExc();}
  double Get(){return generator.randExc();}
  double Get0inc1(){return generator.rand();}
  double HighPrecision(){return generator.rand53();} // for better resolution
  double Ran(){return generator.randExc();}
  double RanInt(int n) {return generator.randInt(n);}
  double RanNorm(double mean, double variance){return generator.randNorm(mean, variance);}
  unsigned long int GetNewSeed(){return generator.randInt();}
 private:
  unsigned long seed;

  RANDOMNUMBERGENERATOR generator;
};

// declare the global random number generator
RandomNumberGenerator RAN(GetRandomSeed());
