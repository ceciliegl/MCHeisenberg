#ifndef VEC3
#define VEC3

#include "cvec3.hpp"

using namespace std;

class vec3

{  friend ostream& operator<<(ostream& os,vec3& c)
    { os << c.x << " " << c.y << " " << c.z; return os; }
public:
  vec3(double x0 = 0.0, double y0 = 0.0, double z0 = 0.0): x(x0), y(y0), z(z0){};
  double x;
  double y;
  double z;
};

vec3 operator+(const vec3 &u, const vec3 &v)
{
  return vec3(u.x+v.x, u.y+v.y, u.z+v.z);
}

vec3 operator-(const vec3 &u, const vec3 &v)
{
  return vec3(u.x-v.x, u.y-v.y, u.z-v.z);
}

vec3 operator*(const vec3 &u, const int &a)
{
  return vec3(a*u.x, a*u.y, a*u.z);
}

vec3 operator*(const vec3 &u, const double &a)
{
  return vec3(a*u.x, a*u.y, a*u.z);
}

vec3 operator*(const int &a, const vec3 &u)
{
  return vec3(a*u.x, a*u.y, a*u.z);
}

vec3 operator*(const double &a, const vec3 &u)
{
  return vec3(a*u.x, a*u.y, a*u.z);
}

vec3 operator/(const vec3 &u, const int &a)
{
  double ainv = 1.0/a;
  return vec3(ainv*u.x, ainv*u.y, ainv*u.z);
}

vec3 operator/(const vec3 &u, const double &a)
{
  double ainv = 1.0/a;
  return vec3(ainv*u.x, ainv*u.y, ainv*u.z);
}

cvec3 operator*(const vec3 &u, const complex<double> &a)
{
  return cvec3(a*complex<double>(u.x), a*complex<double>(u.y), a*complex<double>(u.z));
}

cvec3 operator*(const complex<double> &a, const vec3 &u)
{
  return cvec3(a*complex<double>(u.x), a*complex<double>(u.y), a*complex<double>(u.z));
}

cvec3 operator/(const vec3 &u, const complex<double> &a)
{
  complex<double> ainv = 1.0/a;
  return cvec3(ainv*complex<double>(u.x), ainv*complex<double>(u.y), ainv*complex<double>(u.z));
}

double dot(const vec3 &u, const vec3 &v)
{
  return u.x*v.x + u.y*v.y + u.z*v.z;
}

vec3 cross(const vec3 &u, const vec3 &v)
{
  return vec3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}

#endif
