#ifndef CVEC3
#define CVEC3

using namespace std;

class cvec3

{  friend ostream& operator<<(ostream& os,cvec3& c)
    { os << c.x << " " << c.y << " " << c.z; return os; }
public:
  cvec3(complex<double> x0 = complex<double>(0.0,0.0), complex<double> y0 = complex<double>(0.0,0.0), complex<double> z0 = complex<double>(0.0,0.0)): x(x0), y(y0), z(z0){};
  complex<double> x;
  complex<double> y;
  complex<double> z;
};

cvec3 operator+(const cvec3 &u, const cvec3 &v)
{
  return cvec3(u.x+v.x, u.y+v.y, u.z+v.z);
}

cvec3 operator-(const cvec3 &u, const cvec3 &v)
{
  return cvec3(u.x-v.x, u.y-v.y, u.z-v.z);
}

cvec3 operator*(const cvec3 &u, const complex<double> &a)
{
  return cvec3(a*u.x, a*u.y, a*u.z);
}

cvec3 operator*(const complex<double> &a, const cvec3 &u)
{
  return cvec3(a*u.x, a*u.y, a*u.z);
}

cvec3 operator/(const cvec3 &u, const complex<double> &a)
{
  complex<double> ainv = complex<double>(1.0)/a;
  return cvec3(ainv*u.x, ainv*u.y, ainv*u.z);
}

complex<double> dot(const cvec3 &u, const cvec3 &v)
{
  return conj(u.x)*v.x + conj(u.y)*v.y + conj(u.z)*v.z;
}

/*cvec3 cross(const cvec3 &u, const cvec3 &v)
{
  return cvec3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
}*/

#endif
