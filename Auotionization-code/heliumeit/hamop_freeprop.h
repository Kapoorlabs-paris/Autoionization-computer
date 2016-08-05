#ifndef hamop_h
#define hamop_h hamop_h
#include<complex>
#include<iostream>
#include<grid.h>

using namespace std;

class wavefunction;
class fluid;

class hamop
{
 public:
  hamop(grid g, 
	
	double (*fpsx)(double, double, double, double, int),
	double (*fpsy)(double, double, double, double, int),
	double (*fpsz)(double, double, double, double, int),
	double (*fpixy)(double, double, double, double, int),
	double (*fpimx)(long, long, long, double, grid),
	double (*fpimy)(long, long, long, double, grid)
	
       );
 
  double scalarpot(double x, double y, double z, double time, int me);
  double scalarpotx(double x, double y, double z, double time, int me);
  double scalarpoty(double x, double y, double z, double time, int me);
  double scalarpotz(double x, double y, double z, double time, int me);
  double interactionpotxy(double x, double y, double z, double time, int me);
  double imagpot(long xindex, long yindex, long zindex, double time, grid g);
  double imagpotx(long xindex, long yindex, long zindex, double time, grid g);
  double imagpoty(long xindex, long yindex, long zindex, double time, grid g);
 

 private:
  double delta_x, delta_y, delta_z;
 
  double (*hamopscalarpotx)(double, double, double, double, int);
  double (*hamopscalarpoty)(double, double, double, double, int);
  double (*hamopscalarpotz)(double, double, double, double, int);
  double (*hamopinteractionpotxy)(double, double, double, double, int);
  double (*hamopimagpotx)(long, long, long, double, grid);
  double (*hamopimagpoty)(long, long, long, double, grid);
 


};



#endif // hamop_h




