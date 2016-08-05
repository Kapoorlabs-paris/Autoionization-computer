#ifndef qprop_h
#define qprop_h qprop_h
#include<iostream.h>
#include<complex.h>
#include<math.h>
#include<fstream.h>
#include<string>

class wavefunction;
class fluid;
class grid;
class hamop;
double alpha_y(double time, int me);
double vecpot_x(double time, int me);
double vecpot_y(double time, int me);
double vecpot_z(double time, int me);
double scalarpotx(double x, double y, double z, double time, int me);
double scalarpoty(double x, double y, double z, double time, int me);
double scalarpotz(double x, double y, double z, double time, int me);
double interactionpotxy(double x, double y, double z, double time, int me);
double imagpotx(long xindex, long yindex, long zindex, double time, grid g);
double imagpoty(long xindex, long yindex, long zindex, double time, grid g);
double field(double time, int me);
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &wf_one, const wavefunction &wf_two);



#endif // qprop_h




