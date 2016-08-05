 
#include<TDSE.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

#define pi 3.1415926
double alphahat;
double frequref;
int main(int argc, char **argv)
{
  int me=0;



  FILE *file_wfdat,*file_obser,*file_reading;
     // *** create some files with appropriate appendices
  char string_wfdat[]=       "res/wf_singlet_He+_10_28.dat";  

char string_obser[]=       "res/observ_singlet_He+_10_28.dat"; 
char string_reading[]=       "res/wf_singlet_He+_10_27.dat"; /





  file_wfdat = fopen(string_wfdat,"w");
  file_obser = fopen(string_obser,"w");

file_reading = fopen(string_reading,"r");


   long index, rrindex, xindex, yindex,zindex,index2;
 complex<double> imagi(0.0,1.0);
  double deltx=0.5;
double delty=0.5;
double deltz=1;
  long ngpsx=1500;
long ngpsy=1500;
long ngpsxbig=1500;
long ngpsybig=1500;
  long ngpsz=1;



  // *** declare grid ***
  grid gbig;
  gbig.set_dim(16);            // propagation mode is 3 for xy cartesian
  gbig.set_ngps(ngpsxbig,ngpsybig,ngpsz);    // N_x, N_y, 1 was 100,100,1
  gbig.set_delt(deltx,delty,deltz);  // delta_x, delta_y, deltz was 0.1,0.1,0.1
  gbig.set_offs(ngpsxbig/2,ngpsybig/2,0);    // origin  (usually at N_x/2, N_y/2,N_z/2)






  // *** declare grid ***
  grid g;
  g.set_dim(16);            // propagation mode is 3 for xy cartesian
  g.set_ngps(ngpsx,ngpsy,ngpsz);    // N_x, N_y, 1 was 100,100,1
  g.set_delt(deltx,delty,deltz);  // delta_x, delta_y, deltz was 0.1,0.1,0.1
  g.set_offs(ngpsx/2,ngpsy/2,0);    // origin  (usually at N_x/2, N_y/2,N_z/2)
 
   // *** declare grid for 1D ***
  grid gone;
  gone.set_dim(15);
  gone.set_ngps(ngpsx,1,1);
  gone.set_delt(deltx,1,1);
  gone.set_offs(ngpsx/2,0,0);
  
  


   // *** declare grid for 1D ***
  grid gbigone;
  gbigone.set_dim(15);
  gbigone.set_ngps(ngpsxbig,1,1);
  gbigone.set_delt(deltx,1,1);
  gbigone.set_offs(ngpsxbig/2,0,0);




   
  
 
   // *** declare grid for 1D ***
  grid gbigtwo;
  gbigtwo.set_dim(15);
  gbigtwo.set_ngps(1,ngpsybig,1);
  gbigtwo.set_delt(1,delty,1);
  gbigtwo.set_offs(0,ngpsybig/2,0);



  
  
  // *** declare rest of variables ***
double imag_timestep=0.5;
  double real_timestep=0.25;
 long   no_of_imag_timesteps=0;
   long   no_of_real_timesteps=20000;
int    obs_output_every=1;
  int    wf_output_every=1000;
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=60;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this
  double epsilon=0.00001;
   hamop hamiltonbig(gbig,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotxbig,imagpotybig,field,dftpot);

 hamop hamiltonone(gbigone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
 
 hamop hamiltontwo(gone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
   
  wavefunction wfbig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());        
 wavefunction wfinibig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());


 wavefunction wfsecondbig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());

 
          wavefunction wfprobig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());

  wavefunction wfground(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());


         

  wavefunction wfexcitedheliumplus(g.ngps_x()*g.ngps_y()*g.ngps_z());

 wavefunction wfheliumplus(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 

 wavefunction wfchannel(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 wavefunction wfsecondchannel(gone.ngps_x()*gone.ngps_y()*gone.ngps_z());
 wavefunction wfbigchannel(gbigone.ngps_x()*gbigone.ngps_y()*gbigone.ngps_z());
 wavefunction wfsecondbigchannel(gbigone.ngps_x()*gbigone.ngps_y()*gbigone.ngps_z());




               

wavefunction wfread(g.ngps_x()*g.ngps_y()*g.ngps_z());

wavefunction wfonebigchannel(gbigtwo.ngps_x()*gbigtwo.ngps_y()*gbigtwo.ngps_z());
wavefunction wfreadexcited(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfreadexcited2(g.ngps_x()*g.ngps_y()*g.ngps_z());

wavefunction wfreadexcited3(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfreadexcited4(g.ngps_x()*g.ngps_y()*g.ngps_z());


 wavefunction planewave(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());


   wavefunction wfini(g.ngps_x()*g.ngps_y()*g.ngps_z());
  
 wavefunction wfini2(g.ngps_x()*g.ngps_y()*g.ngps_z());


wavefunction wfeverythingeveniny(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
  wavefunction wfeverythingoddiny(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());

wavefunction wfeverythingeven(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction wfeverythingodd(g.ngps_x()*g.ngps_y()*g.ngps_z());
 
    wavefunction correlation(g.ngps_x()*g.ngps_y()*g.ngps_z()) ;
  wavefunction correlationex(g.ngps_x()*g.ngps_y()*g.ngps_z()) ;
   
wavefunction wfgroundion(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
 
 wavefunction wfgroundionini(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());




wavefunction evenplanewave(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
  wavefunction oddplanewave(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());





//double posofinteresty=2;
// double posofinterestx=3;
    long posofinterest=3200;
  double posofinterestKHx,posofinterestKHy;


  long gpKHx,gplabx,gpKHy,gplaby;

  

  complex<double> timestep;

complex<double> timestep2;

  double time=0.0;
double time2=0.0;  
  long counter_i=0;
  long counter_ii=0;
   
  double frequnull= 1.356;//0.533*0.5;//0.438;  //0.219;
  // initialization
  long outputofinterest=0;
 double energyproject;
 



wfread.init(g,99,0.1,0.0,0.0,file_reading,outputofinterest);
wfbig.nullify();
wfbig.regrid(gbig,g,wfread);

//wfbig.init(gbig,1,0.05,1.38,0.0);

 wfbig*=1.0/sqrt(wfbig.norm(gbig));
 








wfprobig.nullify();

energyproject=-2.2385+frequnull;//-1.8161+frequnull;


//wfprobig*=1.0/sqrt(wfprobig.norm(gbig)); 






complex<double> complenerg;
double norm;
    
  cout << "norm wf    : " << wfprobig.norm(gbig) <<  "\n";

  


 wavefunction staticpot_xbig(gbig.ngps_x());
  staticpot_xbig.calculate_fixed_potential_array_x(gbig,hamiltonbig,0.0,me);

  wavefunction staticpot_ybig(gbig.ngps_y());
  staticpot_ybig.calculate_fixed_potential_array_y(gbig,hamiltonbig,0.0,me);



  wavefunction staticpotbig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
  staticpotbig.calculate_fixed_potential_array(gbig,hamiltonbig,0.0,me);




  wavefunction staticpotdftbig(gbigone.ngps_x()*gbigone.ngps_y()*gbigone.ngps_z());
  staticpotdftbig.calculate_fixed_potential_array(gbigone,hamiltonone,0.0,me);



     wavefunction staticpot_xybig(gbig.ngps_z()*gbig.ngps_y()*gbig.ngps_x());
  staticpot_xybig.calculate_fixed_potential_array_xy(gbig,hamiltonbig,0.0,me);




 
 
 long ts;
  long no_of_timesteps=no_of_imag_timesteps;
 
for (ts=0; ts<no_of_timesteps; ts++)
    {   

      cout << "Imag: " << ts <<  " "  << " energ : " << real(complenerg) <<  endl;

      counter_i++;
      counter_ii++;
      timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep); 

 // wfsecondbig.propagate(timestep,0.0,gbig,hamiltonbig,me,vecpotflag,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge);  

//wfsecondbig=wfsecondbig-((wfbig*wfsecondbig)*g.delt_x()*g.delt_y())*wfbig;


// wfsecondbig*=1.0/sqrt(wfsecondbig.norm(gbig));


complenerg=(wfsecondbig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge));
   

};








    
      counter_i=0;
      counter_ii=0;
 
     
  

   


 // ************* real timeprop
  
           
  
  
  no_of_timesteps=no_of_real_timesteps;
  
  
 timestep=complex<double>(real_timestep,0.0);

complex<double> delta=0.0;

for (ts=0; ts<no_of_timesteps; ts++)
    {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
        cout << "Real: " << ts <<  " "  << " Energie: " << real(complenerg) <<  " "  <<"\n";
    



    wfbig.propagate(timestep,time,gbig,hamiltonbig,me,vecpotflag,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge);  
   
  	
//wfbig*=1.0/sqrt(wfbig.norm(gbig)); 
  
  
  
    
wfprobig=wfprobig+0.5*(1.0-cos(2.0*pi*time/(no_of_timesteps)))*real_timestep*exp(imagi*time*energyproject)*wfbig;

//wfprobig*=0.5*(1.0-cos(2.0*pi*time/no_of_timesteps));

norm=(wfprobig.norm(gbig));


	     
complenerg=(wfprobig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge))/norm;	  
	
//wfprobig*=1.0/sqrt(wfprobig.norm(gbig)); 
 
    
 fprintf(file_obser,"%.14le %.14le %.14le\n ",
	    (time+real(timestep)),
	    real(complenerg),wfprobig.norm(gbig));
   

 

 }


wfprobig*=1.0/sqrt((wfprobig.norm(gbig)));
 
complenerg=(wfprobig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge));
norm=(wfprobig.norm(gbig));
   
 
wfprobig.dump_to_file(gbig,file_wfdat,dumpingstepwidth);



    fclose(file_wfdat);
cout<<" "  << " energ : " << real(complenerg) << "  " << endl;

 cout << me << ": Hasta la vista, ... " << endl;

};


double vvvvecpot_x(double time, int me)
{
  double result=0.0;

  double frequ=frequref;
  double n=100.0;
  double ampl=alphahat*frequ;
  double phi=0.0;
  double dur=n*2.0*M_PI/frequ;
double gap=2.0*M_PI/frequ;
 double ww=0.5*frequ/n;
  

 
   if ((time>0.0)  )
    {
	//result=-0.001;
	// result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time-phi);
 };



  
      

 




  return result;

}  


double vecpot_x(double time, int me)
{

  double result;

  double frequ= frequref;
  double frequ2= 0.0;//0.038;
  double ramping=10.0*2*M_PI/frequ;
  double downramp=10.0*2*M_PI/frequ;
 
  double constperiod=60*2*M_PI/frequ;//392.0*2*M_PI/frequ;
double delay=40.0*2*M_PI/frequ;  

double dur=200;//ramping+constperiod+downramp;

  double ampl = alphahat*frequ;
 double n=80.0;
double n2=3.0;
  double ww=0.5*frequ/n;
double ww2=0.5*frequ2/n2;
double dur2=3.0*2*M_PI/frequ2;

   if (time>0)  
    {

//result=-0.001;
if(time<dur)
{
        result=-0.00;
        // result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time);
} 

if((time>dur))//&& (time<dur+delay))
{
result=0.0;

}
};
/*
if((time>dur+delay)&&(time<dur+delay+dur2))
{
result=ampl*sin(ww2*time)*sin(ww2*time)*sin(frequ2*time);
}

};

/*
 if (time>0)
    {
if ((time<=ramping+constperiod+downramp) )
{
	
	result=ampl*exp(-(2.0*0.693147)*(time-dur)*(time-dur)/(dur*dur/16.0))*cos(frequ*(time));
	
}

if ((time>ramping+constperiod+downramp) )
{
	
	//result=ampl2*time*exp(-(2.0*0.693147)*time*time/(dur*dur))*sin(frequ2*time);
result=0.0;	
}
};
	
	

  
 
 /*

 
  //double n= dur*frequ/(2*M_PI);

//double gap=2.0*2*M_PI/frequ;

if (time>0)
{
	result=-0.001;

}



//double ramping2=2.0*2*M_PI/frequ2;
 // double downramp2=2.0*2*M_PI/frequ2;
 
  //double constperiod2=104.0*2*M_PI/frequ2;
//double ampl2 = 5.0*alphahat*frequ2;
 
 
  if (time>0)
    {

      if (time<ramping){
	result=-ampl/ramping*time*cos(frequ*time) + ampl/(frequ*ramping)*sin(frequ*time);
	}
     

  
 
      if ( (time>=ramping)&&(time<ramping+constperiod))
      {
	result=-ampl*cos(frequ*time);//-ampl*cos(frequ2*time);
      }
     
      if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp)) 
	{
		  result=ampl/downramp*(time-downramp-constperiod-ramping)*cos(frequ*time)-ampl/(frequ*downramp)*sin(frequ*time);
	}
      
      if ( (time >= ramping+constperiod+downramp )) //&& (time<ramping+constperiod+downramp +gap)) 
      {
      result=0.0;
 }
  };
 /*
 
   

      if ((time>ramping+constperiod+downramp+gap) &&(time<ramping2+ ramping+constperiod+downramp+gap )){
	result=-ampl2/(ramping2)*time*cos(frequ2*time) + ampl2/(frequ2*(ramping2))*sin(frequ2*time);
	}
      
      if ( (time>=ramping2+ramping+constperiod+downramp+gap)&&(time<ramping2+constperiod2+ramping+constperiod+downramp+gap))
      {
	result=-ampl2*cos(frequ2*time);//-ampl*cos(frequ2*time);
      }
     
      if ( (time >= ramping2+constperiod2+ramping+constperiod+downramp+gap) && (time < ramping2+constperiod2+downramp2+ramping+constperiod+downramp+gap)) 
	{
		  result=ampl2/downramp2*(time-downramp2-constperiod2-ramping2+ramping+constperiod+downramp+gap)*cos(frequ*time)-ampl2/(frequ2*downramp2)*sin(frequ2*time);
	}
      
      if ( (time >= ramping2+constperiod2+downramp2+ramping+constperiod+downramp+gap ) ) 
      {
      result=0.0;
 }
 */

    
  return 0.0*result;
}  



double vvvvvvecpot_x(double time, int me)
{

  double result=0.0;
  return result;
}  


double alpha_y(double time, int me)
{
  double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

   // <-------- put same in vecpot_y !!!
  vecpotampl = alphahat*frequref;


  result=-vecpotampl/frequref*cos(frequref*time);

  return result;

}  


double alpha_x(double time, int me)
{
 double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

  // <-------- put same in vecpot_y !!!
  vecpotampl = alphahat*frequref;


  result=-vecpotampl/frequref*cos(frequref*time);

  return result;

 

}  

double alpha_z(double time, int me)
{
 

  return alpha_y(time,me);

}  


double vecpot_y(double time, int me)
{
 return vecpot_x(time,me);
}  


double vecpot_z(double time, int me)

{

 
 
    

  return 0.0;
 
 // return vecpot_x(time,me);
}  








double scalarpotx(double x, double y, double z, double time, int me)
{
  double eps=1.0;
 double result;
 
 double slope=0.00522;
 double Static_field=0.0141;
 double field_ramp=27;
 double field_const=270;
 result= -2.0/sqrt(x*x+eps);   //x*Static_field-2.0/sqrt(x*x+eps);
 /*
 if (time>0)
    {
		result=x*Static_field-2.0/sqrt(x*x+eps);

if (time<field_ramp){
	result=x*slope*time-2.0/sqrt(x*x+eps);
	}
	if ((time>=field_ramp) && (time<field_ramp+field_const))
	{
		result=x*Static_field-2.0/sqrt(x*x+eps);
	}
	if(time>=field_ramp+field_const) {
		result= -2.0/sqrt(x*x+eps);
	}
	
   //return -2.0/sqrt(x*x+eps);
        //     return 0;
	};
	 double eps2=0.5034;
     return -3.0/sqrt(x*x+eps2*eps2);
    */
	return result;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double eps=1.0;
    
double result;
 double slope=0.00522;
 double Static_field=0.00141;
 double field_ramp=27;
 double field_const=270;
 
 result= -2.0/sqrt(y*y+eps);   //y*Static_field-2.0/sqrt(y*y+eps);
 /*
 if (time>0)
    {
		

if (time<field_ramp){
	result=y*slope*time-2.0/sqrt(y*y+eps);
	}
	if ((time>=field_ramp) && (time<field_ramp+field_const))
	{
		result=y*Static_field-2.0/sqrt(y*y+eps);
	}
	if(time>=field_ramp+field_const) {
		result= -2.0/sqrt(y*y+eps);
	}
	
   //return -2.0/sqrt(x*x+eps);
        //     return 0;
	};
 //return -2.0/sqrt(y*y+eps);
      //         return 0;
      */
      return result;
// double eps2=0.5034;
  //   return -3.0/sqrt(y*y+eps2*eps2);
    
}

double scalarpotz(double x, double y, double z, double time, int me)
{
double eps=1.0;

 return 0.0;
 //return -3.0/sqrt(z*z+eps);
}

double interactionpotxy(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  
  return 1.0/sqrt((x-y)*(x-y)+eps);
   
   double eps2=0.5034;
  
  // return 1.0/sqrt((x-y)*(x-y)+eps2*eps2);
   // return 0.0;
}

double interactionpotyz(double x, double y, double z, double time, int me)
{
  double eps=1;
   return 0.0;  
 //return 1.0/sqrt((z-y)*(z-y)+eps);
 
}
double interactionpotxz(double x, double y, double z, double time, int me)
{
  double eps=1;
    return 0.0;  
 //return 1.0/sqrt((x-z)*(x-z)+eps);
 
}
  
double field(double time, int me)
{
  double result=0.0;


  return result;


}  







double imagpotx(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  
 double ampl=0.0; // switch imaginary potential off
 //  double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpotxbig(long xindex, long yindex, long zindex, double time, grid gbig)
{
  double x,y,z;

 double ampl=0.0; // switch imaginary potential off
 //double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*gbig.ngps_x())/(0.5*gbig.ngps_x())
        *((double) xindex + 0.5  - 0.5*gbig.ngps_x())/(0.5*gbig.ngps_x());

      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };

}

double imagpoty(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
 double ampl=0.0; // switch imaginary potential off
  
 //  double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 


double imagpotybig(long xindex, long yindex, long zindex, double time, grid gbig)
{
  double y;
  double ampl=0.0; // switch imaginary potential off

  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*gbig.ngps_y())/(0.5*gbig.ngps_y())
            *((double) yindex + 0.5  - 0.5*gbig.ngps_y())/(0.5*gbig.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };

}



double imagpotz(long xindex, long yindex, long zindex, double time, grid g)
{
  double z;
  // double ampl=0.0; // switch imaginary potential off
   double ampl=0.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	    *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
       return ampl*z*z*z*z*z*z*z*z;
    }
  else
    {
       return 0.0;
    };
      
} 






double scalarpotxtwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;

  return -2.0/sqrt(x*x+eps);
         //
          //  double eps=0.5034;
     //return -3.0/sqrt(x*x+eps2*eps2);
    
           
}

double scalarpotytwo(double x, double y, double z, double time, int me)
{
   double result=0.0;
  return result;

}

double scalarpotztwo(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double interactionpotxytwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  double interactionpotyztwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  double interactionpotxztwo(double x, double y, double z, double time, int me)
{
  double eps=1.0;
  //  return 1.0/sqrt((x-y)*(x-y)+eps);
  return 0.0; // attention --- no interaction !!!!!!!!!!!!!!
}

  

double imagpotxtwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double x,y,z;
  
    double ampl=0.0; // switch imaginary potential off
     //double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());
  
      return ampl*x*x*x*x*x*x*x*x;
    }
  else
    {
       return 0.0;
    };
      
}  


double imagpotytwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double y;
    double ampl=0.0; // switch imaginary potential off
  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       y=((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y())
	    *((double) yindex + 0.5  - 0.5*g.ngps_y())/(0.5*g.ngps_y());
       return ampl*y*y*y*y*y*y*y*y;
    }
  else
    {
       return 0.0;
    };
      
} 

double imagpotztwo(long xindex, long yindex, long zindex, double time, grid g)
{
  double z;
    double ampl=0.0; // switch imaginary potential off
  // double ampl=50.0; // switch imaginary potential on


  if (ampl>1.0)
    {
       z=((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z())
	    *((double) zindex + 0.5  - 0.5*g.ngps_z())/(0.5*g.ngps_z());
       return ampl*z*z*z*z*z*z*z*z;
    }
  else
    {
       return 0.0;
    };
      
} 

double dfthartree(grid gone, double x, double y, double z, double time, int me, 
	      const wavefunction & wf, double eps, long box)
{
  double result=0.0;
  return result;
}



// dft is not used
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &v_null, const wavefunction &v_eins)
{
  double result;
  result=0.0;
  return result;
}



