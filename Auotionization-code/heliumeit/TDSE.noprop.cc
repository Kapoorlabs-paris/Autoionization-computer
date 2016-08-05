 
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



  FILE *file_wfdat,*file_alpha,*file_ionized,*file_wfautodat,*file_wfauto2dat;
  FILE *file_auto,*file_autoleinplustest,*file_evenautoprojection,*file_evenautoprojectionsd,*file_oddautoprojection,*file_oddautoprojectionsd,*file_autoleinatpt,*file_autosd;
  FILE *file_reading,*file_reading2,*file_reading3,*file_reading4;    // *** create some files with appropriate appendices
  char string_wfdat[]=       "/mnt/scratch/varun/res/florian/wf_smalllowdensityop.dat";
  char string_reading[]=     "/mnt/scratch/varun/res/florian/wf_tripleafterlaser.dat";
  char string_reading2[]=     "/mnt/scratch/varun/res/florian/wf_ground2500.dat";

  char string_reading3[]=     "/mnt/scratch/varun/res/florian/wf_volkovandion.dat";
  char string_reading4[]=     "/mnt/scratch/varun/res/florian/wfexcited_volkovandion.dat";

char string_evenautoprojection[]=  "/mnt/scratch/varun/res/florian/auto_smalllowdensityopevenproject_boxground.dat";



char string_oddautoprojection[]=  "/mnt/scratch/varun/res/florian/auto_smalllowdensityopoddproject_boxground.dat";


file_evenautoprojection = fopen(string_evenautoprojection,"w");



file_oddautoprojection = fopen(string_oddautoprojection,"w");


  file_reading = fopen(string_reading,"r");
file_reading2 = fopen(string_reading2,"r");



 file_reading3 = fopen(string_reading3,"r");
file_reading4 = fopen(string_reading4,"r");


   long index, rrindex, xindex, yindex,zindex,index2;
 complex<double> imagi(0.0,1.0);
  double deltx=0.5;
double delty=0.5;
double deltz=1;
  long ngpsx=2500;
long ngpsy=2500;
long ngpsxbig=2500;
long ngpsybig=2500;
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

  double real_timestep=0.05;
 
   long  freetime=0;//6800*2;
int    obs_output_every=1;
  long   wf_output_every=300000; 
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=80;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this
  double epsilon=0.00001;
   hamop hamiltonbig(gbig,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotxbig,imagpotybig,field,dftpot);

 hamop hamiltonone(gbigone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
 
 hamop hamiltontwo(gone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
   
  wavefunction wfbig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());        
 wavefunction wfinibig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
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
   long counter_iii=0;
long counter_iv=0;
long counter_v=0;
long counter_vi=0;  
  
  // initialization
  long outputofinterest=0;
 
 wfeverythingeveniny.init(gbig,4,1.0,1.375,0.0);
//  wfeverythingeveniny*=1.0/sqrt(wfeverythingeveniny.norm(gbig));

  wfeverythingoddiny.init(gbig,5,1.0,0.684,0.0);
 //  wfeverythingoddiny*=1.0/sqrt(wfeverythingoddiny.norm(gbig));
 
  
//wf.init(g,1,2.0,0.0,0.0);
      // ground

 wfread.init(g,99,0.1,0.0,0.0,file_reading2,outputofinterest);

wfground.nullify();
wfground.regrid(gbig,g,wfread);



//wfground.init(gbig,99,0.1,0.0,0.0,file_reading2,outputofinterest);
    wfread.init(g,99,0.1,0.0,0.0,file_reading,outputofinterest);
    wfbig.nullify();
 wfinibig.nullify();



    wfbig.regrid(gbig,g,wfread);
wfinibig.regrid(gbig,g,wfread);


wfchannel.init(gone,99,0.1,0.0,0.0,file_reading3,outputofinterest);

wfbigchannel.nullify();

wfbigchannel.regrid(gbigone,gone,wfchannel);




wfonebigchannel=wfbigchannel.rotate(gbigtwo,gbigone);



//wfgroundionini.initwf(gbig,0,1.0,1.375,0.0,wfeverythingeveniny,wfeverythingeveniny);

//wfgroundionini*=1.0/sqrt(wfgroundionini.norm(gbig));

//wfinibig=wfinibig-((wfgroundionini*wfinibig)*gbig.delt_x()*gbig.delt_y())*wfgroundionini;









//evenplanewave.initwfeven(gbig,0,box,1.37473,0.0,wfinibig);
//evenplanewave*=1.0/sqrt(evenplanewave.norm(gbig));




//oddplanewave.initwf(gbig,0,box,0.683633,0.0,wfinibig);

//oddplanewave*=1.0/sqrt(oddplanewave.norm(gbig));






wfsecondchannel.init(gone,99,0.1,0.0,0.0,file_reading4,outputofinterest);

 wfsecondbigchannel.nullify();

wfsecondbigchannel.regrid(gbigone,gone,wfsecondchannel);





//    wf*=1.0/sqrt(wf.norm(g));
  // wfground=wf;
  
  
     
// wfeverythingeveniny.init(gbig,4,1.0,1.375,0.0);
//  wfeverythingeveniny*=1.0/sqrt(wfeverythingeveniny.norm(gbig));

//  wfeverythingoddiny.init(gbig,5,1.0,0.684,0.0);
 //  wfeverythingoddiny*=1.0/sqrt(wfeverythingoddiny.norm(gbig));

//wfgroundion.init(gbig,6,1.0,1.375,0.0);
 
//wfgroundion*=1.0/sqrt(wfgroundion.norm(gbig));
 
 //wf=wfground;

   fclose(file_reading);
 
   fclose(file_reading2);


    
  cout << "norm wf    : " << wfinibig.norm(gbig) <<  "\n";

  


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








 complex<double> complenerg;
 
complex<double>correlationleinatpt, evencorrelationleinplustest, oddcorrelationleinplustest,evencorrelation,oddcorrelation,  evencorrelationfirst,oddcorrelationfirst, evencorrelationsecond,oddcorrelationsecond   ,evencorrelationsurff,oddcorrelationsurff,evencorrelationprojection,oddcorrelationprojection,evencorrelationprojectionini,oddcorrelationprojectionini,  evencorrelationsd,oddcorrelationsd,evencorrelationsdini,oddcorrelationsdini  ; 
    

  
  // ************* real timeprop
  
   



 

complenerg=(wfbig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge));


      


double kint=0.5;
double mint=0.5;
double deltak=0.0001;
double k;

for(k=0;k<5000;k++)
{

evenplanewave.initwfeven(gbig,0,box,kint+mint+k*deltak,0.0,wfground);

evenplanewave*=1.0/sqrt(evenplanewave.norm(gbig));







oddplanewave.initwfodd(gbig,0,box,kint+k*deltak,0.0,wfground);

oddplanewave*=1.0/sqrt(oddplanewave.norm(gbig));











oddcorrelationprojection=wfinibig.projection(gbig,oddplanewave);









evencorrelationprojection=wfinibig.projection(gbig,evenplanewave);



   



      cout << " k"  << "  "  <<  kint+k*deltak  << "  "  << " corr : " << evencorrelationprojection << "  "  << "  "  << " oddcorr : " << oddcorrelationprojection << "  " <<    endl;   





                  fprintf(file_evenautoprojection, " %.14le %.14le %.14le \n ",
            kint+mint+k*deltak,
            real(evencorrelationprojection),
                  imag(evencorrelationprojection)






                  );
                  
                  
                  
               
                  
         
                           

         






  

 

                  fprintf(file_oddautoprojection, " %.14le  %.14le %.14le \n ",
            kint+k*deltak,
            real(oddcorrelationprojection),imag(oddcorrelationprojection)






                  );

         
        


  }; 





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

double dur=ramping+constperiod+downramp;

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
        //result=-0.001;
         result=ampl*sin(ww*time)*sin(ww*time)*sin(frequ*time);
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

    
  return result;
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
  
// double ampl=0.0; // switch imaginary potential off
   double ampl=50.0; // switch imaginary potential on


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

 //double ampl=0.0; // switch imaginary potential off
   double ampl=50.0; // switch imaginary potential on


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
 // double ampl=0.0; // switch imaginary potential off
  
   double ampl=50.0; // switch imaginary potential on


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
 // double ampl=0.0; // switch imaginary potential off

   double ampl=50.0; // switch imaginary potential on


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
  
    //double ampl=0.0; // switch imaginary potential off
     double ampl=50.0; // switch imaginary potential on


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



