


#include<TDSE.h>
#include<stdio.h>
#include<wavefunction.h>
#include<fluid.h>
#include<grid.h>
#include<hamop.h>

double alphahat;
double frequofinterest;
double Mp = 1.0;//836.0;

#define pi 3.1415926
double scalingfactor=sqrt(0.5*Mp);
int main(int argc, char **argv)
{
  int me=0;

  
  FILE *file_nuclpot,*file_elecpot;

 



  FILE  *file_wfdat, *file_KHwf,*file_wf2dat,*file_wf3dat,*file_wf4dat,*file_wf5dat,*file_wf6dat,*file_wf7dat,*file_wf8dat ;
  FILE *file_obser, *file_obser_imag;
  FILE *file_wfatgpoi, *file_wfatgpoiel, *file_dipole_y;
  FILE *file_reading, *file_floquwf,*file_reading2;
  // *** create some files with appropriate appendices
 
  // char string_nuclpot[]=     "/home/theo/kapoor/d34/res/florian/nuclpolthelium.dat";
  // char string_elecpot[]=     "/home/theo/kapoor/d34/res/florian/elecpothelium.dat";
 char string_wfdat[]=      "/home/theo/kapoor/d34/res/florian/wfhelium.dat";
 /* char string_wf2dat[]=      "/home/theo/kapoor/d34/res/florian/wf2localalp5.dat";
 char string_wf3dat[]=      "/home/theo/kapoor/d34/res/florian/wf3localalp5.dat";
 char string_wf4dat[]=      "/home/theo/kapoor/d34/res/florian/wf4localalp5.dat";
 char string_wf5dat[]=      "/home/theo/kapoor/d34/res/florian/wf5localalp5.dat";
 char string_wf6dat[]=      "/home/theo/kapoor/d34/res/florian/wf6localalp5.dat";
 char string_wf7dat[]=      "/home/theo/kapoor/d34/res/florian/wf7localalp5.dat";
 char string_wf8dat[]=      "/home/theo/kapoor/d34/res/florian/wf8localalp5.dat";
 */  
char string_KHwf[]=       "/home/theo/kapoor/d34/res/florian/KHwfhelium.dat";
  char string_floquwf[]=    "/home/theo/kapoor/d34/res/florian/floquwfhelium.dat";
  char string_obser[]=      "/home/theo/kapoor/d34/res/florian/observhelium.dat";
  char string_obser_imag[]= "/home/theo/kapoor/d34/res/florian/observimaghelium.dat";
  char string_wfatgpoi[]=   "/home/theo/kapoor/d34/res/florian/wfatgpoihelium.dat";
    char string_wfatgpoiel[]=   "/home/theo/kapoor/d34/res/florian/wfatgpoielhelium.dat";
  char string_dipole_y[]=   "/home/theo/kapoor/d34/res/florian/dipole_yhelium.dat";
  char string_reading[]=    "/home/theo/kapoor/d34/res/florian/ini/Helium200.dat";
char string_reading2[]=    "/home/theo/kapoor/d34/res/florian/potential400alp1.4.dat";
  
//file_nuclpot = fopen(string_nuclpot,"w");
//file_elecpot = fopen(string_elecpot,"w");
 file_wfatgpoiel = fopen(string_wfatgpoiel,"w");
  file_wfdat = fopen(string_wfdat,"w");
  //  file_wf2dat = fopen(string_wf2dat,"w");
  /*  
file_wf3dat = fopen(string_wf3dat,"w");
  
file_wf4dat = fopen(string_wf4dat,"w");
  
file_wf5dat = fopen(string_wf5dat,"w");
  
file_wf6dat = fopen(string_wf6dat,"w");
  
file_wf7dat = fopen(string_wf7dat,"w");
  
file_wf8dat = fopen(string_wf8dat,"w");
  
  */
 file_KHwf = fopen(string_KHwf,"w");
  file_floquwf = fopen(string_floquwf,"w");
  file_obser = fopen(string_obser,"w");
  file_obser_imag = fopen(string_obser_imag,"w");
  file_wfatgpoi = fopen(string_wfatgpoi,"w");
  file_dipole_y = fopen(string_dipole_y,"w");
  //  file_reading = fopen(string_reading,"r");
  // file_reading2 = fopen(string_reading2,"r");
 // *** declare grid ***

  double deltx=0.2;//0.05*sqrt(0.5*Mp);
double delty=0.2;
 long index, indexKH, xindexKH, rrindex, xindex, yindex,yindexKH,indexel,xindexel,yindexel;
 

 long ngpsx,ngpsy,alp;
 // loop:

 ngpsx=200;
 ngpsy=200;

  grid g;
  g.set_dim(16);            // propagation mode is 16 for xy cartesian
  g.set_ngps(ngpsx,ngpsy,1);    // N_x, N_y, 1 was 100,100,1
  g.set_delt(0.2,0.2,1.0);  // delta_x, delta_y, 1.0 was 0.5,0.5,1.0
  g.set_offs(ngpsx/2,ngpsy/2,0);    // origin  (usually at N_x/2, N_y/2)

  // *** declare smaller grid for reading***
  grid g_small;
  g_small.set_dim(16);
  g_small.set_ngps(ngpsx,ngpsy,1);
  g_small.set_delt(0.2,0.2,1.0);
 // g_small.set_delt(0.05*sqrt(0.5*Mp),0.2,1.0);
  g_small.set_offs(ngpsx/2,ngpsy/2,0);
  /* 
  // *** declare 1D grid for output of effective nuclear potential***
  grid g_only_x;
  g_only_x.set_dim(16);
  g_only_x.set_ngps(g.ngps_x(),1,1);
  g_only_x.set_delt(g.delt_x(),1.0,1.0);
  g_only_x.set_offs(0,0,0);
  
 // *** declare 1D grid for output of effective electric potential***
  grid g_only_y;
  g_only_y.set_dim(16);
  g_only_y.set_ngps(1,g.ngps_y(),1);
  g_only_y.set_delt(1.0,g.delt_y(),1.0);
  g_only_y.set_offs(0,0,0);
  
  */
 
   
  // *** declare rest of variables ***
  double imag_timestep=0.05;
  double real_timestep=0.05;
  long   no_of_imag_timesteps=100;
  long   no_of_real_timesteps=0;
  int    obs_output_every=1;
  long   wf_output_every=100;
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=50;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this

  
   hamop hamilton(g,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
  
   //wavefunction nuclpot(g_only_x.ngps_x());
   //wavefunction elecpot(g_only_y.ngps_y());
 
  wavefunction wfexcited(g.ngps_x()*g.ngps_y()*g.ngps_z());
  // wavefunction wfel(g_only_x.ngps_x()*g_only_y.ngps_y()*g_only_y.ngps_z()); 
  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
 wavefunction wfground(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfeverythingeven(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction wfeverythingodd(g.ngps_x()*g.ngps_y()*g.ngps_z());
  //wavefunction wfeverythingevenel(g_only_y.ngps_x()*g_only_y.ngps_y()*g_only_y.ngps_z());
  //wavefunction wfeverythingoddel(g_only_y.ngps_x()*g_only_y.ngps_y()*g_only_y.ngps_z());
 wavefunction wfini(g.ngps_x()*g.ngps_y()*g.ngps_z());


wavefunction wfread(g_small.ngps_x()*g_small.ngps_y()*g_small.ngps_z());



double posofinteresty=2;
 double posofinterestx=3;//*sqrt(0.5*Mp);
  //  double posofinterest=2;
  double posofinterestKHx,posofinterestKHy;

  long gpKHx,gplabx,gpKHy,gplaby;

  
long noofalphas=1;
  double deltaalpha=0.1;
  double alphanull=0.1;
  
 long alphacounter;
  
  long nooffrequs=1;
  double deltafrequ=0.002;
  double frequnull=12;
  long frequcounter;


  long nooftargetenergs=1;
  double deltatargetenergs=0.02;
  double targetenergnull=-0.45;
  long targetenergcounter;
  double targetenergy;


  // *** grid for floquet wavefunction***
  grid g_floqu;
  g_floqu.set_dim(3);
  g_floqu.set_ngps(ngpsx,ngpsy,nooftargetenergs); // <== z runs over the energies of interest 
  g_floqu.set_delt(deltx,delty,1.0);
  g_floqu.set_offs(ngpsx/2,ngpsy/2,0);

  wavefunction floquwf(g_floqu.ngps_x()*g_floqu.ngps_y()*g_floqu.ngps_z());
  
  // *** grid for floquet wavefunction in KH frame ***
  grid g_KH;
  g_KH.set_dim(3);
  g_KH.set_ngps(ngpsx,ngpsy,nooftargetenergs);
  g_KH.set_delt(deltx,delty,1.0);
  g_KH.set_offs(ngpsx/2,ngpsy/(2),0);

  wavefunction KHwf(g_KH.ngps_x()*g_KH.ngps_y()*g_KH.ngps_z());


  complex<double> timestep;
  double time=0.0;
 complex<double> imagi(0.0,1.0);
    long counter_i=0;
    long counter_ii=0;

 


  //     hamop hamiltonread(g_small.ngps_x()*g_small.ngps_y()*g_small.ngps_z());
  //     hamop hamilton(g.ngps_x()*g.ngps_y()*g.ngps_z());


    // initialization
    long outputofinterest=0;
   
     wf.init(g,1,7.0,0.0,0.0);
    
    //   wfread.init(g_small,99,0.1,0.0,0.0,file_reading,outputofinterest);
     wf*=1.0/sqrt(wf.norm(g));
    // wf.nullify();
    // wf.regrid(g,g_small,wfread);
      // wfexcited.init(g,2,7.0,0.0,0.0);
       


 wfeverythingeven.init(g,4,7.0,0.0,0.0);
  wfeverythingeven*=1.0/sqrt(wfeverythingeven.norm(g));

  wfeverythingodd.init(g,5,7.0,0.0,0.0);
   wfeverythingodd*=1.0/sqrt(wfeverythingodd.norm(g));




	    
   // fclose(file_reading);
 
  	
   


	 wavefunction staticpot_xy(g.ngps_x()*g.ngps_y()*g.ngps_z());
	 staticpot_xy.calculate_fixed_potential_array_xy(g,hamilton,0.0,me);//calculate_fixed_potential_array_xy(g,hamilton,0.0,me);
//initmine(g,99,0.1,0.0,0.0,file_reading2,outputofinterest);//calculate_fixed_potential_array_xy(g,hamilton,0.0,me);
  
	 // staticpot_xy.regrid(g,g_small,staticpot_xy);
	 //	  fclose(file_reading2);
     
 wavefunction staticpot_y(g.ngps_y());
        staticpot_y.calculate_fixed_potential_array_y(g,hamilton,0.0,me);
  cout<<"calculated staticpot_x_y"<<endl;


    wavefunction staticpot_x(g.ngps_x());
       
   staticpot_x.calculate_fixed_potential_array_x(g,hamilton,0.0,me);
   // staticpot_x.regrid(g,g_small,staticpot_x);
   // fclose(file_reading3);
 

	  //  wavefunction staticpot_xy_file(g.ngps_x()*g.ngps_y()*g.ngps_z());
	  //  staticpot_xy_file.calculate_fixed_potential_array_xy_file(g,hamilton,0.0,me);
	
    complex<double> complenerg;
    complex<double> groundstatepop, excitedstatepop1, excitedstatepop2, excitedstatepop3,excitedstatepop4,excitedstatepop5,excitedstatepop6,excitedstatepop7,excitedstatepop8,excitedstatepop9,excitedstatepop10,excitedstatepop11,excitedstatepop12,excitedstatepop13,excitedstatepop14,excitedstatepop15,excitedstatepop16  ;
 complex<double> KHnorm;
  complex<double> complenergex;
  complex<double> correlation;
  complex<double> correlationex;
     complex<double> correlationel;
  complex<double> correlationexel;
   
    // ************* imag timeprop
    long ts;
    long no_of_timesteps=no_of_imag_timesteps;

    for (ts=0; ts<no_of_timesteps; ts++)
      {   

	cout << "Imag: " << ts  << "   "  <<   "norm wf    : " << wf.norm(g) << "   " << "energy   :" <<  "   " << real(complenerg) << "   "  << endl;

	counter_i++;
	counter_ii++;
	timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep);  


	// and now the actual propagation
		wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);
//	wf=wf-((wfground*wf)*g.delt_y()*g.delt_x())*wfground;
	wf*=1.0/sqrt(wf.norm(g));
 
       

       
   

	if (counter_ii==obs_output_every)
      	{
 	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge); 
 
	  fprintf(file_obser_imag,"%i %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
		  ts,real(complenerg),imag(complenerg),wf.norm(g),wf.molecule_ionized(g,box),wf.sing_ionized(g,box),wf.doub_ionized(g,box),wf.expect_x(g),wf.expect_y(g));
	  counter_ii=0;
	};

    };


    fclose(file_obser_imag);
	
        wf.dump_to_file(g,file_wfdat,dumpingstepwidth);




      

    //    wf.nuclenergysurface(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge,nuclpot); 
    // nuclpot.dump_to_file(g_only_x,file_nuclpot,dumpingstepwidth);

  

   

   
  wfini=wf;

 

   
  for (alphacounter=0; alphacounter<noofalphas; alphacounter++)
    {
  for (frequcounter=0; frequcounter<nooffrequs; frequcounter++)
    {
 
      alphahat=alphacounter*deltaalpha+alphanull;
 
      frequofinterest=frequcounter*deltafrequ+frequnull;
   
      wf=wfini;


      counter_i=0;
      counter_ii=0;
    
    
    // ************* real timeprop

      floquwf.nullify();
       KHwf.nullify();
    timestep=complex<double>(real_timestep,0.0);
  
    no_of_timesteps=no_of_real_timesteps;





   
     
 
 for (ts=0; ts<no_of_timesteps; ts++)
      {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
      cout << "Real: " << ts << "  alphahat: " << alphahat <<  "  frequ: " << frequofinterest << " energy: " << real(complenerg) << " norm: " << wf.norm(g) <<  " xgrid:  " <<  g.ngps_x() <<  " ygrid:  "  <<           g.ngps_y() <<    endl;
      

      wf.propagate(timestep,time,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);   

  posofinterestKHx=posofinterestx;
	  gpKHx=g.xindex(posofinterestKHx);
	  

  posofinterestKHy=posofinteresty;
	  gpKHy=g.yindex(posofinterestKHy);
	 
	   for (targetenergcounter=0; targetenergcounter<nooftargetenergs; targetenergcounter++)
	    {

targetenergy=targetenergnull+deltatargetenergs*targetenergcounter;
	      for (xindex=0; xindex<g.ngps_x(); xindex++){
 for (yindex=0; yindex<g.ngps_y(); yindex++)
		
{

		  index=g_floqu.index(xindex,yindex,targetenergcounter);
		  posofinterestKHx=g.x(xindex);
                  posofinterestKHy=g.y(yindex)+alpha_y(time,0);
		  xindexKH=g_KH.xindex(g.x(xindex));
yindexKH=g_KH.yindex(g.y(yindex));
		  floquwf[index]=floquwf[index]+(timestep*exp(imagi*time*targetenergy))*wf[g.index(xindex,yindex,0)];
		  
		  if ((yindexKH>=0) && (yindexKH<g_KH.ngps_y())){
		    if ((xindexKH>=0) && (xindexKH<g_KH.ngps_x()) ){

		  	    
		      gplabx=g.xindex(posofinterestKHx);
                      gplaby=g.yindex(posofinterestKHy);
		      indexKH=g_KH.index(xindexKH,yindexKH,targetenergcounter);
		      KHwf[indexKH]=KHwf[indexKH]+(timestep*exp(imagi*time*targetenergy))*wf[g.index(gplabx,gplaby,0)];
		    };
		};
	    };
	    };
	    };
      if (counter_ii==obs_output_every)
      	{
	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	 
 
  
	      
	      fprintf(file_obser,"%e %e %e %e %e %e %e %e %e\n",
		      time,real(complenerg),imag(complenerg),wf.norm(g),wf.non_ionized(g,box),real(wf[g.index(gpKHx,gpKHy,0)]),imag(wf[g.index(gpKHx,gpKHy,0)]),wf.expect_x(g),wf.expect_y(g));
	      
	      counter_ii=0;
	    };
      
	  //	  correlation=wfini*wf;
	  correlation=wfeverythingeven*wf;
	  correlation*=g.delt_x()*g.delt_y();




	  //	  correlationex=wfexcited*wf;
	  correlationex=wfeverythingodd*wf;
	  correlationex*=g.delt_x()*g.delt_y();


 

      

	  fprintf(file_wfatgpoi,"%e %e %e %e %e %e %e\n",
		  time,real(wf[g.index(gpKHx,gpKHy,0)]),imag(wf[g.index(gpKHx,gpKHy,0)]),
		  real(correlation),imag(correlation),
		  real(correlationex),imag(correlationex));
	  

	 


	  fprintf(file_dipole_y,"%e %e\n",time,wf.expect_y(g));

	    if (counter_i==wf_output_every)
	    {
	    	  counter_i=0;
		  //	    if(ts>100 && ts<1100)
		  // {
	        wf.dump_to_file(g,file_wfdat,dumpingstepwidth);
		//};	  
		/*
               
if(ts>4000 && ts<5000)
{
                wf.dump_to_file(g,file_wf2dat,dumpingstepwidth);
};


	   
if(ts>10000 && ts<11000)
{
                wf.dump_to_file(g,file_wf3dat,dumpingstepwidth);
};


        
   
if(ts>20000 && ts<21000)
{
                wf.dump_to_file(g,file_wf4dat,dumpingstepwidth);
};



 
if(ts>30000 && ts<31000)
{
                wf.dump_to_file(g,file_wf5dat,dumpingstepwidth);
};




if(ts>40000 && ts<41000)
{
                wf.dump_to_file(g,file_wf6dat,dumpingstepwidth);
};



                
if(ts>45000 && ts<46000)
{
                wf.dump_to_file(g,file_wf7dat,dumpingstepwidth);
};




if(ts>48000 && ts<49000)
{
                wf.dump_to_file(g,file_wf8dat,dumpingstepwidth);
};

		*/	     
	    };
	  




         
    
	
	  

		  
  //if (ts=no_of_timesteps) goto loop;	  
      	  
	  
	};
};
	    
    wf.dump_to_file(g,file_wfdat,dumpingstepwidth);
    /*     wf.dump_to_file(g,file_wf2dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf3dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf4dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf5dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf6dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf7dat,dumpingstepwidth);

    wf.dump_to_file(g,file_wf8dat,dumpingstepwidth);
    */
   
      floquwf.dump_to_file(g_floqu,file_floquwf,dumpingstepwidth);
         KHwf.dump_to_file(g_KH,file_KHwf,dumpingstepwidth);
	   
       for (targetenergcounter=0; targetenergcounter<nooftargetenergs; targetenergcounter++)
	    {
	  KHnorm=0.0;
	  for (xindexKH=0; xindexKH<g_KH.ngps_x(); xindexKH++)
  for (yindexKH=0; yindexKH<g_KH.ngps_y(); yindexKH++)
	    {
{
	      indexKH=g_KH.index(xindexKH,yindexKH,targetenergcounter);
	      KHnorm+=conj(KHwf[indexKH])*KHwf[indexKH]*deltx*delty;
	    };

	    };
	  printf("%lf %e \n",real(KHnorm),imag(KHnorm));
	};
    
    
    


};  



   
     fclose(file_obser);
  
 fclose(file_wfdat);
/*  fclose(file_wf2dat);
 fclose(file_wf3dat);
 fclose(file_wf4dat);
 fclose(file_wf5dat);
 fclose(file_wf6dat);
 fclose(file_wf7dat);
 fclose(file_wf8dat);
*/
  fclose(file_KHwf);
  fclose(file_floquwf);
  fclose(file_wfatgpoi);
  fclose(file_dipole_y);
    fclose(file_nuclpot);


    cout << me << ": Hasta la vista, ... " << endl;



}


 double vvvvvvvecpot_y(double time, int me)
{
  
  double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;
  double ww, dur, n;

  frequref=12.0; // <--------- put same in alpha_x !!!

  vecpotampl = alphahat*frequref;


  n=200.0;


  ww=frequref/(2.0*n);
  dur=n*2*M_PI/frequref;
 
  if ((time>0.0) && (time<dur))
    {
      result=-vecpotampl*sin(ww*time)*sin(ww*time)*sin(frequref*time);
    };
  
  return result;
}             

double vecpot_y(double time, int me)
{
   
     double result;
     double vecpotampl, frequref,frequofinterest;
  double ramping, constperiod, downramp;

  frequref=12.0;
  frequofinterest=0.05;
  vecpotampl = alphahat*frequref;



   ramping=2.0*2*pi/frequref;
   constperiod=9551.0*2*pi/frequref;
   downramp=2.0*2*pi/frequref;


   // the first pulse with the reference frequency

   if (time<ramping)
     {
       result=-vecpotampl/ramping*time*cos(frequref*time) + vecpotampl/(frequref*ramping)*sin(frequref*time);
     };
   if ( (time >= ramping) && (time < ramping+constperiod) ) 
     {
       result=-vecpotampl*cos(frequref*time);
     }; 
   if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp) ) 
     {
       result=vecpotampl/downramp*(time-downramp-constperiod-ramping)*cos(frequref*time)-vecpotampl/(frequref*downramp)*sin(frequref*time);
     }; 
   if ( time > ramping+constperiod+downramp )
     {
       result=0.0;
     };



//    // second pulse (same ramping periods as first pulse)

   /*   vecpotampl=0.01/frequofinterest;

  if (time<ramping)
      {
        result+=-vecpotampl/ramping*time*cos(frequofinterest*time) + vecpotampl/(frequofinterest*ramping)*sin(frequofinterest*time);
      };
    if ( (time >= ramping) && (time < ramping+constperiod) ) 
      {
        result+=-vecpotampl*cos(frequofinterest*time);
      }; 
    if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp) ) 
      {
        result+=vecpotampl/downramp*(time-downramp-constperiod-ramping)*cos(frequofinterest*time)-vecpotampl/(frequofinterest*downramp)*sin(frequofinterest*time);
      }; 
    if ( time > ramping+constperiod+downramp )
      {
        result+=0.0;
      };

   */
   
  
   return result;
  
 
}  


double alpha_y(double time, int me)
{
  double result=0.0;
  double vecpotampl, frequref;
  double ramping, constperiod, downramp;

  frequref=12; // <-------- put same in vecpot_y !!!
  vecpotampl = alphahat*frequref;


  result=-vecpotampl/frequref*sin(frequref*time);

  return result;

}  




double vecpot_x(double time, int me)
{
   double result=0.0;
   return result;
}
   /* double ampl,frequ;
   double ehat,w,ww,phi,n,dur;

   // n=1000.0;
   //	 phi=0.0;
	 // w=3.14;
	  // ehat=2.5*w*w;
         frequ=3.14;
     ampl=2.5*frequ*frequ;
          double ramping, constperiod, downramp;
     constperiod=230.0*2*M_PI/frequ;
 ramping=35.0*2*M_PI/frequ;
 downramp=35.0*2*M_PI/frequ;

   if (time<ramping)
     {
       result=-(ampl/ramping)*time*cos(frequ*time)+ ampl/(frequ*ramping)*sin(frequ*time);
     
     };
   if ( (time >= ramping) && (time < ramping+constperiod) ) 
     {
       result=-ampl*cos(frequ*time);
     }; 
   if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp) ) 
     {
       result=ampl/downramp*(time-downramp-constperiod-ramping)*cos(frequ*time)-ampl/(frequ*downramp)*sin(frequ*time);
     }; 
   if ( time > ramping+constperiod+downramp )
     {
       result=0.0;
     };

   return result;

}


  ww=w/(2.0*n);
  dur=n*2*M_PI/w;
 
  if ((time>0.0) && (time<dur))
    {
      result=-(ehat)*sin(w*time)*sin(ww*time)*sin(ww*time);
    };
  
  return result;
 

} 

  

       
double vvecpot_z(double time, int me)
{
  double result=0.0;
  double ampl,frequ;
  double ramping, constperiod, downramp;

  frequ=0.375;
  ampl = 0.037589/frequ;

   constperiod=47.0*2*pi/frequ;
   ramping=2.0*2*pi/frequ;
   downramp=2.0*2*pi/frequ;


   if (time<ramping)
     {
       result=-ampl/ramping*time*cos(frequ*time) + ampl/(frequ*ramping)*sin(frequ*time);
     };
   if ( (time >= ramping) && (time < ramping+constperiod) ) 
     {
       result=-ampl*cos(frequ*time);
     }; 
   if ( (time >= ramping+constperiod) && (time < ramping+constperiod+downramp) ) 
     {
       result=ampl/downramp*(time-downramp-constperiod-ramping)*cos(frequ*time)-ampl/(frequ*downramp)*sin(frequ*time);
     }; 
   if ( time > ramping+constperiod+downramp )
     {
       result=0.0;
     };

   return result;

}  

*/
    
double vecpot_z(double time, int me)

{
  return 0.0;
}  

double scalarpotx(double x, double y, double z, double time, int me)
{
  double result;  
   double eps_n=0.4; 

//result = 1.0/sqrt(2.0*80*80/Mp+eps_n);
   result = -2.0/sqrt(x*x+eps_n);
//-1.0/(cosh(x*sqrt(2/Mp))*cosh(x*sqrt(2/Mp)));//1.0/sqrt(2.0*x*x/Mp+eps_n);
   
  return result;
}

double scalarpoty(double x, double y, double z, double time, int me)
{
  double result;
   double eps_n=0.4;
   // double eps_e=1.0;
   result = -2.0/sqrt(y*y+eps_n);
  //result =    -1.0/(sqrt((y-2.66)*(y-2.66)+eps_e));
 //   result +=   -1.0/(sqrt((y+2.66)*(y+2.66)+eps_e));
  
       return result;

}


double scalarpotz(double x, double y, double z, double time, int me)
{
  double result=0.0;
  return result;
}

double interactionpotxy(double x, double y, double z,double time, int me)
{
  

double  result=0.0; 

 double eps_e=0;
   result =    1.0/(sqrt(y*y+x*x+eps_e));
 //  result =    -1.0/(sqrt((y-x/sqrt(2.0*Mp))*(y-x/sqrt(2.0*Mp))+eps_e));
 // result +=   -1.0/(sqrt((y+x/sqrt(2.0*Mp))*(y+x/sqrt(2.0*Mp))+eps_e));
 
 //result =    -1.0/(sqrt((y-3)*(y-3)+eps_e));
 //result +=   -1.0/(sqrt((y+3)*(y+3)+eps_e));
 //result = -1.0/(cosh((y-x*sqrt(0.5/Mp)))*cosh((y-x*sqrt(0.5/Mp))));
 //result+=-1.0/(cosh((y+x*sqrt(0.5/Mp)))*cosh((y+x*sqrt(0.5/Mp))));
       return result;

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
       double ampl=30.0; // switch imaginary potential on


  if (ampl>1.0)
    {
      //      x=((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x())
      //	*((double) xindex + 0.5  - 0.5*g.ngps_x())/(0.5*g.ngps_x());

    x=((double) xindex + 0.5)/(1.0*g.ngps_x())
	*((double) xindex + 0.5)/(1.0*g.ngps_x());
  
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
        double ampl=30.0; // switch imaginary potential on


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



// dft is not used
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &v_null, const wavefunction &v_eins)
{
  double result;
  result=0.0;
  return result;
}

