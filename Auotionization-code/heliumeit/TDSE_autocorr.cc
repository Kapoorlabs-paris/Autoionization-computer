 
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
  FILE *file_realdftpotdat,*file_densityapprox,*file_auto;
  FILE *file_realpot,*file_realpotnophase,*file_tddftpot, *file_kohnshamorbital, *file_wfdftexact,*file_wfdftexactnophase,*file_wfdftapprox;
  FILE *file_obser, *file_obser_imag,*file_alphapart,*file_justalpha;
  FILE *file_reading,*file_reading2,*file_reading3,*file_reading4;    // *** create some files with appropriate appendices
  char string_wfdat[]=       "/data/vk060/wf_tripleafterlaserhigh.dat";


  char string_obser[]=       "/data/vk060/observ_tripleafterlaserhigh.dat";
  char string_obser_imag[]=  "/data/vk060/observimag_tripleafterlaserhigh.dat";
 
  char string_reading[]=     "/data/vk060/wf_tripleafterlaser.dat";
  char string_reading2[]=     "/data/vk060/wf_ground3500.dat";

  char string_reading3[]=     "/data/vk060/wf_volkovandion.dat";
  char string_reading4[]=     "/data/vk060/wfexcited_volkovandion.dat";

char string_auto[]=  "/data/vk060/auto_tripleafterlaserhigh.dat";


  file_auto = fopen(string_auto,"w");

  file_wfdat = fopen(string_wfdat,"w");


  file_obser = fopen(string_obser,"w");
  file_obser_imag = fopen(string_obser_imag,"w");
   file_reading = fopen(string_reading,"r");
file_reading2 = fopen(string_reading2,"r");



 file_reading3 = fopen(string_reading3,"r");
file_reading4 = fopen(string_reading4,"r");


   long index, rrindex, xindex, yindex,zindex,index2;
 complex<double> imagi(0.0,1.0);
  double deltx=0.5;
double delty=0.5;
double deltz=1;
  long ngpsx=3000;
long ngpsy=3000;
long ngpsxbig=3500;
long ngpsybig=3500;
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




   
  
 
  
  
  // *** declare rest of variables ***
  double imag_timestep=0.2;
  double real_timestep=0.05;
  long   no_of_imag_timesteps=100;
  long   no_of_real_timesteps=7400;//6800*2;
   long  freetime=0;//6800*2;
int    obs_output_every=1;
  long   wf_output_every=300000; 
  int    dumpingstepwidth=1;
  int    vecpotflag=1;   // don't touch this
  int    box=50;          //was 50
  double masses[]={1.0,1.0};    // don't touch this
  double charge=0.0;    // don't touch this
  double epsilon=0.00001;
  hamop interactionhamil(gone,vecpot_x,vecpot_y,vecpot_z,interactionpotxy,scalarpotytwo,scalarpotztwo,interactionpotxy,imagpotxtwo,imagpotytwo,field,dftpot);
  hamop hamilton(g,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotx,imagpoty,field,dftpot);
 
  hamop hamiltonbig(gbig,vecpot_x,vecpot_y,vecpot_z,scalarpotx,scalarpoty,scalarpotz,interactionpotxy,imagpotxbig,imagpotybig,field,dftpot);

 
 
 hamop hamiltonone(gbigone,vecpot_x,vecpot_y,vecpot_z,scalarpotxtwo,scalarpotytwo,scalarpotztwo,interactionpotxytwo,imagpotxtwo,imagpotytwo,field,dftpot);
 
  wavefunction wf(g.ngps_x()*g.ngps_y()*g.ngps_z());
   wavefunction wfex(g.ngps_x()*g.ngps_y()*g.ngps_z());
   wavefunction wfex2(g.ngps_x()*g.ngps_y()*g.ngps_z());
   wavefunction wfexcited(g.ngps_x()*g.ngps_y()*g.ngps_z());
    wavefunction wfexcited2(g.ngps_x()*g.ngps_y()*g.ngps_z());
      
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

wavefunction wfreadexcited(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfreadexcited2(g.ngps_x()*g.ngps_y()*g.ngps_z());

wavefunction wfreadexcited3(g.ngps_x()*g.ngps_y()*g.ngps_z());
wavefunction wfreadexcited4(g.ngps_x()*g.ngps_y()*g.ngps_z());




   wavefunction wfini(g.ngps_x()*g.ngps_y()*g.ngps_z());
  
 wavefunction wfini2(g.ngps_x()*g.ngps_y()*g.ngps_z());


wavefunction wfeverythingeveniny(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
  wavefunction wfeverythingoddiny(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());

wavefunction wfeverythingeven(g.ngps_x()*g.ngps_y()*g.ngps_z());
  wavefunction wfeverythingodd(g.ngps_x()*g.ngps_y()*g.ngps_z());
 
    wavefunction correlation(g.ngps_x()*g.ngps_y()*g.ngps_z()) ;
  wavefunction correlationex(g.ngps_x()*g.ngps_y()*g.ngps_z()) ;
   
 
 
//double posofinteresty=2;
// double posofinterestx=3;
    long posofinterest=2500;
  double posofinterestKHx,posofinterestKHy;

  long gpKHx,gplabx,gpKHy,gplaby;

  
long noofalphas=1;
  double deltaalpha=0.02;
  double alphanull=0.1;
  
 long alphacounter;
  
  long nooffrequs=1;
  double deltafrequ=0.006;
  double frequnull=1.7;//1.356;//1.356;
  long frequcounter;


  long nooftargetenergs=1;
  double deltatargetenergs=0.002;
  double targetenergnull=-0.881875;
  long targetenergcounter;
  double targetenergy;

  complex<double> timestep;

complex<double> timestep2;

  double time=0.0;
double time2=0.0;  
  long counter_i=0;
  long counter_ii=0;
   long counter_iii=0;
long counter_iv=0;
long counter_v=0;
  
  
  // initialization
  long outputofinterest=0;
  
  
wf.init(g,1,2.0,0.0,0.0);
      // ground


//wfground.init(gbig,99,0.1,0.0,0.0,file_reading2,outputofinterest);
//    wfread.init(g,99,0.1,0.0,0.0,file_reading,outputofinterest);
  //  wfbig.nullify();
 //wfinibig.nullify();



   // wfbig.regrid(gbig,g,wfread);
//wfinibig.regrid(gbig,g,wfread);


//wfchannel.init(gone,99,0.1,0.0,0.0,file_reading3,outputofinterest);

// wfbigchannel.nullify();

//wfbigchannel.regrid(gbigone,gone,wfchannel);




//wfsecondchannel.init(gone,99,0.1,0.0,0.0,file_reading4,outputofinterest);

 //wfsecondbigchannel.nullify();

//wfsecondbigchannel.regrid(gbigone,gone,wfsecondchannel);







    wf*=1.0/sqrt(wf.norm(g));
  // wfground=wf;
  
  
     
 wfeverythingeveniny.init(gbig,4,1.0,1.375,0.0);
  wfeverythingeveniny*=1.0/sqrt(wfeverythingeveniny.norm(gbig));

  wfeverythingoddiny.init(gbig,5,1.0,0.684,0.0);
   wfeverythingoddiny*=1.0/sqrt(wfeverythingoddiny.norm(gbig));

;
  
 //wf=wfground;

   fclose(file_reading);
 
   fclose(file_reading2);


    
  cout << "norm wf    : " << wf.norm(g) <<  "\n";

  

  wavefunction staticpot_x(g.ngps_x());
  staticpot_x.calculate_fixed_potential_array_x(g,hamilton,0.0,me);
  
  wavefunction staticpot_y(g.ngps_y());
  staticpot_y.calculate_fixed_potential_array_y(g,hamilton,0.0,me);
  
  
  
  wavefunction staticpot(g.ngps_x()*g.ngps_y()*g.ngps_z());
  staticpot.calculate_fixed_potential_array(g,hamilton,0.0,me);
   

  

  
     wavefunction staticpot_xy(g.ngps_z()*g.ngps_y()*g.ngps_x());
  staticpot_xy.calculate_fixed_potential_array_xy(g,hamilton,0.0,me);




 wavefunction staticpot_xbig(gbig.ngps_x());
  staticpot_xbig.calculate_fixed_potential_array_x(gbig,hamiltonbig,0.0,me);

  wavefunction staticpot_ybig(gbig.ngps_y());
  staticpot_ybig.calculate_fixed_potential_array_y(gbig,hamiltonbig,0.0,me);



  wavefunction staticpotbig(gbig.ngps_x()*gbig.ngps_y()*gbig.ngps_z());
  staticpotbig.calculate_fixed_potential_array(gbig,hamiltonbig,0.0,me);





     wavefunction staticpot_xybig(gbig.ngps_z()*gbig.ngps_y()*gbig.ngps_x());
  staticpot_xybig.calculate_fixed_potential_array_xy(gbig,hamiltonbig,0.0,me);



  
  wavefunction staticpotdft(gtwo.ngps_x()*gtwo.ngps_y()*gtwo.ngps_z());
  staticpotdft.calculate_fixed_potential_array(gtwo,hamiltontwo,0.0,me);




  wavefunction staticpotdftbig(gbigone.ngps_x()*gbigone.ngps_y()*gbigone.ngps_z());
  staticpotdft.calculate_fixed_potential_array(gbigone,hamiltonone,0.0,me);

  


 complex<double> dftapproxenerg;
 complex<double> complenerg,complenerg0,dftenerg,complenerg2,complenerg3,complenerg4,complenerg5,complenerg6,complenerg7,complenerg8,complenerg9,complenerg10,complenerg11,complenerg12,complenerg13,complenerg14,complenerg15,complenerg16,complenerg17,complenerg18,complenerg19,complenerg20,complenerg21;
  complex<double> groundstatepop, excitedstatepop,excitedstate1dpop,groundstatepopdft;
 complex<double> KHnorm;
  complex<double> complenergex,evencorrelation,oddcorrelation;
  complex<double>correlationint,correlationexint,correlationdftint,correlationexdftint;
  complex<double>correlationdftapproxint,correlationexdftapproxint,autocorrelation;
  


  // ************* imag timeprop
  long ts;
  long no_of_timesteps=no_of_imag_timesteps;
  for (ts=0; ts<no_of_timesteps; ts++)
    {   
                                                                                  
      cout << "Imag: " << ts <<"  "  <<  " wf  "  << real(complenerg0) <<   endl; 
           cout<< real(dftapproxenerg) << endl;
      counter_i++;
      counter_ii++;
      timestep=complex<double>(0.0*real_timestep,-1.0*imag_timestep);  
      time=-imag(timestep*(complex<double>)(ts));
      
      // and now the actual propagation

 

  wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);

wf*=1.0/sqrt(wf.norm(g));

//wfexcited50=wfexcited50-((wfexcitedheliumplus*wfexcited50)*g.delt_x()*g.delt_y())*wfexcitedheliumplus;

/* 
     wf.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);
 

//wf=wf-((wfexcitedheliumplus*wf)*g.delt_x()*g.delt_y())*wfexcitedheliumplus;


      wf*=1.0/sqrt(wf.norm(g));



      wfexcited.propagate(timestep,0.0,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);  

wfexcited=wfexcited-((wf*wfexcited)*g.delt_x()*g.delt_y())*wf;

 wfexcited*=1.0/sqrt(wfexcited.norm(g));
 
 */









     
      

      
      if (counter_ii==obs_output_every)
      	{
			 complenerg0=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	
	  
	  	  fprintf(file_obser_imag,"%li %.14le %.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
			  ts,real(complenerg),imag(complenerg),wf.norm(g),wf.non_ionized(g,box),wf.sing_ionized(g,box),wf.doub_ionized(g,box),wf.expect_x(g),wf.expect_y(g));
	  counter_ii=0;
	 
	};
      
    };
  
  
  fclose(file_obser_imag);
  
  
 wfini=wf;
 
 

  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
 	 
       
    

    for (alphacounter=0; alphacounter<noofalphas; alphacounter++)
   {
   for (frequcounter=0; frequcounter<nooffrequs; frequcounter++)
   {
    
      counter_i=0;
      counter_ii=0;
counter_iv=0;
  counter_v=0;   

 alphahat=alphacounter*deltaalpha+alphanull;
 
      frequref=frequcounter*deltafrequ+frequnull;
   
      wf=wfini;
  
  // ************* real timeprop
  
           
  timestep=complex<double>(real_timestep,0.0);
  
  no_of_timesteps=no_of_real_timesteps;
  
  
 



  for (ts=0; ts<no_of_timesteps; ts++)
    {
      counter_i++;
      counter_ii++;
      time=real(timestep*(complex<double>)(ts));
        cout << "Real: " << ts << " alphahat: " << alphahat <<  " frequ: " << frequref << "  "  << " norm wf : " << wf.norm(g) << "  "  <<  " wfdftexact  "  << wfdftexact.norm(gone) << "  "  <<  " wfdftexactnophase  "  << wfdftexactnophase.norm(gone) << "\n";
    


      wf.propagate(timestep,time+0.5*real(timestep),g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);   
 
 

   if (counter_ii==obs_output_every)
      	{
	  complenerg=wf.energy(0.0,g,hamilton,me,masses,staticpot_x,staticpot_y,staticpot_xy,charge);
	       
	       groundstatepop=wf*wfini*g.delt_x()*g.delt_y();
excitedstatepop=wf*wfexcited*g.delt_x()*g.delt_y();
 
	  	  fprintf(file_obser,"  %.14le  %.14le %.14le %.14le %.14le %.14le %.14le %.14le  %.14le   \n ",
	    (time+real(timestep)),
	    real(complenerg),
		  imag(complenerg),
		  wf.norm(g),

 
		  
		
		  wf.expect_x(g),
		  wf.expect_y(g),
		hamilton.vecpot_x(time,me),
			
		 
		  real(conj(excitedstatepop)*excitedstatepop),
			  
	
                  real(conj(groundstatepop)*groundstatepop)	  

		  );
	  
	  counter_ii=0;
	};
    

};
 
 
 };
};

 


   wfini2=wf;
    
    wf=wfini2;
    timestep2=complex<double>(real_timestep,0.0);

//wfinibig=wfinibig-((wfground*wfinibig)*gbig.delt_x()*gbig.delt_y())*wfground;
//wfbig=wfbig-((wfground*wfbig)*gbig.delt_x()*gbig.delt_y())*wfground;


  wfini2.dump_to_file(g,file_wfdat,dumpingstepwidth);

complenerg=(wfbig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge));

//wfbig.regrid(gbig,g,wfini2);
//wfinibig.regrid(gbig,g,wfini2);    
//wfprobig=wfbig;
//wfbig=wfprobig;
  for (ts=no_of_timesteps; ts<no_of_timesteps+freetime; ts++)
    {
      counter_iv++;
      counter_v++;
      time2=real(timestep2*(complex<double>)(ts));
      
      cout << "Real: " << ts << " "  << " corr : " << evencorrelation << "  " << " odd :" << oddcorrelation <<      endl;  
 



//wf.propagate(timestep2,time2,g,hamilton,me,vecpotflag,staticpot_x,staticpot_y,staticpot_xy,charge);  
 wfbig.propagate(timestep2,time2,gbig,hamiltonbig,me,vecpotflag,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge);   
//wfprobig=wfbig-((wfbig*wfinibig)*g.delt_x()*g.delt_y())*wfinibig;

 //wf.oddinz(g);
//wf*=1.0/sqrt(wf.norm(g));

// wfchannel.propagatedft(timestep2,time2,gbigone,hamiltonone,me,vecpotflag,staticpotdftbig,0.0*(staticpot),charge);

//wfsecondchannel.propagatedft(timestep2,time2,gbigone,hamiltonone,me,vecpotflag,staticpotdftbig,0.0*(staticpot),charge);


complenerg=(wfbig.energy(0.0,gbig,hamiltonbig,me,masses,staticpot_xbig,staticpot_ybig,staticpot_xybig,charge));

//evencorrelation=wfbig.autocorrelationsecond(gbig,wfeverythingeveniny,posofinterest);  
    

  //oddcorrelation=wfbig.autocorrelationsecond(gbig,wfeverythingoddiny,posofinterest);


//evencorrelation=wfbig.autocorrelationnew(gbig,wfeverythingeveniny,wfinibig); 

//oddcorrelation=wfbig.autocorrelationnew(gbig,wfeverythingoddiny,wfinibig);



evencorrelation=wfeverythingeveniny*wfbig;
evencorrelation*=gbig.delt_x()*gbig.delt_y();

oddcorrelation=wfeverythingoddiny*wfbig;
oddcorrelation*=gbig.delt_x()*gbig.delt_y();



//evencorrelation=wfbig.tsurff(gbig,gbigone,wfbigchannel,1.37473,posofinterest);
//oddcorrelation=wfbig.tsurff(gbig,gbigone,wfsecondbigchannel,0.683633,posofinterest);


   if (counter_iv==obs_output_every)
  	{
	 
//	   if (counter_iv==wf_output_every)
 //       {

   //       counter_iv=0;
     //     wfbig.dump_to_file(gbig,file_wfdat,dumpingstepwidth);
//  wfprobig.dump_to_file(gbig,file_wfprodat,dumpingstepwidth);


 
	     
	

	  	  fprintf(file_auto, " %.14le %.14le %.14le %.14le %.14le   \n ",
	    (time2+real(timestep)),
	    real(evencorrelation),
		  imag(evencorrelation),real(oddcorrelation),imag(oddcorrelation)
		 
		 
			
	
                  

		  );
	  
	  counter_iv=0;
	};



 if (counter_v==wf_output_every)
        {

         counter_v=0;
          wfbig.dump_to_file(g,file_wfdat,dumpingstepwidth);

};

};      	

    




    fclose(file_obser);
   
   fclose(file_auto);

    fclose(file_wfdat);

    cout << me << ": Hasta la vista, ... " << endl;

 

 }



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
 
  double constperiod=80*2*M_PI/frequ;//392.0*2*M_PI/frequ;
double delay=40.0*2*M_PI/frequ;  

double dur=ramping+constperiod+downramp;

  double ampl = alphahat*frequ;
 double n=100.0;
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


// dft is not used
double dftpot(grid g, double x, double y, double z, double time, int me, 
	      const fluid &v_null, const wavefunction &v_eins)
{
  double result;
  result=0.0;
  return result;
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


