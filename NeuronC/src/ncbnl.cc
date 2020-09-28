/* research simulation of A-> B process with discrete time steps */

#include <stdlib.h>
//#include <math.h>
#include "nconst.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
}
#endif

int cumrand;
int rseed=13245;

int binomdev(double pp, register int n);

int  main(){
        double tau=.55;
        double time, timestep=0.5;
        double subtime, subtimestep;
        double p, p_corr, p_small;                       
        double a1 ,a2, a_uncor;
                                
        subtimestep= tau/100;
	        
        p_corr= 1-exp(-timestep/tau);           
        p = timestep/tau;
        p_small = subtimestep/tau;
        
        a1= 1e6;
        a2= 1e6;
        a_uncor= 1e6;
        
        for (time = 0; time <=5; time+= timestep){
                a1 -= binomdev(p_corr, (int)a1);
                a_uncor -= binomdev(p, (int)a_uncor);
                for (subtime=time; subtime <= time+timestep; 
						subtime+= subtimestep){
                        a2 -= binomdev(p_small, (int)a2);                                       
                }
                printf("%g %g %g %g\n", time+timestep/2, a1, a_uncor, a2 );
        }
        exit(0);
}


