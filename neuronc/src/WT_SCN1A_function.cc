/*
 *  SCN1A_function.cc
 *  
 *
 *  Created by Colleen E. Clancy on Wed Jan 29 2003.
 *  Copyright (c) 2003 Columbia University. All rights reserved.
 *
 */ 

#ifndef WT_SCN1A_function_cc
#define WT_SCN1A_function_cc 

#include "SCN1A_global_function.h"  		//header file for Na+ channel parameters



//**************************************************************************************
//	FUNCTION TO COMPUTE SCN1A NA CURRENT--Markov Model INCEPTION 09/17/2002
//**************************************************************************************
	 double WT_SCN1A_function ()
	
{	
	
	Ena=(r*T/F)*log(Naout/Nain);

double a11, a12, a13, b11, b12, b13, a2, a3, a4, b4, b1, b2, b3, a5, b5, mu1, mu2, aa2, bb2, bb3, aa3;
double err1, err2, err3,err4, err5, err6, err7, err8, err9, err10, err11, err12, err13, err14; 

double  y1, z1,  x1, x2, y2, z2, e1, e2,  q1, q2,  w1, w2,  d1, d2, c1, c2, f1, j2, f2, tap1, tap2, jj1, u2, u1, p1, p2, t1, zz1, zz2;

double dny, dnz, dne1, dne, dne2, dny1, dnz1, dnx, dnx1,  dnq1, dnq,dnw, dnw1, dnd, 
 dnd1, dnc, dnc1, dnu, dnu1,  dnp, dnp1,  dntap, dntap1,  dnf1, dnj1, dnf, dnj, dnzz, dnzz1;

int i4;

 Nain=10.0; 
 Naout= 140.0;

			


if (which_channels ==1)   // WT rates
{
    a11= 2.802/(0.21*exp(-v/17.0)+0.23*exp(-v/150));
    a12= 2.802/(0.23*exp(-v/15.0)+0.25*exp(-v/150));
    a13= 2.802/(0.25*exp(-v/12.0)+0.27*exp(-v/150));
  

     b11= 0.4*exp(-v/20.3);
     b12= 0.4*exp(-(v-5)/20.3);
     b13= 0.4*exp(-(v-10)/20.3)/4.5;
     a3= (3.7933e-7*exp(-v/7.6))*3;
     b3= (0.0084+.00002*v);
    a2= ((9.178*exp(v/29.68))/4.5);
     b2= ((a13*a2*a3)/(b13*b3));	
     a4= (a2/100)*1.5;
      b4 = a3/5;
     a5= (a2/95000)*80.0;
      b5 = (a3/30)/10.0;

mu1= 0.0;
mu2= 0.0;


if (t<dt)
     {
     
ina = 0.0;	
			//initialize the Na+ current
Na_MC3=0.997014;
Na_MO=3.40959e-13;
Na_MI=2.30533e-11;
Na_MIs=3.94247e-14;
Na_MC2=0.000173255;
IC3=0.00281263;
IC2=4.88762e-07;
Na_MIs2=2.27107e-16;
M_C3=0;
M_C2=0;
M_C=0;
M_O=0;
M_IF=0;


//-65


/*Na_MC3=0.439163;
Na_MO=2.15859e-06;
Na_MI=0.000188108;
Na_MIs=0.000517741;
Na_MC2=0.0125585;
IC3=0.528912;
IC2=0.015125;
Na_MIs2=0.00337774;
M_C3=0;
M_C2=0;
M_C=0;
M_O=0;
M_IF=0;

*/

   Na_MC= 1-(Na_MO+Na_MI+Na_MC2+Na_MC3+M_C3+Na_MIs+M_C2+M_C+M_O+IC3+IC2+Na_MIs2);
   
        
 	}

} 

else if (which_channels ==2)   // RH rates
{
    a11= 2.802/(0.21*exp(-v/17.0)+0.23*exp(-v/150));
    a12= 2.802/(0.23*exp(-v/15.0)+0.25*exp(-v/150));
    a13= 2.802/(0.25*exp(-v/12.0)+0.27*exp(-v/150));
  

     b11= 0.4*exp(-v/20.3);
     b12= 0.4*exp(-(v-5)/20.3);
     b13= 0.4*exp(-(v-10)/20.3)/4.5;
     a3= (3.7933e-7*exp(-v/7.6))*3;
     b3= (0.0084+.00002*v);
    a2= ((9.178*exp(v/29.68))/4.5);
     b2= ((a13*a2*a3)/(b13*b3));	
     
     a4= (a2/100)*1.5;
      b4 = a3/5;
      
     a5= (a2/95000)*80.0;
      b5 = (a3/30)/10.0;

// flicker inactivation rates

aa2= ((9.178*exp(v/29.68))/4.5);
bb3= (0.0084+.00002*v);
aa3= (3.7933e-7*exp(-v/7.7))*83;

bb2= ((a13*aa2*aa3)/(b13*bb3));

mu1= 2e-5;
mu2= 2e-4;


if (t<dt)
     {
     
ina = 0.0;				//initialize the Na+ current

Na_MC3=0.908416;
Na_MO=1.4419e-14;
Na_MI=1.11474e-12;
Na_MIs=3.65127e-16;
Na_MC2=7.04203e-05;
IC3=0.000665298;
IC2=5.15739e-08;
Na_MIs2=4.02847e-19;
M_C3=0.0908416;
M_C2=7.04203e-06;
M_C=1.5221e-10;
M_O=1.4419e-15;
M_IF=4.94622e-15;


   Na_MC= 1-(Na_MO+Na_MI+Na_MC2+Na_MC3+M_C3+Na_MIs+M_C2+M_C+M_O+IC3+IC2+Na_MIs2+M_IF);
   
        
 	}

} 


  	if (t>=60000.0-dt&&Get_steady_state==1)
 	{
 	
 	cout<<"Na_MC3="<<Na_MC3<<";"<<endl;
 	cout<<"Na_MO="<<Na_MO<<";"<<endl;
 	cout<<"Na_MI="<<Na_MI<<";"<<endl;
 	cout<<"Na_MIs="<<Na_MIs<<";"<<endl;
 	cout<<"Na_MC2="<<Na_MC2<<";"<<endl;
 	cout<<"IC3="<<IC3<<";"<<endl;
 	cout<<"IC2="<<IC2<<";"<<endl;
	cout<<"Na_MIs2="<<Na_MIs2<<";"<<endl;
 	
 	cout<<"M_C3="<<M_C3<<";"<<endl;
 	cout<<"M_C2="<<M_C2<<";"<<endl;
 	cout<<"M_C="<<M_C<<";"<<endl;
 	cout<<"M_O="<<M_O<<";"<<endl;
        cout<<"M_IF="<<M_IF<<";"<<endl;
 	
 	exit (0);
 	}
        
  
   
 if (t>=dt)
   {	
     y1=Na_MO;     z1=Na_MI;	x1= Na_MC2;	q1= Na_MC3;	w1= Na_MIs;
     d1= M_C3; 	   c1= M_C2;	u1= M_C;	p1= M_O;	tap1= Na_MC;
     jj1= IC3;	   f1= IC2;	e1 = Na_MIs2;	zz1= M_IF;
     
     
     //UPPER STATES
     dny=((Na_MC*a13+Na_MI*b2+M_O*mu2)- (Na_MO*(b13+a2+mu1)))*dt;
     dnz=((IC2*a12+Na_MO*a2+Na_MC*b3+Na_MIs*b4)- (Na_MI*(b2+a3+a4+b12)))*dt;
     dnx= ((IC2*a3+Na_MC*b12+Na_MC3*a11+M_C2*mu2)-(Na_MC2*(a12+b11+mu1+b3)))*dt;
     dnq= (IC3*a3+Na_MC2*b11+M_C3*mu2-(Na_MC3*(a11+mu1+b3)))*dt;
     dnw= ((Na_MI*a4+Na_MIs2*b5)-Na_MIs*(b4+a5))*dt;
     dntap= (Na_MC2*a12+Na_MO*b13+M_C*mu2+Na_MI*a3-(Na_MC*(a13+b12+mu1+b3)))*dt;
     dnj= ((Na_MC3*b3+IC2*b11)-IC3*(a11+a3))*dt;
     dnf= ((Na_MC2*b3+IC3*a11+Na_MI*b12)-IC2*(a12+b11+a3))*dt;
     dne= (Na_MIs*a5-Na_MIs2*b5)*dt;
     
     //LOWER STATES
     dnd= (Na_MC3*mu1+M_C2*b11-(M_C3*(a11+mu2)))*dt;
     dnc=  (Na_MC2*mu1+M_C3*a11+M_C*b12-(M_C2*(a12+mu2+b11)))*dt;
     dnu= (Na_MC*mu1+M_C2*a12+M_O*b13+M_IF*aa3-(M_C*(a13+mu2+b12+bb3)))*dt;
     dnp= ((M_C*a13+Na_MO*mu1+M_IF*bb2)-(M_O*(b13+mu2+aa2)))*dt;
     dnzz= ((M_O*aa2+M_C*bb3)- (M_IF*(bb2+aa3)))*dt;
     
 //  cout<<dny+dnz+dnx+dnq+dnw+dntap+dnj+dnf+dnd+dnc+dnu+dnp+dne+dnzz<<endl; 		//Use this as a test that derivatives sum to zero
     
   
    
   
   y2=Na_MO+dny;
    z2=Na_MI+dnz;
    x2= Na_MC2+dnx;
    q2= Na_MC3+dnq;
    w2= Na_MIs+dnw;
    tap2 = Na_MC+dntap;
    j2= IC3+dnj;
    f2= IC2+dnf;
    e2 = Na_MIs2+dne;
    
       d2= M_C3+dnd;
       c2= M_C2+dnc;
       u2= M_C+dnu;
       p2= M_O+dnp;
       zz2= M_IF+dnzz;
       
    
   
    
  // Na_MC= 1-(y2+z2+x2+q2+w2+d2+c2+u2+p2);

err1= y2-y1;	err2= z2-z1;	err3=x2-x1;	err10= tap2-tap1;
err4= q2-q1;	err5= w2-w1;	err6=d2-d1;	err11= j2-jj1;
err7= c2-c1;	err8= u2-u1;	err9= p2-p1;	err12= f2- f1;
err13= e2-e1;	err14= zz2-zz1;

	dny1= dny;
	dnz1=dnz;
	dnx1= dnx;		//UPPER STATES
	dnq1= dnq;
	dnw1=dnw;
	dntap1=dntap;
	dnj1= dnj;
	dnf1= dnf;
	dne1= dne;
	
	
	dnd1= dnd;
	dnc1= dnc;		//LOWER STATES
	dnu1= dnu;
	dnp1= dnp;
        dnzz1= dnzz;
	
i4=0;
while (((err1>1e-5)||(err1<-1e-5)||(err2>1e-5)||(err2<-1e-5)||(err3>1e-5)
||(err3<-1e-5)||(err4>1e-5)||(err4<-1e-5)||(err5>1e-5)||(err5<-1e-5)||
(err6>1e-5)||(err6<-1e-5)||(err7>1e-5)||(err7<-1e-5)||(err8>1e-5)||
(err8<-1e-5)||(err9>1e-5)||(err9<-1e-5)||(err10>1e-5)||(err10<-1e-5)||
(err11>1e-5)||(err11<-1e-5)||(err12>1e-5)||(err12<-1e-5)||(err13>1e-5)||(err13<-1e-5)||(err14>1e-5)||(err14<-1e-5))
&&(i4<40))
 {
	y1=y2; 	 	z1=z2;	 	 x1=x2;		q1=q2;		w1= w2;
	d1= d2;		c1= c2;		u1= u2;		p1= p2;		tap1= tap2;
	jj1= j2;	f1= f2;		e1= e2;		zz1= zz2;
	
    dny=((tap1*a13+z1*b2+p1*mu2)- (y1*(b13+a2+mu1)))*dt;
    dnz=((f1*a12+y1*a2+tap1*b3+w1*b4)- (z1*(b2+a3+a4+b12)))*dt;
    dnx=  (((f1*a3+tap1*b12+q1*a11+c1*mu2)-(x1*(a12+b11+mu1+b3)))*dt);
    dnq= ((jj1*a3+x1*b11+d1*mu2-q1*(a11+mu1+b3))*dt);
    dnw= ((z1*a4+e1*b5)-w1*(b4+a5))*dt;
    dntap= (x1*a12+y1*b13+u1*mu2+z1*a3-(tap1*(a13+b12+mu1+b3)))*dt;
    dnj= ((q1*b3+f1*b11)-jj1*(a11+a3))*dt;
    dnf= ((x1*b3+jj1*a11+z1*b12)-f1*(a12+b11+a3))*dt;
    dne= (w1*a5- e1*b5)*dt;
    
  // cout<<dntap<<setw(15)<<t<<endl;
    
    
    dnd= (q1*mu1+c1*b11-(d1*(a11+mu2)))*dt;
    dnc= (x1*mu1+d1*a11+u1*b12-(c1*(a12+mu2+b11)))*dt;
    dnu= (tap1*mu1+c1*a12+p1*b13+zz1*aa3-(u1*(a13+mu2+b12+bb3)))*dt;
    dnp= (u1*a13+y1*mu1+zz1*bb2-(p1*(mu2+b13+aa2)))*dt;
    dnzz= (p1*aa2+u1*bb3-(zz1*(aa3+bb2)))*dt;
    
    
    
    dny= (dny+dny1)/2;
    dnz= (dnz+dnz1)/2;
    dnx= (dnx+dnx1)/2;
    dnq= (dnq+dnq1)/2;
    dnw= (dnw+dnw1)/2;
    dntap= (dntap+dntap1)/2;
    dnj= (dnj+dnj1)/2;
    dnf= (dnf+dnf1)/2;
    dne= (dne+dne1)/2;
    
    dnd= (dnd+dnd1)/2;
    dnc= (dnc+dnc1)/2;
    dnu= (dnu+dnu1)/2;
    dnp= (dnp+dnp1)/2;
    dnzz= (dnzz+dnzz1)/2;
    
     y2= Na_MO+dny;
     z2= Na_MI+dnz;
     x2= Na_MC2+dnx;
     q2= Na_MC3+dnq;
     w2= Na_MIs+dnw;
     tap2= Na_MC+dntap;
     j2= IC3+dnj;
     f2= IC2+dnf;
     e2= Na_MIs2+dne;
     
     d2= M_C3+dnd;
     c2= M_C2+dnc;
     u2= M_C+dnu;
     p2= M_O+dnp;
     zz2= M_IF+dnzz;
     
     dny1=dny;
     dnz1=dnz;
     dnx1= dnx;
     dnq1=dnq;
     dnw1= dnw;
     dntap1= dntap;
     dnf1= dnf;
     dnj1= dnj;
     dne1= dne;
     
     dnd1= dnd;
     dnc1= dnc;
     dnu1= dnu;
     dnp1= dnp;
     dnzz1= dnzz;
     
    // Na_MC= 1-(y2+z2+x2+q2+w2+d2+c2+u2+p2);
     
	err1=y2-y1;
 	err2= z2-z1;
 	err3= x2-x1;
 	err4= q2-q1;
 	err5= w2-w1;
 	err6= d2-d1;
 	err7= c2-c1;
 	err8= u2-u1;
 	err9= p2-p1;
 	err10= tap2-tap1;
 	err11= j2-jj1;
 	err12= f2-f1;
 	err13= e2-e1;
        err14= zz2-zz1;
 	
 	i4++;
	
   
 }
 
 
 if (i4<40)	
    { 
   Na_MO=y2; 	 Na_MI=z2;	 Na_MC2= x2;	Na_MC3= q2;	Na_MIs=w2;
   M_C3= d2;   	 M_C2= c2;	 M_C= u2;	M_O= p2;	Na_MC=tap2;
   IC3= j2;	IC2=f2;		Na_MIs2= e2; 	M_IF= zz2;
   
  //  Na_MC=1-(Na_MO+Na_MI+ Na_MC2+Na_MIs+Na_MC3+M_C3+M_C2+M_C+M_O);
  
    }
    
   
    
  else 
  { 
   cout<<"error_ina"<<setw(15)<<i4<<setw(15)<<t<<endl;
   exit(0);
    }
  
    
  
   
   } 
    
    M_OP= Na_MO+M_O;
   M_OP= (M_OP*1.0);
  Na_OP= Na_O*0;
    
    OPs= M_OP+Na_OP;
     ina= 10.5*(OPs)*(v-Ena);		
    
	
	return ina;
}


#endif


