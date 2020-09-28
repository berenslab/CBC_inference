
/* From Diamond, 2000

neuronal transporters were configured by modifying existing transporter models
(Wadiche and Kavanaugh, 1998; Otis and Kavanaugh, 2000), such that glutamate
was removed from the simulation during transport and reverse transport did not
occur. The transporter Markov model comprised four states, connected by the
follow- ing rates (forward, backward; units are sec ?1 or M ?1 sec ?1): Tout ??
GTout (2 ? 10 7 ? [Glu], 300); GTout ?? GTin (500, 0); GTin ?? Tin (2000, 0);
and Tin ?? Tout (40, 0). In an effort to detect transporter-mediated currents
in pyramidal cells, Bergles and Jahr (1998) recorded synaptic responses with
SCN ? as the major internal anion to maximize the conductance through the
transporter-associated anion conductance (Wadiche et al., 1995). Those
experiments were simulated here by convolving the time course over which
transport events occurred during the simulation with an exponential wave- form
(peak, 0.427 fA; ? of 3 msec, i.e., 8 qe/cycle) to reflect the deactivation of
transporter currents in patches (Bergles and Jahr, 1998) and the en- hanced
current with SCN ? (Watzke et al., 2001). To simulate currents with an
impermanent anion, the exponential waveform was scaled down to reflect only the
net inward flux resulting from the electrogenic transport cycle (2 qe/cycle)
(Zerangue and Kavanaugh, 1996).
*/


/*

Rates: /M/s or /s

Tout <> GTout (2e7 * [Glu], 300); 
GTout <> GTin (500, 0); 
GTin <> Tin (2000, 0); and
Tin <> Tout (40, 0). 

F1 = 2e7*[Glu]
R1 = 300
F2 = 500
R2 = 0
F3 = 2000
R3 = 0
F4 = 40
R4 = 0

         F1
 Tout   --->   GTout 
         R1
   /|             /| 
F4  | R4        R2 | F2
    |/             |/   
         F3
  Tin  <---   GTin
         R3 
*/



/* 

See:
- - - - - - - - - - - - - - - -

Wadiche JI, Kavanaugh MP. (1998) Macroscopic and microscopic properties of a cloned glutamate transporter/chloride channel. J Neurosci. 1998 Oct 1;18(19):7650-61

Figure 8. Computer simulation. A, K inetic model of L-glutamate and D-aspartate
transport and anion conductance. Model parameters were obtained by least
squares fitting of data and fit an average pulse of L-glutamate of D-aspartate.
The microscopic rates were as follows: kon out = 6.8e6/M/sec, koff
out = 30.6/sec; k1 = 16.0/ sec; k-1 = 2.9/sec; kon in = 6.8e6/ M/
sec; koff in = 37.2/sec; k2 = 885/sec; k-2 = 200/se; a1 = 8094/sec; 
b1 = 100/sec; a2 = 1260/sec; b2 = 70 sec. D-Aspartate data were
fit with identical rates for agonist independent states and koff out = 7.6/sec; 
k1 = 7.3/sec; k-1 = 1.0/sec; koff(in) = 165/sec; a2 = 978/sec,
and b2 = 70/sec. DHK binding was assigned as 6.8e6/M/sec, whereas
DHK unbinding (k dhkout) = 97/sec. 

B, Probability of occupancy for each state in the kinetic scheme shown in A
during a 250 msec pulse of 10 mM L-glutamate. The top traces show the
nonconducting states: the unliganded states are represented by dashed lines
(Tout and Tin; bold), whereas the liganded states are represented by a solid
line (ToutGlu; bold and TinGlu). The bottom traces show the occupancy of the
anion conducting states. Note the different scale bars. C, Simulation of a 250
msec pulse of 10 mM L-glutamate or D-aspartate ( A). The channel?s steady-state
open probability was determined from nonstationary noise analysis (Fig. 4) and
the DHK-blocked currents (Fig.  6). The fraction of transporters in either
conducting state TGluopen or Topen are plotted as a f unction of time. D,
Concentration dependence of the time constant for activation of L-glutamate
currents for the kinetic scheme shown in A. The time constants for the
activation and deactivation were calculated by fitting the current records to a
single exponential.  The limiting slope for the activation rate equals
6.8e6/M/sec. Inset, L-Glutamate concentration dependence of the open
probability (1 uM, 10 uM, 100 uM, and 1 mM). The model?s apparent affinity at
steady state is 7 uM.

Table 2. Kinetic parameters
	                     L-Glutamate L-Glutamate model   D-Aspartate D-Aspartate model
--------------------------------------------------------------------------------
     tau activation   (msec) 0.96 +- 0.1 (19)  0.83     2.66 +- 0.2 (19)    1.5
     tau deactivation (msec) 22.8 +- 3.9 (9)   22.5    75.1 +- 10.5 (9)    72.4
     tau inactivation (msec) 14.1 +- 2.4 (18)  18.3    85.4 +- 15.7 (11)   94.7
    Peak/steady-state ratio  1.56 +- 0.1 (11)   1.5     1.05 +- 0.02 (17)   1.05
--------------------------------------------------------------------------------

 Kinetic parameters were determined at ?80 mV with a 10 mM application of L-glutamate or D-aspartate.



kon(out)  = 6.8e6 /M/s, 
koff(out) = 30.6/s
k1        = 16.0/s
k-1       = 2.9/s
kon(in)   = 6.8e6/M/s
koff(in)  = 37.2/s
k2        = 885/s
k-2       = 200/s;
a1        = 8094/s
b1        = 100/s
a2        = 1260/s
b2        = 70/s

D-Aspartate data were fit with identical rates for agonist independent states and 
koff(out) = 7.6/s
k1        = 7.3/s 
k-1       = 1.0/s
koff(in)  = 165/s
a2        = 978/s
b2        = 70/s

DHK binding = 6.8e6/M/s, 
DHK unbinding = (k dhkout) = 97/s. 


   Tout*         ToutGlu*
    /|              /|
     |               |
  b1 | a1         b2 | a2
     |               |
     |/  konout      |/
   Tout  ----->   ToutGlu
    /|   koffout    /|
     |               |
  k2 | k-2       k-1 | k1
     |               |
     |/   konin      |/
    Tin  ----->   TinGlu
         koffin



- - - - - - - - - - - - - - - - - - - - - - - - - -

Otis TS, Kavanaugh MP. (2000) Isolation of current components and partial reaction cycles in the glial glutamate transporter EAAT2. J Neurosci. 2000 Apr 15;20(8):2749-57.


Figure 8. A simple kinetic model that accounts for the two components
of transporter current. A, The model consists of three sets of states
connected by the following rate constants: k1 = 0 or 2000/sec; k-1 = 4/sec;
k2 = 480/sec * exp[(-VF )/(2 RT )]; k-2 = 400/sec * exp[(VF )/(2 RT )];
k3 = 10/sec * exp[(-VF )/(2 RT )]; and k-3 = 1/sec * exp[(VF )/(2 RT )]. V
of -90 mV for simulations in B and C and ?100 mV in D. B, Simulated
stoichiometric current (the sum of the net fluxes for the two voltage-
dependent transitions [conducting] 3 [nonconductingin] and [noncon-
ductingin] 3 [nonconductingout]) in response to a 100 msec duration jump
in the value of k1 from 0 to 2000/sec (to simulate the Glu pulse). C,
Simulated anion current (the occupancy of [conducting]) in response to
the Glu pulse. D, Anion current in response to the Glu pulse with the
effects of high [Na]in and [Glu]in simulated by increasing the rate constant
factor in k?2 from 400/sec to 6000/sec.



             nonconducting(out)

             /|                /|
         k3   |  k-3        k-1 |  k1
              |/                |/
                        k-2
 nonconducting(in)   <-----    conducting
                        k2

k1  = 0 or 2000/sec
k-1 = 4/sec;
k2  = 480/sec * exp[(-VF )/(2 RT )]
k-2 = 400/sec * exp[(VF )/(2 RT )]
k3  = 10/sec * exp[(-VF )/(2 RT )]
k-3 = 1/sec * exp[(VF )/(2 RT )]

of ?90 mV for simulations in B and C and ?100 mV in D. B, Simulated
stoichiometric current (the sum of the net fluxes for the two voltage-

- - - - - - - - - - - - - - - -


*/

