/* control.h */

int nrecrd;
int nstim;
int vidmode;			/* sets video output mode */
int silent;			/* =1 -> info=0, turns off info printout */
int info;			/* sets level of information in printout */
int disp;			/* sets display mode */

int prmap;
int stim_elem;
int fread_expr;			/* fread files with expressions (no spaces) */

/* default values for some constants */

double vk;
double vna;
double vcl;
double vglu;
double dnai;             	/* internal [Na] = 12 mM, used for Na vrev */
double dnao;             	/* external [Na] = 130 mM, used for Na vrev */
double dpkna;              	/* K perm / Na perm = .09, used for Na vrev */
double dpcana;              	/* Ca perm / Na perm = .09, used for Na vrev */
double dko;              	/* external [K] = 6.5 mM, used for K vrev */
double dki;              	/* internal [K] = 130 mM, used for K,Ca vrev */
double dpnak;              	/* Na perm / K perm = .02, used for K vrev */
double dpcak;              	/* Ca perm / K perm = .01, used for K vrev */
double dpnaca;              	/* Na perm / Ca perm = 1e-3, for Ca vrev */
double dpkca;              	/* K perm / Ca perm =.33e-3, for Ca vrev */
double dpcaampa;              	/* rel Ca perm for ampa synapse, for Ca flux */
double dpcacgmp;              	/* rel Ca perm for cGMP synapse, for Ca flux */
double dpcanmda;              	/* rel Ca perm for NMDA synapse, for Ca flux */
double dpcasyn2;              	/* rel Ca perm for syn2 synapse, for Ca flux */
double dclo;             	/* external [Cl] = 125 mM, used for Cl vrev */
double dcli;             	/* internal [Cl] = 10 mM, used for Cl vrev */
double deleccap;             	/* electrode capacitance */
double dcm;             	/* membrane capacitance */
double ddia;
double dlen;
double dri;			/* axial resistance */
double drg;			/* gap junction resistance */
double drm;			/* membrane resistance */
double drs;			/* electrode resistance */
double dsfa;			/* synaptic filters 1 */
double dfta;			/* synaptic time constant */
double dftah;			/* high pass synaptic time constant */
double dsfb;			/* synaptic filters 2 */
double dftb;			/* synaptic time constant */
double dscu;			/* unitary conductance of chans */
double dsvn;			/* number of vesicle release sites */
double dscavg;			/* synaptic ca vesicle release gain */
double dscaeg;			/* synaptic ca power for release gain */
double dvg;			/* vesicle release gain */
double dvsz;			/* size of vesicles released */
double dst;			/* synaptic threshold */
double dsc;			/* synaptic curve */
double dmaxsyn;			/* maximum cond for synapse */
double dskd;			/* half-max saturation gain */
double dshc;			/* synaptic postsyn Hill coeff */
double dstr;			/* transmitter conc for markov st receptors */
double dsmsgc;			/* transmitter conc for 2nd mesg receptors */
double dsms;			/* maximum sustained release rate */
double dsmrrp;			/* maximum readily releasible pool */
double dsrrg;			/* readily releasible pool gain mult */
double dsrrp;			/* readily releasible pool */
double dcakd;			/* half-max Ca saturation at chan */
double dcahc;			/* hill coeff for Ca binding at chan */
double dckd;			/* half-max cG saturation */
double dchc;			/* hill coeff for cGMP binding at chan */
double dsn;			/* synaptic neurotransmitter gain */
double dsg;			/* synaptic cycG gain */
double dsei;			/* synaptic input gain for expon */
double dmaxrod;			/* maximum cond for rod          */
double dmaxcon;			/* maximum cond for cone         */
double dnadens;			/* density for Na chan (S/cm2)   */
double dkdens;			/* density for K chan (S/cm2)    */
double dcadens;			/* density for Ca chan (S/cm2)   */
double dmaxna;			/* maximum cond for Na channel   */ 
double dnau;			/* unitary cond for Na channel   */       
double dku;			/* unitry cond for K channel     */
double dkau;			/* unitry cond for KA channel    */
double dkcabu;			/* unitry cond for Kca BK channel*/
double dkcasu;			/* unitry cond for Kca SK channel*/
double dclcasu;			/* unitry cond for Clca channel*/
double dkihu;			/* unitry cond for Ih channel    */
double dhcnu;			/* unitry cond for HCN channel    */
double dcalu;			/* unitry cond for Ca T channel  */
double dcatu;			/* unitry cond for Ca L channel  */
double dampau;			/* unitry cond for AMPA channel  */
double dcgmpu;			/* unitry cond for cGMP channel  */
double dgabau;			/* unitry cond for GABA channel  */
double dglyu;			/* unitry cond for GLY channel   */
double dpru;			/* unitry cond for photoreceptor channel   */
double dmaxampa;		/* maximum cond for AMPA channel */ 
double dmaxnmda;		/* maximum cond for NMDA channel */ 
double dmaxgaba;		/* maximum cond for GABA channel */ 
double dnaoffsm;		/* offset activation voltage for Na channel */ 
double dnaoffsh;		/* offset inactivation voltage for Na channel */ 
double dmaxk;			/* maximum cond for K channel */ 
double dmaxclca;		/* maximum cond for Clca channel */ 
double dclcavs;			/* ClCa chan voltage sens over 120 mV */ 
double dclcavsc;		/* ClCa chan voltage sens over 120 mV, core */ 
double dsclcac;			/* ClCa chan core Ca sens */ 
double dkoffsn;			/* offset activation   voltage for K channel */ 
double dkoffsh;			/* offset inactivation voltage for K channel */ 
double dnatauf;			/* tau for Na noise in 2-state model */
double dktauf;			/* tau for K noise in 2-state model */
double dcatauf;			/* tau for Ca noise in 2-state model */
double dnataum;			/* act. tau for Na channel */ 
double dnatauh;			/* inact. tau for Na channel */ 
double dktaum;			/* act. tau for K channel    */ 
double dktauh;			/* inact. tau for K channel  */ 
double dcataum;			/* act. tau for Ca channel  */ 
double dcatauh;			/* inact. tau for Ca channel */ 
double dkcataum;		/* act. tau for Kca channel    */ 
double dkcatauh;		/* inact. tau for Kca channel  */ 
double dsyntau;			/* tau for syn2 channel  */ 
double dsyntauf;		/* tau for syn2 channel noise */ 
double dnmdamg;			/* mg conc for NMDA channel */ 
double dqm;               	/* Q10 for am, bm */
double dqh;               	/* Q10 for ah, bh, normally 3 */
double dqn;               	/* Q10 for an, bn, normally 3 */
double dqna;               	/* Q10 for an, normally 3.2 */
double dqnb;               	/* Q10 for bn, normally 2.8 */
double dqca;               	/* Q10 for ac, bc (Ca chan) */
double dqd;               	/* Q10 for ad, bd (KA h) */
double dqsyn;             	/* Q10 for synaptic chans */
double dqcab;             	/* Q10 for Ca buffer */
double dqkca;             	/* Q10 for Kca chan */
double dqc;             		/* Q10 for channel conductance */
double dqcavmax;             	/* Q10 for Ca pump */
double dqri;             	/* Q10 for axial conductance */
double dqrm;             	/* Q10 for membrane conductance */
double dqdc;             	/* Q10 for diffusion constant */
double dqrec;             	/* Q10 for photoreceptor transduction shape */
double dqcrec;             	/* Q10 for photoreceptor transduction cond */
double dgjoff;			/* gap junction voltage offset */
double dgjtau;			/* activation time const for gap junction */
double dgjatpdec;		/* decrement for atp per time step */
double dgjnv;			/* non-volt-sens fraction of gj conductance */
double dcavoff;          	/* Q10 for voffset from cao */ 
double dcaspvrev;          	/* Ca surface potential for vrev(mult dcavoff)*/
double dcaphg;          	/* Ca rel conductance from pH */

                                /* Dan Emerson's CICR parameters: */
double betaIP3;          	/* init fractional rate of constant IP3 store release, default 0.31 */
double v1IP3;            	/* init constant IP3 store flux to 7.3e-6 MS-1 (Goldbeter 1990) */
double kfCICR;           	/* passive leak rate const from ryanodine store, 1/s (Goldbeter 1990) */
double vm2CICR;          	/* max rate of Ca pump into ryanodine store, 65e-6/ms (Goldbeter 1990) */
double vm3CICR;          	/* max rate of Ca pump from ryanodine store,500e-6/ms (Goldbeter 1990) */
double k1CICR;           	/* assoc thresh pump const, ryan store uptake,5e-6 M (Goldbeter 1990) */
double k2CICR;           	/* assoc thresh pump const, ryan store uptake,1e-6 M (Goldbeter 1990) */
double krCICR;           	/* assoc thresh pump const, ryan store release,2e-6 M (Goldbeter 1990) */
double kaCICR;           	/* assoc thresh pump const, ryan store release,0.9e-6 M (Goldbeter 1990) */
double nCICR;            	/* Hill coeff, coop pump binding ryan store uptake, 1 (Goldbeter 1990 */
double mCICR;            	/* Hill coeff, coop pump binding ryan store release, 1 (Goldbeter 1990 */
double pCICR;            	/* Hill coeff, coop pump binding ryan store release, 4 (Goldbeter 1990 */
double c1CICR;                  /* ratio of ER vol to cytoplasmic vol, default is 0.185 (Rinzel & Li 1994) */

double oIP3;                    /* Hill coeff, coop pump binding ip3 store uptake, 2 (Rinzel and Li) */
double a3IP3;                   /* IP3 store constant for IP3 Tau H, 2e5 1/MS (Rinzel and Li 1994) */
double mtypeIP3;                /* 0 => static (def = use equilibrated inf m value),(Rinzel and Li, 1994) */
double d1IP3;                   /* IP3 store const for Q value for h gate , 0.13e-6 M (Rinzel and Li, 1994) */
double d2IP3;                   /* IP3 store const for Q value for h gate, 1.049e-6 M (Rinzel and Li, 1994) */
double d3IP3;                   /* IP3 store const for Q value for h gate, 0.9434e-6 M (Rinzel and Li, 1994) */
double d4IP3;                   /* IP3 store const for m gate at infinity, 0.82e-6 M (Rinzel and Li, 1994) */
double a2IP3;                   /* IP3 store forward const for Tau M value, */
double b2IP3;                   /* IP3 store backward const for Tau M value */
double k3IP3;                   /* IP3 store kd for Ca uptake, 0.1 S-1 (Rinzel and Li, 1994) */
double v3IP3;                   /* IP3 store constant for Ca uptake, 0.9 S-1 (Rinzel and Li, 1994) */
double v2IP3;                   /* IP3 store rate for Ca release, 0.11 M S-1 (Rinzel and Li, 1994) */
double v4IP3;                   /* IP3 store rate for Ca release, 6 M S-1 (Rinzel and Li, 1994) */

double ddca;			/* Diffusion cons for internal calcium cm2/s */
double ddcao;			/* Diffusion const for external calcium cm2/s */
double dcoo;			/* external cobalt concentration (for cavoff) */
double dmgo;			/* external magnesium conc (for cavoff) */
double dcao;			/* external calcium concentration [Ca]o */
double dcai;			/* internal calcium concentration [Ca]i */
double dtcai;			/* threshold [Ca]i for pump */
double dcapkm;			/* 1/2 sat. conc. for Ca pump */
double dcavmax;			/* Vmax for Ca pump */
double dicafrac;		/* Fraction of Ica added to comp Itot */
double dcakex;			/* Rate for Na/Ca exchanger */
double diexchfrac;		/* Fraction of Icapump added to comp Itot */
double dcabnd;			/* Ratio of bound to free calcium */
double dcabt;			/* Ca buffer total in shells */
double dcabti;			/* Ca buffer total in first shell */
double dcabf;			/* forward rate to Ca.B */
double dcabr;			/* reverse rate to Ca.B */
double dmaxca;			/* maximum cond for Ca channel */ 
double dcaoffs;			/* Voltage activation offset for ca chan */
double dbd1;			/* default d1 for KCa BK channel    */ 
double dbd2;			/* default d2 for KCa BK channel    */ 
double dbk1;			/* default k1 for KCa BK channel    */ 
double dbk2;			/* default k2 for KCa BK channel    */ 
double dsd1;			/* default d1 for KCa SK channel    */ 
double dsd2;			/* default d2 for KCa SK channel    */ 
double dsk1;			/* default k1 for KCa SK channel    */ 
double dsk2;			/* default k2 for KCa SK channel    */ 
double dbasetc;			/* default base temp for chan kinetics, cond */ 
double dbasetsyn;		/* default base temp for synaptic kinetics */ 
double dbasetca;		/* default base temp for calcium pumps */ 
double dbasetdc;		/* default base temp for diffusion constant */ 
double dratehhm;		/* default base multiplier for HH m rate func */
double dratehhh;		/* default base multiplier for HH h rate func */
double dratehhna;        	/* default base multiplier for HH alphan func */
double dratehhnb;        	/* default base multiplier for HH betan func */
double dqeff;        		/* default quantum efficiency, R* per/photon */
double dsintinc;       		/* default min time incr for sine waves */
double dsintres;       		/* default time res for sine waves */

double complam;
double lamcrit;
double simtime;
double timinc;
double stiminc;
double tempcel;
double dpki;
double frt;
double f2rt;
double r2ft;
double ktemp;
double ktrb;
double djnoise;
double endexp;
double psplen;
double recrd;
double ploti;
double rate;
double stloc;
double crit;
double relax;
double relincr;
double synaptau;
double calcnernst;	/* =1 -> realc Nernst pots., not internal ions */
double plmax;
double plmin;
double setxmax;		/* overrides xmax for plotinit() */
double setxmin;		/* overrides xmin for plotinit() */
double setymax;		/* overrides ymax for plotinit() */
double setymin;		/* overrides ymin for plotinit() */
double unmaskthr;	/* threshold for unmasking light   */
double stimonh;		/* high threshold for transducer vclamp */
double stimonl;		/* low threshold for transducer vclamp */
double implfk;		/* implicit factor, for C-N -> 0.5 */
double ddiacao;		/* dia of external ca core */
double dcasarea;	/* area of int shell for Ca flux */
double srtimestep;	/* time step for reading stimuli, def 1e-4 */

int    dcashell;	/* Number of int Ca diffusion shells */
int    dcaoshell;	/* Number of ext Ca diffusion shells */
int    debug;
int    debugz;
int    poisfl;
int    plsep;
int    allowrev;
int    dashfl;
int    euler;
int    implicit;
int    stimelem;
int    vcolor;
int    lcolor;
int    rcolor;
int    cacolor;
int    nacolor;
int    sgcolor;
int    srcolor;
int    nozerochan;
int    oldsynapz;
int    use_ghki;
int    numstimchan;
int rseed;
int srseed;
int chrseed;
int dkrseed;
int phrseed;
extern int drseed;
double ncversion;
int makestim;
int prcomps;
