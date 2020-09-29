/* control.h */

extern int nrecrd;
extern int nstim;
extern int vidmode;		/* sets video output mode */
extern int silent;		/* =1 -> info=0, turns off info printout */
extern int info;		/* sets level of information in printout */
extern int disp;		/* sets display mode */
extern int prmap;
extern int stim_elem;
extern int fread_expr;		/* fread files with expressions (no spaces) */

/* default values for some constants */

extern double vk;
extern double vna;
extern double vcl;
extern double vglu;
extern double dnai;		/* internal [Na] = 12 mM, used for Na vrev */
extern double dnao;		/* external [Na] = 145 mM, used for Na vrev */
extern double dpkna;		/* K perm / Na perm = .09, used for Na vrev */
extern double dpcana;		/* Ca perm / Na perm = .09, used for Na vrev */
extern double dko;		/* internal [K] = 5 mM, used for K vrev */
extern double dki;		/* internal [K] = 130 mM, used for K,Ca vrev */
extern double dpnak;		/* Na perm / K perm = .02, used for K vrev */
extern double dpcak;		/* Ca perm / K perm = .01, used for K vrev */
extern double dpnaca;		/* Na permeability rel to Ca, for Ca vrev */
extern double dpkca;		/* K permeability rel to Ca, for Ca vrev */
extern double dpcaampa;		/* rel Ca perm for AMPA synapse, for Ca flux */
extern double dpcacgmp;		/* rel Ca perm for cGMP synapse, for Ca flux */
extern double dpcanmda;		/* rel Ca perm for NMDA synapse, for Ca flux */
extern double dpcasyn2;		/* rel Ca perm for syn2 synapse, for Ca flux */
extern double dclo;		/* external [Cl] = 125 mM, used for Cl vrev */
extern double dcli;		/* internal [Cl] = 10 mM, used for Cl vrev */
extern double deleccap;		/* electrode capacitance */
extern double dcm;		/* membrane capacitance */
extern double ddia;
extern double dlen;
extern double dri;		/* axial resistance */
extern double drg;		/* gap junction resistance */
extern double drm;		/* membrane resistance */
extern double drs;		/* electrode resistance */
extern double dsfa;		/* synaptic filters 1 */
extern double dfta;		/* synaptic time constant */
extern double dftah;		/* high pass synaptic time constant */
extern double dsfb;		/* synaptic filters 2 */
extern double dftb;		/* synaptic time constant */
extern double dscu;		/* unitary conductance of chans */
extern double dsvn;		/* number of vesicle release sites */
extern double dscavg;		/* synaptic ca vesicle release gain */
extern double dscaeg;		/* synaptic ca power for release gain */
extern double dvg;		/* vesicle release gain */
extern double dvsz;		/* size of vesicles released */
extern double dst;		/* synaptic threshold */
extern double dsc;		/* synaptic curve */
extern double dskd;		/* half-max saturation gain */
extern double dshc;		/* synaptic postsyn Hill coeff */
extern double dstr;		/* transmitter conc for markov st receptors */
extern double dsmsgc;		/* transmitter conc for 2nd mesg receptors */
extern double dsms;		/* maximum sustained release rate */
extern double dsmrrp;		/* maximum readily releasible pool */
extern double dsrrg;		/* readily releasible pool gain mult */
extern double dsrrp;		/* readily releasible pool */
extern double dckd;		/* half-max cG saturation */
extern double dchc;		/* hill coeff for cGMP binding at chan */
extern double dcakd;		/* half-max Ca saturation at chan */
extern double dcahc;		/* hill coeff for Ca binding at chan */
extern double dmaxsyn;		/* maximum cond for synapse      */
extern double dsn;		/* synaptic neurotransmitter gain */
extern double dsg;		/* synaptic cycG gain */
extern double dsei;		/* synaptic input gain for expon  */
extern double dmaxrod;		/* maximum cond for rod           */
extern double dmaxcon;		/* maximum cond for cone          */
extern double dnadens;		/* density for Na chan (S/cm2)    */
extern double dkdens;		/* density for K chan (S/cm2)     */
extern double dcadens;		/* density for Ca chan (S/cm2)    */
extern double dmaxna;		/* maximum cond for Na channel    */ 
extern double dnau;		/* unitary cond for Na channel    */ 
extern double dku;		/* unitary cond for K channel      */ 
extern double dkau;		/* unitary cond for KA channel     */ 
extern double dkcabu;		/* unitary cond for Kca BK channel */ 
extern double dkcasu;		/* unitary cond for Kca SK channel */ 
extern double dclcasu;		/* unitary cond for Clca channel */ 
extern double dkihu;		/* unitary cond for Ih channel     */ 
extern double dhcnu;		/* unitary cond for HCN channel     */ 
extern double dcalu;		/* unitary cond for Ca L channel   */ 
extern double dcatu;		/* unitary cond for Ca T channel   */ 
extern double dampau;		/* unitary cond for AMPA channel   */ 
extern double dcgmpu;		/* unitary cond for cGMP channel   */ 
extern double dgabau;		/* unitary cond for GABA channel   */ 
extern double dglyu;		/* unitary cond for GLY channel   */ 
extern double dpru;		/* unitary cond for photoreceptor channel   */ 
extern double dmaxampa;		/* maximum cond for AMPA channel  */ 
extern double dmaxnmda;		/* maximum cond for NMDA channel  */ 
extern double dmaxgaba;		/* maximum cond for GABA channel  */ 
extern double dnaoffsm;		/* offset activation   voltage for NA chan */
extern double dnaoffsh;		/* offset inactivation voltage for NA chan */
extern double dmaxk;		/* maximum cond for K channel     */ 
extern double dmaxclca;		/* maximum cond for ClCa channel     */ 
extern double dclcavs;		/* ClCa chan voltage sens over 120 mV */ 
extern double dclcavsc;		/* ClCa chan voltage sens over 120 mV, core */ 
extern double dsclcac;		/* ClCa chan core Ca sens */ 
extern double dkoffsn;		/* offset activation   voltage for K chan  */
extern double dkoffsh;		/* offset inactivation voltage for K chan  */
extern double dnatauf;		/* tau for Na noise in 2-state model */
extern double dktauf;		/* tau for K noise in 2-state model */
extern double dcatauf;		/* tau for Ca noise in 2-state model */
extern double dnataum;		/* act. tau for NA chan */
extern double dnatauh;		/* inact. tau for NA chan  */
extern double dktaum;		/* act. tau for K chan */
extern double dktauh;		/* inact. tau for K chan */
extern double dcataum;		/* act. tau for CA chan */
extern double dcatauh;		/* inact. tau for CA chan */
extern double dkcataum;		/* act. tau for Kca chan */
extern double dkcatauh;		/* inact. tau for Kca chan */
extern double dsyntau;		/* tau for syn2 channel */
extern double dsyntauf;		/* tau for syn2 channel noise */
extern double dnmdamg;		/* mg conc for NMDA chan */
extern double dqm;		/* Q10 for am, bm */
extern double dqh;		/* Q10 for ah, bh, normally 3 */
extern double dqn;		/* Q10 for an, bn, normally 3 */
extern double dqna;		/* Q10 for an, normally 3.2 */
extern double dqnb;		/* Q10 for bn, normally 2.8 */
extern double dqca;		/* Q10 for ac, bc (Ca chan) */
extern double dqd;		/* Q10 for ad, bd (KA h) */
extern double dqsyn;		/* Q10 for synaptic chans */
extern double dqcab;		/* Q10 for Ca buffer */
extern double dqkca;		/* Q10 for Kca chan */
extern double dqc;		/* Q10 for channel conductance */
extern double dqcavmax;		/* Q10 for Ca pump */
extern double dqri;		/* Q10 for axial conductance */
extern double dqrm;		/* Q10 for membrane conductance */
extern double dqdc;		/* Q10 for diffusion constant */
extern double dqrec;		/* Q10 for photoreceptor transduction shape */
extern double dqcrec;		/* Q10 for photoreceptor transduction cond */
extern double dgjoff;		/* gap junction voltage offset */
extern double dgjtau;		/* activation time const for gap junction */
extern double dgjatpdec;	/* decrement for atp per time step */
extern double dgjnv;		/* non-volt-sens fraction of gj conductance */
extern double dcavoff;		/* Q10 for voffset from cao */
extern double dcaspvrev;	/* Ca surface potential for vrev(mult dcavoff) */
extern double dcaphg;	        /* Ca rel conductance from pH (Barnes & Bui, 1991) */

				/* Dan Emerson's CICR parameters: */
extern double betaIP3;		/* init fractional rate of constant IP3 store release, default 0.31 */
extern double v1IP3;		/* init constant IP3 store flux to 7.3e-6 MS-1 (Goldbeter 1990) */
extern double kfCICR;		/* passive leak rate const from ryanodine store, 1/s (Goldbeter 1990) */
extern double vm2CICR;		/* max rate of Ca pump into ryanodine store, 65e-6/ms (Goldbeter 1990) */
extern double vm3CICR;		/* max rate of Ca pump from ryanodine store,500e-6/ms (Goldbeter 1990) */
extern double k1CICR;		/* assoc thresh pump const, ryan store uptake,5e-6 M (Goldbeter 1990) */
extern double k2CICR;		/* assoc thresh pump const, ryan store uptake,1e-6 M (Goldbeter 1990) */
extern double krCICR;		/* assoc thresh pump const, ryan store release,2e-6 M (Goldbeter 1990) */
extern double kaCICR;		/* assoc thresh pump const, ryan store release,0.9e-6 M (Goldbeter 1990) */
extern double nCICR;		/* Hill coeff, coop pump binding ryan store uptake, 1 (Goldbeter 1990 */
extern double mCICR;		/* Hill coeff, coop pump binding ryan store release, 1 (Goldbeter 1990 */
extern double pCICR;		/* Hill coeff, coop pump binding ryan store release, 4 (Goldbeter 1990 */
extern double c1CICR;           /* ratio of ER vol to cytoplasmic vol, default is 0.185 (Rinzel & Li 1994) */

extern double oIP3;             /* Hill coeff, coop pump binding ip3 store uptake, 2 (Rinzel and Li) */
extern double a3IP3;            /* IP3 store constant for IP3 Tau H, 2e5 1/MS (Rinzel and Li 1994) */
extern double mtypeIP3;         /* 0 => static (def = use equilibrated inf m value),(Rinzel and Li, 1994) */
extern double d1IP3;            /* IP3 store const for Q value for h gate , 0.13e-6 M (Rinzel and Li, 1994) */
extern double d2IP3;            /* IP3 store const for Q value for h gate, 1.049e-6 M (Rinzel and Li, 1994) */
extern double d3IP3;            /* IP3 store const for Q value for h gate, 0.9434e-6 M (Rinzel and Li, 1994) */
extern double d4IP3;            /* IP3 store const for m gate at infinity, 0.82e-6 M (Rinzel and Li, 1994) */
extern double a2IP3;            /* IP3 store forward const for Tau M value, */
extern double b2IP3;            /* IP3 store backward const for Tau M value */
extern double k3IP3;            /* IP3 store kd for Ca uptake, 0.1 S-1 (Rinzel and Li, 1994) */
extern double v3IP3;            /* IP3 store constant for Ca uptake, 0.9 S-1 (Rinzel and Li, 1994) */
extern double v2IP3;            /* IP3 store rate for Ca release, 0.11 S-1 (Rinzel and Li, 1994) */ 
extern double v4IP3;            /* IP3 store rate for Ca release, 6 M S-1 (Rinzel and Li, 1994) */ 


extern double ddca;		/* Diffusion const for internal calcium cm2/s */
extern double ddcao;		/* Diffusion const for external calcium cm2/s */
extern double dcoo;		/* external cobalt concentration (for cavoff) */
extern double dmgo;		/* external magnesium conc (for cavoff) */
extern double dcao;		/* external calcium concentration [Ca]o */
extern double dcai;		/* internal calcium concentration [Ca]i */
extern double dtcai;		/* threshold [Ca]i for pump */
extern double dcapkm;		/* 1/2 sat Ca conc for pump */
extern double dcavmax;		/* Vmax for Ca pump (ma/cm2) */
extern double dicafrac;		/* Fraction of Ica added to comp Itot */
extern double dcakex;		/* Rate for Na/Ca exchanger */
extern double diexchfrac;	/* Fraction of Icapump added to comp Itot */
extern double dcabnd;		/* Ratio of bound to free calcium */
extern double dcabt;		/* Ca buffer total in shells */
extern double dcabti;		/* Ca buffer total in first shell */
extern double dcabf;		/* forward rate to  Ca.B */
extern double dcabr;		/* reverse rate to Ca.B */
extern double dcaoffs;		/* Voltage offset for ca chan */
extern double dmaxca;		/* maximum cond for Ca channel    */ 
extern double dbd1;		/* default d1 for KCa BK channel    */ 
extern double dbd2;		/* default d2 for KCa BK channel    */ 
extern double dbk1;		/* default k1 for KCa BK channel    */ 
extern double dbk2;		/* default k2 for KCa BK channel    */ 
extern double dsd1;		/* default d1 for KCa SK channel    */ 
extern double dsd2;		/* default d2 for KCa SK channel    */ 
extern double dsk1;		/* default k1 for KCa SK channel    */ 
extern double dsk2;		/* default k2 for KCa SK channel    */ 
extern double dbasetc;		/* default base temp for chan kinetics, cond */
extern double dbasetsyn;	/* default base temp for synaptic kinet */
extern double dbasetca;		/* default base temp for ca pumps */
extern double dbasetdc;		/* default base temp for diffusion constant */
extern double dratehhm;		/* default base multiplier for HH m rate func */
extern double dratehhh;		/* default base multiplier for HH h rate func */
extern double dratehhna;	/* default base multiplier for HH alphan func */
extern double dratehhnb;	/* default base multiplier for HH betan func */
extern double dqeff;		/* default quantum efficiency R* per photon */
extern double dsintinc;		/* default min time incr for sine waves */
extern double dsintres;		/* default time res for sine waves */

extern double complam;
extern double lamcrit;
extern double simtime;
extern double timinc;
extern double stiminc;
extern double tempcel;
extern double ktemp;
extern double ktrb;
extern double dpki;
extern double frt;
extern double f2rt;
extern double r2ft;
extern double djnoise;
extern double endexp;
extern double psplen;
extern double recrd;
extern double ploti;
extern double rate;
extern double stloc;
extern double crit;
extern double relax;
extern double relincr;
extern double synaptau;
extern double calcnernst;	/* =1 -> calc Nernst pots., not internal ions */
extern double plmax;
extern double plmin;
extern double setxmax;		/* overrides xmax for plotinit() */
extern double setxmin;		/* overrides xmin for plotinit() */
extern double setymax;		/* overrides ymax for plotinit() */
extern double setymin;		/* overrides ymin for plotinit() */
extern double unmaskthr;	/* threshold for unmasking light   */
extern double stimonh;		/* high threshold for transducer vclamp */
extern double stimonl;		/* low threshold for transducer vclamp */
extern double implfk;		/* implicit factor, for C-N -> 0.5 */
extern double ddiacao;		/* dia of external ca core */
extern double dcasarea;		/* int area of shell for Ca flux */
extern double lm_usertol;	/* user-set epsilon for lmfit() */
extern double srtimestep;	/* time step for reading stimuli, def 1e-4 */

extern int    dcashell;		/* Number of int Ca diffusion shells */
extern int    dcaoshell;	/* Number of ext Ca diffusion shells */
extern int    debug;
extern int    debugz;
extern int    poisfl;
extern int    plsep;
extern int    allowrev;
extern int    dashfl;
extern int    euler;
extern int    implicit;
extern int    stimelem;
extern int    vcolor;
extern int    lcolor;
extern int    rcolor;
extern int    nacolor;
extern int    sgcolor;
extern int    srcolor;
extern int    cacolor;
extern int    nozerochan;
extern int    oldsynapz;
extern int    use_ghki;
extern int    numstimchan;	/* number of stimulus channels */
extern int rseed;
extern int srseed;
extern int phrseed;
extern int chrseed;
extern int dkrseed;
extern int drseed;
extern double ncversion;
extern int makestim;
extern int prcomps;
