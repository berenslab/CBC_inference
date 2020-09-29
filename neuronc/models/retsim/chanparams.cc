
/* chanparams.cc */
/* Variable defs for chanparams file */

/*--------------------------------------------------------------------*/

/* defs for chanparams */

int _columns = 0;

int _NA0_offm = C_NA0*NCH+CHOFFM;     /* NaV1.2, no noise, type 0 */
int _NA0_offh = C_NA0*NCH+CHOFFH;
int _NA0_taua = C_NA0*NCH+CHTAUA;
int _NA0_taub = C_NA0*NCH+CHTAUB;
int _NA0_tauc = C_NA0*NCH+CHTAUC;
int _NA0_taud = C_NA0*NCH+CHTAUD;

int _NA1_offm = C_NA1*NCH+CHOFFM;     /* NaV1.2, simple noise */
int _NA1_offh = C_NA1*NCH+CHOFFH;
int _NA1_taua = C_NA1*NCH+CHTAUA;
int _NA1_taub = C_NA1*NCH+CHTAUB;
int _NA1_tauc = C_NA1*NCH+CHTAUC;
int _NA1_taud = C_NA1*NCH+CHTAUD;

int _NA2_offm = C_NA2*NCH+CHOFFM;     /* NaV1,2 same as type 1, but markov */
int _NA2_offh = C_NA2*NCH+CHOFFH;
int _NA2_taua = C_NA2*NCH+CHTAUA;
int _NA2_taub = C_NA2*NCH+CHTAUB;
int _NA2_tauc = C_NA2*NCH+CHTAUC;
int _NA2_taud = C_NA2*NCH+CHTAUD;

int _NA3_offm = C_NA3*NCH+CHOFFM;     /* Na: type 3 */
int _NA3_offh = C_NA3*NCH+CHOFFH;
int _NA3_taua = C_NA3*NCH+CHTAUA;
int _NA3_taub = C_NA3*NCH+CHTAUB;
int _NA3_tauc = C_NA3*NCH+CHTAUC;
int _NA3_taud = C_NA3*NCH+CHTAUD;

int _NA4_offm = C_NA4*NCH+CHOFFM;     /* Na: type 4 */
int _NA4_offh = C_NA4*NCH+CHOFFH;
int _NA4_taua = C_NA4*NCH+CHTAUA;
int _NA4_taub = C_NA4*NCH+CHTAUB;
int _NA4_tauc = C_NA4*NCH+CHTAUC;
int _NA4_taud = C_NA4*NCH+CHTAUD;

int _NA5_offm = C_NA5*NCH+CHOFFM;     /* slowly inactivating Na: type 5 */
int _NA5_offh = C_NA5*NCH+CHOFFH;
int _NA5_taua = C_NA5*NCH+CHTAUA;
int _NA5_taub = C_NA5*NCH+CHTAUB;
int _NA5_tauc = C_NA5*NCH+CHTAUC;
int _NA5_taud = C_NA5*NCH+CHTAUD;

int _NA6_offm = C_NA6*NCH+CHOFFM;     /* persistent Na: type 6 */
int _NA6_offh = C_NA6*NCH+CHOFFH;
int _NA6_taua = C_NA6*NCH+CHTAUA;
int _NA6_taub = C_NA6*NCH+CHTAUB;
int _NA6_tauc = C_NA6*NCH+CHTAUC;
int _NA6_taud = C_NA6*NCH+CHTAUD;

int _NA8_offm = C_NA8*NCH+CHOFFM;     /* Nav1.8, TTX-resistant Na: type 8 */
int _NA8_offh = C_NA8*NCH+CHOFFH;
int _NA8_taua = C_NA8*NCH+CHTAUA;
int _NA8_taub = C_NA8*NCH+CHTAUB;
int _NA8_tauc = C_NA8*NCH+CHTAUC;
int _NA8_taud = C_NA8*NCH+CHTAUD;



int _K0_offm = C_K0*NCH+CHOFFM;      /* type 0 , no noise */
int _K0_offh = C_K0*NCH+CHOFFH;
int _K0_taua = C_K0*NCH+CHTAUA;
int _K0_taub = C_K0*NCH+CHTAUB;
int _K0_tauc = C_K0*NCH+CHTAUC;
int _K0_taud = C_K0*NCH+CHTAUD;

int _K1_offm = C_K1*NCH+CHOFFM;      /* type 0, noise: type 1 */
int _K1_offh = C_K1*NCH+CHOFFH;
int _K1_taua = C_K1*NCH+CHTAUA;
int _K1_taub = C_K1*NCH+CHTAUB;
int _K1_tauc = C_K1*NCH+CHTAUC;
int _K1_taud = C_K1*NCH+CHTAUD;

int _K2_offm = C_K2*NCH+CHOFFM;      /* type 2, noise: type 1 */
int _K2_offh = C_K2*NCH+CHOFFH;
int _K2_taua = C_K2*NCH+CHTAUA;
int _K2_taub = C_K2*NCH+CHTAUB;
int _K2_tauc = C_K2*NCH+CHTAUC;
int _K2_taud = C_K2*NCH+CHTAUD;

int _K3_offm = C_K3*NCH+CHOFFM;      /* KA, type 3 */
int _K3_offh = C_K3*NCH+CHOFFH;
int _K3_taua = C_K3*NCH+CHTAUA;
int _K3_taub = C_K3*NCH+CHTAUB;
int _K3_tauc = C_K3*NCH+CHTAUC;
int _K3_taud = C_K3*NCH+CHTAUD;

int _K4_offm = C_K4*NCH+CHOFFM;      /* Ih, type 4 */
int _K4_offh = C_K4*NCH+CHOFFH;
int _K4_taua = C_K4*NCH+CHTAUA;
int _K4_taub = C_K4*NCH+CHTAUB;
int _K4_tauc = C_K4*NCH+CHTAUC;
int _K4_taud = C_K4*NCH+CHTAUD;

// joesterle begin

int _HCN1_offm = C_HCN1*NCH+CHOFFM;
int _HCN1_offh = C_HCN1*NCH+CHOFFH;
int _HCN1_taua = C_HCN1*NCH+CHTAUA;
int _HCN1_taub = C_HCN1*NCH+CHTAUB;
int _HCN1_tauc = C_HCN1*NCH+CHTAUC;
int _HCN1_taud = C_HCN1*NCH+CHTAUD;

int _HCN2_offm = C_HCN2*NCH+CHOFFM;
int _HCN2_offh = C_HCN2*NCH+CHOFFH;
int _HCN2_taua = C_HCN2*NCH+CHTAUA;
int _HCN2_taub = C_HCN2*NCH+CHTAUB;
int _HCN2_tauc = C_HCN2*NCH+CHTAUC;
int _HCN2_taud = C_HCN2*NCH+CHTAUD;

int _HCN3_offm = C_HCN3*NCH+CHOFFM;
int _HCN3_offh = C_HCN3*NCH+CHOFFH;
int _HCN3_taua = C_HCN3*NCH+CHTAUA;
int _HCN3_taub = C_HCN3*NCH+CHTAUB;
int _HCN3_tauc = C_HCN3*NCH+CHTAUC;
int _HCN3_taud = C_HCN3*NCH+CHTAUD;

int _HCN4_offm = C_HCN4*NCH+CHOFFM;
int _HCN4_offh = C_HCN4*NCH+CHOFFH;
int _HCN4_taua = C_HCN4*NCH+CHTAUA;
int _HCN4_taub = C_HCN4*NCH+CHTAUB;
int _HCN4_tauc = C_HCN4*NCH+CHTAUC;
int _HCN4_taud = C_HCN4*NCH+CHTAUD;
// joesterle end

int _K5_offm = C_K5*NCH+CHOFFM;      /* Kir, Hz cell, type 5 */
int _K5_offh = C_K5*NCH+CHOFFH;
int _K5_taua = C_K5*NCH+CHTAUA;
int _K5_taub = C_K5*NCH+CHTAUB;
int _K5_tauc = C_K5*NCH+CHTAUC;
int _K5_taud = C_K5*NCH+CHTAUD;

int _K6_offm = C_K6*NCH+CHOFFM;      /* Kv3, SBAC, type 6 */
int _K6_offh = C_K6*NCH+CHOFFH;
int _K6_taua = C_K6*NCH+CHTAUA;
int _K6_taub = C_K6*NCH+CHTAUB;
int _K6_tauc = C_K6*NCH+CHTAUC;
int _K6_taud = C_K6*NCH+CHTAUD;

int _K7_offm = C_K7*NCH+CHOFFM;      /* Kv3b, SBAC, type 7 */
int _K7_offh = C_K7*NCH+CHOFFH;
int _K7_taua = C_K7*NCH+CHTAUA;
int _K7_taub = C_K7*NCH+CHTAUB;
int _K7_tauc = C_K7*NCH+CHTAUC;
int _K7_taud = C_K7*NCH+CHTAUD;

int _KCA0_offm = C_KCA0*NCH+CHOFFM;     /* BK, HH-type */
int _KCA0_offh = C_KCA0*NCH+CHOFFH;
int _KCA0_taua = C_KCA0*NCH+CHTAUA;
int _KCA0_taub = C_KCA0*NCH+CHTAUB;
int _KCA0_tauc = C_KCA0*NCH+CHTAUC;
int _KCA0_taud = C_KCA0*NCH+CHTAUD;

int _KCA1_offm = C_KCA1*NCH+CHOFFM;     /* SK, no voltage sensitivity, Markov */
int _KCA1_offh = C_KCA1*NCH+CHOFFH;
int _KCA1_taua = C_KCA1*NCH+CHTAUA;
int _KCA1_taub = C_KCA1*NCH+CHTAUB;
int _KCA1_tauc = C_KCA1*NCH+CHTAUC;
int _KCA1_taud = C_KCA1*NCH+CHTAUD;

int _KCA2_offm = C_KCA2*NCH+CHOFFM;     /* BK, Markov */
int _KCA2_offh = C_KCA2*NCH+CHOFFH;
int _KCA2_taua = C_KCA2*NCH+CHTAUA;
int _KCA2_taub = C_KCA2*NCH+CHTAUB;
int _KCA2_tauc = C_KCA2*NCH+CHTAUC;
int _KCA2_taud = C_KCA2*NCH+CHTAUD;

int _KCA3_offm = C_KCA3*NCH+CHOFFM;     /* BK, fast, Markov */
int _KCA3_offh = C_KCA3*NCH+CHOFFH;
int _KCA3_taua = C_KCA3*NCH+CHTAUA;
int _KCA3_taub = C_KCA3*NCH+CHTAUB;
int _KCA3_tauc = C_KCA3*NCH+CHTAUC;
int _KCA3_taud = C_KCA3*NCH+CHTAUD;

int _KCA4_offm = C_KCA4*NCH+CHOFFM;     /* SK 1, moderately slow, Markov */
int _KCA4_offh = C_KCA4*NCH+CHOFFH;
int _KCA4_taua = C_KCA4*NCH+CHTAUA;
int _KCA4_taub = C_KCA4*NCH+CHTAUB;
int _KCA4_tauc = C_KCA4*NCH+CHTAUC;
int _KCA4_taud = C_KCA4*NCH+CHTAUD;

int _KCA5_offm = C_KCA5*NCH+CHOFFM;     /* SK 2,slow, Markov */
int _KCA5_offh = C_KCA5*NCH+CHOFFH;
int _KCA5_taua = C_KCA5*NCH+CHTAUA;
int _KCA5_taub = C_KCA5*NCH+CHTAUB;
int _KCA5_tauc = C_KCA5*NCH+CHTAUC;
int _KCA5_taud = C_KCA5*NCH+CHTAUD;

int _KCA6_offm = C_KCA6*NCH+CHOFFM;     /* SK 2,slow, Markov */
int _KCA6_offh = C_KCA6*NCH+CHOFFH;
int _KCA6_taua = C_KCA6*NCH+CHTAUA;
int _KCA6_taub = C_KCA6*NCH+CHTAUB;
int _KCA6_tauc = C_KCA6*NCH+CHTAUC;
int _KCA6_taud = C_KCA6*NCH+CHTAUD;



int _CLCA1_offm = C_CLCA1*NCH+CHOFFM;     /* Calcium-actived chloride current */
int _CLCA1_offh = C_CLCA1*NCH+CHOFFH;
int _CLCA1_taua = C_CLCA1*NCH+CHTAUA;
int _CLCA1_taub = C_CLCA1*NCH+CHTAUB;
int _CLCA1_tauc = C_CLCA1*NCH+CHTAUC;
int _CLCA1_taud = C_CLCA1*NCH+CHTAUD;

int _CLCA2_offm = C_CLCA2*NCH+CHOFFM;     /* Calcium-actived chloride current, core Ca */
int _CLCA2_offh = C_CLCA2*NCH+CHOFFH;
int _CLCA2_taua = C_CLCA2*NCH+CHTAUA;
int _CLCA2_taub = C_CLCA2*NCH+CHTAUB;
int _CLCA2_tauc = C_CLCA2*NCH+CHTAUC;
int _CLCA2_taud = C_CLCA2*NCH+CHTAUD;



int _CA0_offm = C_CA0*NCH+CHOFFM;     /* L-type Ca */
int _CA0_offh = C_CA0*NCH+CHOFFH;
int _CA0_taua = C_CA0*NCH+CHTAUA;
int _CA0_taub = C_CA0*NCH+CHTAUB;
int _CA0_tauc = C_CA0*NCH+CHTAUC;
int _CA0_taud = C_CA0*NCH+CHTAUD;

int _CA1_offm = C_CA1*NCH+CHOFFM;     /* L-type Ca, Markov */
int _CA1_offh = C_CA1*NCH+CHOFFH;
int _CA1_taua = C_CA1*NCH+CHTAUA;
int _CA1_taub = C_CA1*NCH+CHTAUB;
int _CA1_tauc = C_CA1*NCH+CHTAUC;
int _CA1_taud = C_CA1*NCH+CHTAUD;

int _CA2_offm = C_CA2*NCH+CHOFFM;     /* T-type Ca */
int _CA2_offh = C_CA2*NCH+CHOFFH;
int _CA2_taua = C_CA2*NCH+CHTAUA;
int _CA2_taub = C_CA2*NCH+CHTAUB;
int _CA2_tauc = C_CA2*NCH+CHTAUC;
int _CA2_taud = C_CA2*NCH+CHTAUD;

int _CA3_offm = C_CA3*NCH+CHOFFM;     /* T-type Ca, Markov */
int _CA3_offh = C_CA3*NCH+CHOFFH;
int _CA3_taua = C_CA3*NCH+CHTAUA;
int _CA3_taub = C_CA3*NCH+CHTAUB;
int _CA3_tauc = C_CA3*NCH+CHTAUC;
int _CA3_taud = C_CA3*NCH+CHTAUD;

int _CA4_offm = C_CA4*NCH+CHOFFM;     /* T-type Ca */
int _CA4_offh = C_CA4*NCH+CHOFFH;
int _CA4_taua = C_CA4*NCH+CHTAUA;
int _CA4_taub = C_CA4*NCH+CHTAUB;
int _CA4_tauc = C_CA4*NCH+CHTAUC;
int _CA4_taud = C_CA4*NCH+CHTAUD;

int _CA5_offm = C_CA5*NCH+CHOFFM;     /* T-type Ca, Markov */
int _CA5_offh = C_CA5*NCH+CHOFFH;
int _CA5_taua = C_CA5*NCH+CHTAUA;
int _CA5_taub = C_CA5*NCH+CHTAUB;
int _CA5_tauc = C_CA5*NCH+CHTAUC;
int _CA5_taud = C_CA5*NCH+CHTAUD;

int _CA6_offm = C_CA6*NCH+CHOFFM;     /* T-type Ca, low thresh */
int _CA6_offh = C_CA6*NCH+CHOFFH;
int _CA6_taua = C_CA6*NCH+CHTAUA;
int _CA6_taub = C_CA6*NCH+CHTAUB;
int _CA6_tauc = C_CA6*NCH+CHTAUC;
int _CA6_taud = C_CA6*NCH+CHTAUD;

int _CA7_offm = C_CA7*NCH+CHOFFM;     /* T-type Ca, low thresh */
int _CA7_offh = C_CA7*NCH+CHOFFH;
int _CA7_taua = C_CA7*NCH+CHTAUA;
int _CA7_taub = C_CA7*NCH+CHTAUB;
int _CA7_tauc = C_CA7*NCH+CHTAUC;
int _CA7_taud = C_CA7*NCH+CHTAUD;


int _NMDA1_offm = C_NMDA1*NCH+CHOFFM;     /* NMDA simple */
int _NMDA1_offh = C_NMDA1*NCH+CHOFFH;
int _NMDA1_taua = C_NMDA1*NCH+CHTAUA;
int _NMDA1_taub = C_NMDA1*NCH+CHTAUB;
int _NMDA1_tauc = C_NMDA1*NCH+CHTAUC;
int _NMDA1_taud = C_NMDA1*NCH+CHTAUD;

int _NMDA2_offm = C_NMDA2*NCH+CHOFFM;     /* NMDA complex */
int _NMDA2_offh = C_NMDA2*NCH+CHOFFH;
int _NMDA2_taua = C_NMDA2*NCH+CHTAUA;
int _NMDA2_taub = C_NMDA2*NCH+CHTAUB;
int _NMDA2_tauc = C_NMDA2*NCH+CHTAUC;
int _NMDA2_taud = C_NMDA2*NCH+CHTAUD;

// joesterle

int _MGLUR_offm = C_MGLUR*NCH+CHOFFM;
int _MGLUR_offh = C_MGLUR*NCH+CHOFFH;
int _MGLUR_taua = C_MGLUR*NCH+CHTAUA;
int _MGLUR_taub = C_MGLUR*NCH+CHTAUB;
int _MGLUR_tauc = C_MGLUR*NCH+CHTAUC;
int _MGLUR_taud = C_MGLUR*NCH+CHTAUD;

// joesterle

int _AMPA1_offm = C_AMPA1*NCH+CHOFFM;     /* AMPA type 1 */
int _AMPA1_offh = C_AMPA1*NCH+CHOFFH;
int _AMPA1_taua = C_AMPA1*NCH+CHTAUA;
int _AMPA1_taub = C_AMPA1*NCH+CHTAUB;
int _AMPA1_tauc = C_AMPA1*NCH+CHTAUC;
int _AMPA1_taud = C_AMPA1*NCH+CHTAUD;

int _AMPA2_offm = C_AMPA2*NCH+CHOFFM;     /* AMPA type 2 */
int _AMPA2_offh = C_AMPA2*NCH+CHOFFH;
int _AMPA2_taua = C_AMPA2*NCH+CHTAUA;
int _AMPA2_taub = C_AMPA2*NCH+CHTAUB;
int _AMPA2_tauc = C_AMPA2*NCH+CHTAUC;
int _AMPA2_taud = C_AMPA2*NCH+CHTAUD;

int _AMPA3_offm = C_AMPA3*NCH+CHOFFM;     /* AMPA type 3 */
int _AMPA3_offh = C_AMPA3*NCH+CHOFFH;
int _AMPA3_taua = C_AMPA3*NCH+CHTAUA;
int _AMPA3_taub = C_AMPA3*NCH+CHTAUB;
int _AMPA3_tauc = C_AMPA3*NCH+CHTAUC;
int _AMPA3_taud = C_AMPA3*NCH+CHTAUD;

int _AMPA4_offm = C_AMPA4*NCH+CHOFFM;     /* AMPA type 4 */
int _AMPA4_offh = C_AMPA4*NCH+CHOFFH;
int _AMPA4_taua = C_AMPA4*NCH+CHTAUA;
int _AMPA4_taub = C_AMPA4*NCH+CHTAUB;
int _AMPA4_tauc = C_AMPA4*NCH+CHTAUC;
int _AMPA4_taud = C_AMPA4*NCH+CHTAUD;

int _AMPA5_offm = C_AMPA5*NCH+CHOFFM;     /* AMPA type 5 */
int _AMPA5_offh = C_AMPA5*NCH+CHOFFH;
int _AMPA5_taua = C_AMPA5*NCH+CHTAUA;
int _AMPA5_taub = C_AMPA5*NCH+CHTAUB;
int _AMPA5_tauc = C_AMPA5*NCH+CHTAUC;
int _AMPA5_taud = C_AMPA5*NCH+CHTAUD;


int _GABA1_offm = C_GABA1*NCH+CHOFFM;     /* GABA type 1 */
int _GABA1_offh = C_GABA1*NCH+CHOFFH;
int _GABA1_taua = C_GABA1*NCH+CHTAUA;
int _GABA1_taub = C_GABA1*NCH+CHTAUB;
int _GABA1_tauc = C_GABA1*NCH+CHTAUC;
int _GABA1_taud = C_GABA1*NCH+CHTAUD;

int _GABA2_offm = C_GABA2*NCH+CHOFFM;     /* GABA type 2 */
int _GABA2_offh = C_GABA2*NCH+CHOFFH;
int _GABA2_taua = C_GABA2*NCH+CHTAUA;
int _GABA2_taub = C_GABA2*NCH+CHTAUB;
int _GABA2_tauc = C_GABA2*NCH+CHTAUC;
int _GABA2_taud = C_GABA2*NCH+CHTAUD;

int _GABA3_offm = C_GABA3*NCH+CHOFFM;     /* GABA type 3 */
int _GABA3_offh = C_GABA3*NCH+CHOFFH;
int _GABA3_taua = C_GABA3*NCH+CHTAUA;
int _GABA3_taub = C_GABA3*NCH+CHTAUB;
int _GABA3_tauc = C_GABA3*NCH+CHTAUC;
int _GABA3_taud = C_GABA3*NCH+CHTAUD;

int _GABA4_offm = C_GABA4*NCH+CHOFFM;     /* GABA type 4 */
int _GABA4_offh = C_GABA4*NCH+CHOFFH;
int _GABA4_taua = C_GABA4*NCH+CHTAUA;
int _GABA4_taub = C_GABA4*NCH+CHTAUB;
int _GABA4_tauc = C_GABA4*NCH+CHTAUC;
int _GABA4_taud = C_GABA4*NCH+CHTAUD;


int _GLY_offm = C_GLY*NCH+CHOFFM;     /* glycine type 1 */
int _GLY_offh = C_GLY*NCH+CHOFFH;
int _GLY_taua = C_GLY*NCH+CHTAUA;
int _GLY_taub = C_GLY*NCH+CHTAUB;
int _GLY_tauc = C_GLY*NCH+CHTAUC;
int _GLY_taud = C_GLY*NCH+CHTAUD;

int _CGMP1_offm = C_CGMP1*NCH+CHOFFM;     /* GGMP type 1 */
int _CGMP1_offh = C_CGMP1*NCH+CHOFFH;
int _CGMP1_taua = C_CGMP1*NCH+CHTAUA;
int _CGMP1_taub = C_CGMP1*NCH+CHTAUB;
int _CGMP1_tauc = C_CGMP1*NCH+CHTAUC;
int _CGMP1_taud = C_CGMP1*NCH+CHTAUD;

int _CGMP2_offm = C_CGMP2*NCH+CHOFFM;     /* GGMP type 2 */
int _CGMP2_offh = C_CGMP2*NCH+CHOFFH;
int _CGMP2_taua = C_CGMP2*NCH+CHTAUA;
int _CGMP2_taub = C_CGMP2*NCH+CHTAUB;
int _CGMP2_tauc = C_CGMP2*NCH+CHTAUC;
int _CGMP2_taud = C_CGMP2*NCH+CHTAUD;

int _CGMP3_offm = C_CGMP3*NCH+CHOFFM;     /* GGMP type 3 */
int _CGMP3_offh = C_CGMP3*NCH+CHOFFH;
int _CGMP3_taua = C_CGMP3*NCH+CHTAUA;
int _CGMP3_taub = C_CGMP3*NCH+CHTAUB;
int _CGMP3_tauc = C_CGMP3*NCH+CHTAUC;
int _CGMP3_taud = C_CGMP3*NCH+CHTAUD;

int _CGMP4_offm = C_CGMP4*NCH+CHOFFM;     /* GGMP type 4 */
int _CGMP4_offh = C_CGMP4*NCH+CHOFFH;
int _CGMP4_taua = C_CGMP4*NCH+CHTAUA;
int _CGMP4_taub = C_CGMP4*NCH+CHTAUB;
int _CGMP4_tauc = C_CGMP4*NCH+CHTAUC;
int _CGMP4_taud = C_CGMP4*NCH+CHTAUD;

int _CGMP5_offm = C_CGMP5*NCH+CHOFFM;     /* GGMP type 5 */
int _CGMP5_offh = C_CGMP5*NCH+CHOFFH;
int _CGMP5_taua = C_CGMP5*NCH+CHTAUA;
int _CGMP5_taub = C_CGMP5*NCH+CHTAUB;
int _CGMP5_tauc = C_CGMP5*NCH+CHTAUC;
int _CGMP5_taud = C_CGMP5*NCH+CHTAUD;

int _CGMP6_offm = C_CGMP6*NCH+CHOFFM;     /* GGMP type 6 */
int _CGMP6_offh = C_CGMP6*NCH+CHOFFH;
int _CGMP6_taua = C_CGMP6*NCH+CHTAUA;
int _CGMP6_taub = C_CGMP6*NCH+CHTAUB;
int _CGMP6_tauc = C_CGMP6*NCH+CHTAUC;
int _CGMP6_taud = C_CGMP6*NCH+CHTAUD;

void setchanptrs(void) 

/* set ptrs to chanparam row names */

{
   setptrn ("columns", &_columns);

   setptrn ("NA0.offm", &_NA0_offm);
   setptrn ("NA0.offh", &_NA0_offh);
   setptrn ("NA0.taua", &_NA0_taua);
   setptrn ("NA0.taub", &_NA0_taub);
   setptrn ("NA0.tauc", &_NA0_tauc);
   setptrn ("NA0.taud", &_NA0_taud);

   setptrn ("NA1.offm", &_NA1_offm);
   setptrn ("NA1.offh", &_NA1_offh);
   setptrn ("NA1.taua", &_NA1_taua);
   setptrn ("NA1.taub", &_NA1_taub);
   setptrn ("NA1.tauc", &_NA1_tauc);
   setptrn ("NA1.taud", &_NA1_taud);

   setptrn ("NA2.offm", &_NA2_offm);
   setptrn ("NA2.offh", &_NA2_offh);
   setptrn ("NA2.taua", &_NA2_taua);
   setptrn ("NA2.taub", &_NA2_taub);
   setptrn ("NA2.tauc", &_NA2_tauc);
   setptrn ("NA2.taud", &_NA2_taud);

   setptrn ("NA3.offm", &_NA3_offm);
   setptrn ("NA3.offh", &_NA3_offh);
   setptrn ("NA3.taua", &_NA3_taua);
   setptrn ("NA3.taub", &_NA3_taub);
   setptrn ("NA3.tauc", &_NA3_tauc);
   setptrn ("NA3.taud", &_NA3_taud);

   setptrn ("NA4.offm", &_NA4_offm);
   setptrn ("NA4.offh", &_NA4_offh);
   setptrn ("NA4.taua", &_NA4_taua);
   setptrn ("NA4.taub", &_NA4_taub);
   setptrn ("NA4.tauc", &_NA4_tauc);
   setptrn ("NA4.taud", &_NA4_taud);

   setptrn ("NA5.offm", &_NA5_offm);
   setptrn ("NA5.offh", &_NA5_offh);
   setptrn ("NA5.taua", &_NA5_taua);
   setptrn ("NA5.taub", &_NA5_taub);
   setptrn ("NA5.tauc", &_NA5_tauc);
   setptrn ("NA5.taud", &_NA5_taud);

   setptrn ("NA6.offm", &_NA6_offm);
   setptrn ("NA6.offh", &_NA6_offh);
   setptrn ("NA6.taua", &_NA6_taua);
   setptrn ("NA6.taub", &_NA6_taub);
   setptrn ("NA6.tauc", &_NA6_tauc);
   setptrn ("NA6.taud", &_NA6_taud);

   setptrn ("NA8.offm", &_NA8_offm);
   setptrn ("NA8.offh", &_NA8_offh);
   setptrn ("NA8.taua", &_NA8_taua);
   setptrn ("NA8.taub", &_NA8_taub);
   setptrn ("NA8.tauc", &_NA8_tauc);
   setptrn ("NA8.taud", &_NA8_taud);


   setptrn ("K0.offm", &_K0_offm);
   setptrn ("K0.offh", &_K0_offh);
   setptrn ("K0.taua", &_K0_taua);
   setptrn ("K0.taub", &_K0_taub);
   setptrn ("K0.tauc", &_K0_tauc);
   setptrn ("K0.taud", &_K0_taud);

   setptrn ("K1.offm", &_K1_offm);
   setptrn ("K1.offh", &_K1_offh);
   setptrn ("K1.taua", &_K1_taua);
   setptrn ("K1.taub", &_K1_taub);
   setptrn ("K1.tauc", &_K1_tauc);
   setptrn ("K1.taud", &_K1_taud);

   setptrn ("K2.offm", &_K2_offm);
   setptrn ("K2.offh", &_K2_offh);
   setptrn ("K2.taua", &_K2_taua);
   setptrn ("K2.taub", &_K2_taub);
   setptrn ("K2.tauc", &_K2_tauc);
   setptrn ("K2.taud", &_K2_taud);

   setptrn ("K3.offm", &_K3_offm);
   setptrn ("K3.offh", &_K3_offh);
   setptrn ("K3.taua", &_K3_taua);
   setptrn ("K3.taub", &_K3_taub);
   setptrn ("K3.tauc", &_K3_tauc);
   setptrn ("K3.taud", &_K3_taud);

   setptrn ("K4.offm", &_K4_offm);
   setptrn ("K4.offh", &_K4_offh);
   setptrn ("K4.taua", &_K4_taua);
   setptrn ("K4.taub", &_K4_taub);
   setptrn ("K4.tauc", &_K4_tauc);
   setptrn ("K4.taud", &_K4_taud);
   
   // joesterle begin
   setptrn ("HCN1.offm", &_HCN1_offm);
   setptrn ("HCN1.offh", &_HCN1_offh);
   setptrn ("HCN1.taua", &_HCN1_taua);
   setptrn ("HCN1.taub", &_HCN1_taub);
   setptrn ("HCN1.tauc", &_HCN1_tauc);
   setptrn ("HCN1.taud", &_HCN1_taud);
   
   setptrn ("HCN2.offm", &_HCN2_offm);
   setptrn ("HCN2.offh", &_HCN2_offh);
   setptrn ("HCN2.taua", &_HCN2_taua);
   setptrn ("HCN2.taub", &_HCN2_taub);
   setptrn ("HCN2.tauc", &_HCN2_tauc);
   setptrn ("HCN2.taud", &_HCN2_taud);
   
   setptrn ("HCN3.offm", &_HCN3_offm);
   setptrn ("HCN3.offh", &_HCN3_offh);
   setptrn ("HCN3.taua", &_HCN3_taua);
   setptrn ("HCN3.taub", &_HCN3_taub);
   setptrn ("HCN3.tauc", &_HCN3_tauc);
   setptrn ("HCN3.taud", &_HCN3_taud);
   
   setptrn ("HCN4.offm", &_HCN4_offm);
   setptrn ("HCN4.offh", &_HCN4_offh);
   setptrn ("HCN4.taua", &_HCN4_taua);
   setptrn ("HCN4.taub", &_HCN4_taub);
   setptrn ("HCN4.tauc", &_HCN4_tauc);
   setptrn ("HCN4.taud", &_HCN4_taud);
   // joesterle end
   
   setptrn ("K5.offm", &_K5_offm);
   setptrn ("K5.offh", &_K5_offh);
   setptrn ("K5.taua", &_K5_taua);
   setptrn ("K5.taub", &_K5_taub);
   setptrn ("K5.tauc", &_K5_tauc);
   setptrn ("K5.taud", &_K5_taud);

   setptrn ("K6.offm", &_K6_offm);
   setptrn ("K6.offh", &_K6_offh);
   setptrn ("K6.taua", &_K6_taua);
   setptrn ("K6.taub", &_K6_taub);
   setptrn ("K6.tauc", &_K6_tauc);
   setptrn ("K6.taud", &_K6_taud);

   setptrn ("K7.offm", &_K7_offm);
   setptrn ("K7.offh", &_K7_offh);
   setptrn ("K7.taua", &_K7_taua);
   setptrn ("K7.taub", &_K7_taub);
   setptrn ("K7.tauc", &_K7_tauc);
   setptrn ("K7.taud", &_K7_taud);



   setptrn ("KCA0.offm", &_KCA0_offm);
   setptrn ("KCA0.offh", &_KCA0_offh);
   setptrn ("KCA0.taua", &_KCA0_taua);
   setptrn ("KCA0.taub", &_KCA0_taub);
   setptrn ("KCA0.tauc", &_KCA0_tauc);
   setptrn ("KCA0.taud", &_KCA0_taud);

   setptrn ("KCA1.offm", &_KCA1_offm);
   setptrn ("KCA1.offh", &_KCA1_offh);
   setptrn ("KCA1.taua", &_KCA1_taua);
   setptrn ("KCA1.taub", &_KCA1_taub);
   setptrn ("KCA1.tauc", &_KCA1_tauc);
   setptrn ("KCA1.taud", &_KCA1_taud);

   setptrn ("KCA2.offm", &_KCA2_offm);
   setptrn ("KCA2.offh", &_KCA2_offh);
   setptrn ("KCA2.taua", &_KCA2_taua);
   setptrn ("KCA2.taub", &_KCA2_taub);
   setptrn ("KCA2.tauc", &_KCA2_tauc);
   setptrn ("KCA2.taud", &_KCA2_taud);

   setptrn ("KCA3.offm", &_KCA3_offm);
   setptrn ("KCA3.offh", &_KCA3_offh);
   setptrn ("KCA3.taua", &_KCA3_taua);
   setptrn ("KCA3.taub", &_KCA3_taub);
   setptrn ("KCA3.tauc", &_KCA3_tauc);
   setptrn ("KCA3.taud", &_KCA3_taud);

   setptrn ("KCA4.offm", &_KCA4_offm);
   setptrn ("KCA4.offh", &_KCA4_offh);
   setptrn ("KCA4.taua", &_KCA4_taua);
   setptrn ("KCA4.taub", &_KCA4_taub);
   setptrn ("KCA4.tauc", &_KCA4_tauc);
   setptrn ("KCA4.taud", &_KCA4_taud);

   setptrn ("KCA5.offm", &_KCA5_offm);
   setptrn ("KCA5.offh", &_KCA5_offh);
   setptrn ("KCA5.taua", &_KCA5_taua);
   setptrn ("KCA5.taub", &_KCA5_taub);
   setptrn ("KCA5.tauc", &_KCA5_tauc);
   setptrn ("KCA5.taud", &_KCA5_taud);

   setptrn ("KCA6.offm", &_KCA6_offm);
   setptrn ("KCA6.offh", &_KCA6_offh);
   setptrn ("KCA6.taua", &_KCA6_taua);
   setptrn ("KCA6.taub", &_KCA6_taub);
   setptrn ("KCA6.tauc", &_KCA6_tauc);
   setptrn ("KCA6.taud", &_KCA6_taud);



   setptrn ("CLCA1.offm", &_CLCA1_offm);
   setptrn ("CLCA1.offh", &_CLCA1_offh);
   setptrn ("CLCA1.taua", &_CLCA1_taua);
   setptrn ("CLCA1.taub", &_CLCA1_taub);
   setptrn ("CLCA1.tauc", &_CLCA1_tauc);
   setptrn ("CLCA1.taud", &_CLCA1_taud);

   setptrn ("CLCA2.offm", &_CLCA2_offm);
   setptrn ("CLCA2.offh", &_CLCA2_offh);
   setptrn ("CLCA2.taua", &_CLCA2_taua);
   setptrn ("CLCA2.taub", &_CLCA2_taub);
   setptrn ("CLCA2.tauc", &_CLCA2_tauc);
   setptrn ("CLCA2.taud", &_CLCA2_taud);



   setptrn ("CA0.offm", &_CA0_offm);
   setptrn ("CA0.offh", &_CA0_offh);
   setptrn ("CA0.taua", &_CA0_taua);
   setptrn ("CA0.taub", &_CA0_taub);
   setptrn ("CA0.tauc", &_CA0_tauc);
   setptrn ("CA0.taud", &_CA0_taud);

   setptrn ("CA1.offm", &_CA1_offm);
   setptrn ("CA1.offh", &_CA1_offh);
   setptrn ("CA1.taua", &_CA1_taua);
   setptrn ("CA1.taub", &_CA1_taub);
   setptrn ("CA1.tauc", &_CA1_tauc);
   setptrn ("CA1.taud", &_CA1_taud);

   setptrn ("CA2.offm", &_CA2_offm);
   setptrn ("CA2.offh", &_CA2_offh);
   setptrn ("CA2.taua", &_CA2_taua);
   setptrn ("CA2.taub", &_CA2_taub);
   setptrn ("CA2.tauc", &_CA2_tauc);
   setptrn ("CA2.taud", &_CA2_taud);

   setptrn ("CA3.offm", &_CA3_offm);
   setptrn ("CA3.offh", &_CA3_offh);
   setptrn ("CA3.taua", &_CA3_taua);
   setptrn ("CA3.taub", &_CA3_taub);
   setptrn ("CA3.tauc", &_CA3_tauc);
   setptrn ("CA3.taud", &_CA3_taud);

   setptrn ("CA4.offm", &_CA4_offm);
   setptrn ("CA4.offh", &_CA4_offh);
   setptrn ("CA4.taua", &_CA4_taua);
   setptrn ("CA4.taub", &_CA4_taub);
   setptrn ("CA4.tauc", &_CA4_tauc);
   setptrn ("CA4.taud", &_CA4_taud);

   setptrn ("CA5.offm", &_CA5_offm);
   setptrn ("CA5.offh", &_CA5_offh);
   setptrn ("CA5.taua", &_CA5_taua);
   setptrn ("CA5.taub", &_CA5_taub);
   setptrn ("CA5.tauc", &_CA5_tauc);
   setptrn ("CA5.taud", &_CA5_taud);

   setptrn ("CA6.offm", &_CA6_offm);
   setptrn ("CA6.offh", &_CA6_offh);
   setptrn ("CA6.taua", &_CA6_taua);
   setptrn ("CA6.taub", &_CA6_taub);
   setptrn ("CA6.tauc", &_CA6_tauc);
   setptrn ("CA6.taud", &_CA6_taud);

   setptrn ("CA7.offm", &_CA7_offm);
   setptrn ("CA7.offh", &_CA7_offh);
   setptrn ("CA7.taua", &_CA7_taua);
   setptrn ("CA7.taub", &_CA7_taub);
   setptrn ("CA7.tauc", &_CA7_tauc);
   setptrn ("CA7.taud", &_CA7_taud);

   // joesterle
   setptrn ("MGLUR.offm", &_MGLUR_offm);
   setptrn ("MGLUR.offh", &_MGLUR_offh);
   setptrn ("MGLUR.taua", &_MGLUR_taua);
   setptrn ("MGLUR.taub", &_MGLUR_taub);
   setptrn ("MGLUR.tauc", &_MGLUR_tauc);
   setptrn ("MGLUR.taud", &_MGLUR_taud);
   // joesterle

   setptrn ("NMDA1.offm", &_NMDA1_offm);
   setptrn ("NMDA1.offh", &_NMDA1_offh);
   setptrn ("NMDA1.taua", &_NMDA1_taua);
   setptrn ("NMDA1.taub", &_NMDA1_taub);
   setptrn ("NMDA1.tauc", &_NMDA1_tauc);
   setptrn ("NMDA1.taud", &_NMDA1_taud);

   setptrn ("NMDA2.offm", &_NMDA2_offm);
   setptrn ("NMDA2.offh", &_NMDA2_offh);
   setptrn ("NMDA2.taua", &_NMDA2_taua);
   setptrn ("NMDA2.taub", &_NMDA2_taub);
   setptrn ("NMDA2.tauc", &_NMDA2_tauc);
   setptrn ("NMDA2.taud", &_NMDA2_taud);

   setptrn ("AMPA1.offm", &_AMPA1_offm);
   setptrn ("AMPA1.offh", &_AMPA1_offh);
   setptrn ("AMPA1.taua", &_AMPA1_taua);
   setptrn ("AMPA1.taub", &_AMPA1_taub);
   setptrn ("AMPA1.tauc", &_AMPA1_tauc);
   setptrn ("AMPA1.taud", &_AMPA1_taud);

   setptrn ("AMPA2.offm", &_AMPA2_offm);
   setptrn ("AMPA2.offh", &_AMPA2_offh);
   setptrn ("AMPA2.taua", &_AMPA2_taua);
   setptrn ("AMPA2.taub", &_AMPA2_taub);
   setptrn ("AMPA2.tauc", &_AMPA2_tauc);
   setptrn ("AMPA2.taud", &_AMPA2_taud);

   setptrn ("AMPA3.offm", &_AMPA3_offm);
   setptrn ("AMPA3.offh", &_AMPA3_offh);
   setptrn ("AMPA3.taua", &_AMPA3_taua);
   setptrn ("AMPA3.taub", &_AMPA3_taub);
   setptrn ("AMPA3.tauc", &_AMPA3_tauc);
   setptrn ("AMPA3.taud", &_AMPA3_taud);

   setptrn ("AMPA4.offm", &_AMPA4_offm);
   setptrn ("AMPA4.offh", &_AMPA4_offh);
   setptrn ("AMPA4.taua", &_AMPA4_taua);
   setptrn ("AMPA4.taub", &_AMPA4_taub);
   setptrn ("AMPA4.tauc", &_AMPA4_tauc);
   setptrn ("AMPA4.taud", &_AMPA4_taud);

   setptrn ("AMPA5.offm", &_AMPA5_offm);
   setptrn ("AMPA5.offh", &_AMPA5_offh);
   setptrn ("AMPA5.taua", &_AMPA5_taua);
   setptrn ("AMPA5.taub", &_AMPA5_taub);
   setptrn ("AMPA5.tauc", &_AMPA5_tauc);
   setptrn ("AMPA5.taud", &_AMPA5_taud);



   setptrn ("GABA1.offm", &_GABA1_offm);
   setptrn ("GABA1.offh", &_GABA1_offh);
   setptrn ("GABA1.taua", &_GABA1_taua);
   setptrn ("GABA1.taub", &_GABA1_taub);
   setptrn ("GABA1.tauc", &_GABA1_tauc);
   setptrn ("GABA1.taud", &_GABA1_taud);

   setptrn ("GABA2.offm", &_GABA2_offm);
   setptrn ("GABA2.offh", &_GABA2_offh);
   setptrn ("GABA2.taua", &_GABA2_taua);
   setptrn ("GABA2.taub", &_GABA2_taub);
   setptrn ("GABA2.tauc", &_GABA2_tauc);
   setptrn ("GABA2.taud", &_GABA2_taud);

   setptrn ("GABA3.offm", &_GABA3_offm);
   setptrn ("GABA3.offh", &_GABA3_offh);
   setptrn ("GABA3.taua", &_GABA3_taua);
   setptrn ("GABA3.taub", &_GABA3_taub);
   setptrn ("GABA3.tauc", &_GABA3_tauc);
   setptrn ("GABA3.taud", &_GABA3_taud);

   setptrn ("GABA4.offm", &_GABA4_offm);
   setptrn ("GABA4.offh", &_GABA4_offh);
   setptrn ("GABA4.taua", &_GABA4_taua);
   setptrn ("GABA4.taub", &_GABA4_taub);
   setptrn ("GABA4.tauc", &_GABA4_tauc);
   setptrn ("GABA4.taud", &_GABA4_taud);



   setptrn ("GLU.offm", &_GLY_offm);
   setptrn ("GLU.offh", &_GLY_offh);
   setptrn ("GLU.taua", &_GLY_taua);
   setptrn ("GLU.taub", &_GLY_taub);
   setptrn ("GLU.tauc", &_GLY_tauc);
   setptrn ("GLU.taud", &_GLY_taud);



   setptrn ("CGMP1.offm", &_CGMP1_offm);
   setptrn ("CGMP1.offh", &_CGMP1_offh);
   setptrn ("CGMP1.taua", &_CGMP1_taua);
   setptrn ("CGMP1.taub", &_CGMP1_taub);
   setptrn ("CGMP1.tauc", &_CGMP1_tauc);
   setptrn ("CGMP1.taud", &_CGMP1_taud);

   setptrn ("CGMP2.offm", &_CGMP2_offm);
   setptrn ("CGMP2.offh", &_CGMP2_offh);
   setptrn ("CGMP2.taua", &_CGMP2_taua);
   setptrn ("CGMP2.taub", &_CGMP2_taub);
   setptrn ("CGMP2.tauc", &_CGMP2_tauc);
   setptrn ("CGMP2.taud", &_CGMP2_taud);

   setptrn ("CGMP3.offm", &_CGMP3_offm);
   setptrn ("CGMP3.offh", &_CGMP3_offh);
   setptrn ("CGMP3.taua", &_CGMP3_taua);
   setptrn ("CGMP3.taub", &_CGMP3_taub);
   setptrn ("CGMP3.tauc", &_CGMP3_tauc);
   setptrn ("CGMP3.taud", &_CGMP3_taud);

   setptrn ("CGMP4.offm", &_CGMP4_offm);
   setptrn ("CGMP4.offh", &_CGMP4_offh);
   setptrn ("CGMP4.taua", &_CGMP4_taua);
   setptrn ("CGMP4.taub", &_CGMP4_taub);
   setptrn ("CGMP4.tauc", &_CGMP4_tauc);
   setptrn ("CGMP4.taud", &_CGMP4_taud);

   setptrn ("CGMP5.offm", &_CGMP5_offm);
   setptrn ("CGMP5.offh", &_CGMP5_offh);
   setptrn ("CGMP5.taua", &_CGMP5_taua);
   setptrn ("CGMP5.taub", &_CGMP5_taub);
   setptrn ("CGMP5.tauc", &_CGMP5_tauc);
   setptrn ("CGMP5.taud", &_CGMP5_taud);

   setptrn ("CGMP6.offm", &_CGMP6_offm);
   setptrn ("CGMP6.offh", &_CGMP6_offh);
   setptrn ("CGMP6.taua", &_CGMP6_taua);
   setptrn ("CGMP6.taub", &_CGMP6_taub);
   setptrn ("CGMP6.tauc", &_CGMP6_tauc);
   setptrn ("CGMP6.taud", &_CGMP6_taud);

}

