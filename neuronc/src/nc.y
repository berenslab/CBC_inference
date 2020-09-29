/* nc.y */
%{
#include <stdio.h>
#include "nc.h"
#define code2(c1,c2)	code(c1); code(c2)
#define code3(c1,c2,c3)	code(c1); code(c2); code(c3)
#define code4(c1,c2,c3,c4) code(c1); code(c2); code(c3); code(c4)
extern int indef, inlocal, savlocal, argcount, arrcount, formal, starg;
%}
%union { Symbol	*sym;	/* symbol table pointer */
	Inst	*inst;	/* machine instructions */
	int	narg;	/* number of arguments */
}
%token	<sym>	NUMBER STRING PRINT PRINTF SPRINTF FPRINTF NUMSTR STRNUM
%token	<sym>	SCANF FSCANF SSCANF FOPEN FCLOSE FGETS FGETC FPUTC GETFLDS
%token	<sym>	CONST VAR BLTIN UNDEF WHILE IF ELSE DEBUGF
%token	<sym>	FUNCTION PROCEDURE RETURN EXIT FUNC PROC 
%token  <sym>   READ FREAD FREADS FWRITE UNLINK
%token	<sym>   FOR BREAK CONTINUE EDIT INCLUDE COMNDLINE SYSTEM 
%token	<sym>   FOREACH WITHIN2D WITHIN3D OF
%token	<sym>	INCROP DECROP PFIELD ARG LOCAL GLOBAL LOCALVAR TEMP
%token	<sym>	DIM DIMS ARRAY ARRAYVAL ARRAYCONST
%token	<sym>	ARGARRVAL LOCALARR LOCALARRVAL SIZEOF
%token	<sym>	FFT IFFT PFFT ACOV LMFIT LMFIT2D
%token  <sym>   SFILE PLOT V VM I IM L S MAX MIN CHAR CCHAR VPEN
%token  <sym>   LITCHAR 
%token  <sym>   DISPLAY GRAPH INIT PEN RESTART
%token	<sym>	X Y Z XLOC YLOC ZLOC
%token	<sym>	FA0 FA1 FA2 FA3 FA4 FA8 FA9 FH0 FH1 FH2 FH3 FH4 
%token	<sym>	FB0 FB1 FB2 FB3 FB4 FC0 FC1 FC2 FC3 FC4 FC9
%token	<sym>	G M H G0 G1 G2 G3
%token	<sym>	GMOVE GDRAW GRMOVE GRDRAW GPEN GFRAME GSIZ GROT GCROT GCWID
%token	<sym>	GDASH GTEXT GLABEL GORIG GCIRC GRECT RMOVE GWINDOW GVPEN GPURGE
%token	<sym>	ELEMENT EXCEPT RANGE MATCHING SIZE WINDOW DSCALE COLOR HIDE COMPS
%token	<sym>	XROT YROT ZROT CALIBLIN PLNAME PLNUM PLSIZE PLVAL PLARR CMAP NEWPAGE 
%token	<sym>	CABLE SPHERE GJ PNX SYNAPSE LOAD RESISTOR DIODE ELECTRODE AXIALRES
%token	<sym>	MAXCOND CENTER WAVEL XENON SUN TUNGSTEN PIGM PATHL ATTF
%token	<sym>	ROD CONE PHOTREC ITRANSDUCER VTRANSDUCER LINIT CHR 
%token	<sym>	RRPOOL RRPOOLG MRRPOOL MAXSRATE
%token	<sym>	DYAD SDYAD SCURVE NUMDYAD SPOST NUMSPRE SENS
%token	<sym>	RSD INITRAND
%token	<sym>	STIM BAR SPOT SINE GABOR SINEANN WINDMILL SIMAGE WIDTH LOC DUR
%token	<sym>	CHECKERBOARD RECT	
%token	<sym>	TFREQ ORIENT SPHASE CONTRAST SSCALE DRIFT XENV YENV SQ
%token	<sym>	STIMCHAN MASK SECTOR
%token	<sym>	START STOP INTEN BACKGR BLUR SCATTER FILT SAVE RESTORE
%token	<sym>	LENGTH DIA DIA2 RADIUS VCLAMP CCLAMP GNDBATT BATT 
%token	<sym>	GNDCAP CAP NODE NNODE LABEL
%token	<sym>	CHAN HH NA K KCa ClCa MG
%token  <sym>	CA CAO CAI TCAI CACOMP IP IE IPE BTOT BTOTI 
%token  <sym>	CAPUMP VMAX KM CAEXCH KEX CASHELL CAOSHELL CAPERM CBOUND CABUF CABUFB
%token	<sym>	CICR CAS CAS2 VM2 VM3 P KA KF KR C1
%token	<sym>	IP3 IP3I BIP3 VIP3 V2 V3 V4 B2 OI A3 /* D1 D2*/ D3 D4 A2 K3 MTYPE MIP3 HIP3
%token	<sym>	TYPE NTYPE STYPE REGION DENSITY NDENSITY NSTATE
%token	<sym>	UNIT RM RI RG CM RS
%token	<sym>	OPAMP BUF NBUF PUFF JNOISE CONNECT AT TO 
%token	<sym>	EEDIST E3DIST E2DIST EZDIST EFRAC N3DIST N2DIST NZDIST
%token  <sym>	ONLY NUMCONN NUMSYN EXIST ELNUM
%token	<sym>	NODE1A NODE1B NODE1C NODE1D NODE2A NODE2B NODE2C NODE2D
%token	<sym>	CCHNOISE PHOTNOISE CHANNOISE DARKNOISE VESNOISE N VSIZE VCOV NUNIT
%token  <sym>   GLU AMPA KAINATE NMDA CNQX GABA BIC PTX GLY CHRC STRY CAMP CGMP PH ATP
%token  <sym>   LPFILT LP HP MFILT BESSFILT
%token	<sym>	GAUSNN NND NNSTD GINFO REG COV REFR 
%token	<sym>	NOTINIT VARNUM VARSTR VARCHR 
%token	<sym>	LINEAR EXPON NFILT1 NFILT1H NFILT2 NFILT3 
%token	<sym>	RESP MESGIN MESGOUT SYN2 LOOPG
%token	<sym>	TIMEC1 TIMEC1H TIMEC2 TIMEC3 TFALL2 TFALL3 
%token	<sym>	OPEN CLOSE THRESH VVREV VREST DELAY HGAIN CGAIN COFF KD HCOF CKD CHC
%token	<sym>	TRCONC MESGCONC TAU TAUM TAUH TAUN TAUA TAUB TAUC TAUD TAUE TAUF
%token	<sym>	D1 D2 K1 K2 AREA CAKD CAHC
%token	<sym>	BGAIN BOFFSET VGAIN EGAIN GMAX GNV REV
%token	<sym>	MODIFY ENAME ELABL PLACE CPLAM OFFSET OFFSETM OFFSETH PUT
%token	<sym>	RUN STEP ERASE MODEL ELIMIT
%type	<inst>	cexpr commaexpr setexpr expr stmt asgn prlist stmtlist stmtls
%type	<inst>	cond while if begin end defn filename comndline
%type	<inst>  for forexpr var foreach elemtype elemtypee within
%type	<inst>	elemlist geometry buf nbuf
%type	<inst>	stimulus stimtype stiminten stimfile stimwave
%type	<inst>  graphlist graphtype maxmin
%type	<inst>	cableparam gjtype sphereparam elecparam membrtype connect 
%type	<inst>	chantype noisetype jnoisetype stimparm nfilt
%type	<inst>	synaptype synapresp synapmesng synapsens receptype transducer
%type	<inst>	plotype plotypee plotypeb plotypef plotyped plotypeh
%type	<inst>	gprim displist disptype rotatype 
%type	<inst>	cachan cacomp cacomparm catype caexch capump cabuf cicr ip3
%type	<inst>	psca pscachan
%type	<inst>	limitlist limitype rinitlist chanfield
%type	<inst>	elemfield nodefield elemnum plotlist gnnlist gnnterm 
%type	<inst>	grltyp grphlst plltyp
%type	<sym>	procname listerm undefarrname arrayname
%type	<sym>	synparm memparm fftparm
%type	<sym>	chanparm receparm caparm timectype gjparam secondmsg ntrans
%type	<sym>	caplot ggplot localarrdef
%type	<narg>	arglist parglist formarg localvar vararg prflist sprflist 
%type	<narg>	scnflist scnfargs timecspec localdef prfargs
%type	<narg>	dimlist initlist nodenum nodenumex snode plotypec plotypeg
%right	'='
%right	ADDEQ SUBEQ MULEQ DIVEQ ANDEQ OREQ
%left	OR AND XOR BITAND BITOR 
%left	GT GE LT LE EQ NE
%left	'+' '-'
%left	'*' '/' '%'
%left	UNARYMINUS UNARYPLUS NOT
%right	'^'
%%


list:	  /* nothing */
	| list listerm		{ code(STOPC); return 1; }
/*	| list defn listerm	{ code(STOPC); return 1; } */
	| list asgn listerm	{ code2(expop,STOPC); return 1; }
	| list stmt listerm	{ code(STOPC); return 1; }
	| list expr listerm	{ code2(print, STOPC); return 1; }
	| list error listerm	{ yyerrok; }
	;
listerm: ';'			{}
	;
elemlist: connect		{}
	| elemlist membrtype	{}
	| elemlist geometry	{}
	| elemlist ELABL expr	{ code2(xmod,(Inst)$2); }
	| elemlist ENAME var	{ code2(xmod,(Inst)$2); }
	;
geometry: CABLE cableparam	{ $$=$2; }
	| SPHERE sphereparam	{ $$=$2; } 
	| synaptype		{ $$=$1; }
	| CHAN chantype		{ $$=$2; }
	| CHAN cachan		{ $$=$2; }
	| cacomp		{ $$=$1; }
	| receptype		{ $$=$1; }
	| ELECTRODE elecparam 	{ $$=$2; }
	| GJ gjtype		{ $$=$2; code2(xgj,(Inst)$1);} 
	| PNX gjtype		{ $$=$2; code2(xgj,(Inst)$1);} 
	| LOAD expr		{ $$=$2; code(rload); }
	| RESISTOR expr		{ $$=$2; code(xresistor); }
	| DIODE expr		{ $$=$2; code(xdiode); }
	| CAP  expr		{ $$=$2; code(rcap); }
	| GNDCAP  expr		{ $$=$2; code(xgcap); }
	| BATT expr		{ $$=$2; code(rbatt); }
	| GNDBATT expr		{ $$=$2; code(xgbatt); }
	| buf			{ } 
	| nbuf			{ } 
	;
buf:	  BUF			{ $$=(Inst *)code(xvbuf); } 
	| buf DELAY setexpr 	{ $$=$1; code(xvbufd); }
	| buf BGAIN setexpr 	{ $$=$1; code(xvbufg); }
	| buf BOFFSET setexpr 	{ $$=$1; code(xvbufo); }
 	;
nbuf:	  NBUF			{ $$=(Inst *)code(xnbuf); } 
	| buf BGAIN setexpr 	{ $$=$1; code(xnbufg); }
	| buf BOFFSET setexpr 	{ $$=$1; code(xnbufo); }
	| buf OFFSET setexpr 	{ $$=$1; code(xnbuft); }
	| nbuf ntrans 		{ $$=$1; code(xnbufn); }
 	;
cableparam: LENGTH setexpr		{ $$=$2; code2(xcable,(Inst)$1); }
	| DIA setexpr			{ $$=$2; code2(xcable,(Inst)$1); }
	| DIA2 setexpr			{ $$=$2; code2(xcable,(Inst)$1); }
	| CPLAM setexpr			{ $$=$2; code2(xcable,(Inst)$1); }
	| cableparam LENGTH setexpr	{ code2(xcable,(Inst)$2); }
	| cableparam DIA setexpr	{ code2(xcable,(Inst)$2); }
	| cableparam DIA2 setexpr	{ code2(xcable,(Inst)$2); }
	| cableparam CPLAM setexpr	{ code2(xcable,(Inst)$2); }
	;
sphereparam: expr 		{ code2(xsphere,(Inst)0); }
	| DIA setexpr		{ $$=$2; code2(xsphere,(Inst)$1); }
	;
elecparam: /* nothing */	{ $$=(Inst*)code2(xelec,(Inst)0); }
	| DIA setexpr		{ $$=$2;  code2(xelec,(Inst)$1); }
	| LENGTH setexpr	{ $$=$2;  code2(xelec,(Inst)$1); }
	| RS setexpr		{ $$=$2;  code2(xelec,(Inst)$1); }
	| CAP setexpr		{ $$=$2;  code2(xelec,(Inst)$1); }
	| elecparam DIA setexpr		{ code2(xelec,(Inst)$2); }
	| elecparam LENGTH setexpr	{ code2(xelec,(Inst)$2); }
	| elecparam RS setexpr		{ code2(xelec,(Inst)$2); }
	| elecparam VREST setexpr	{ code2(xelec,(Inst)$2); }
	;
gjtype:  expr 			{ code2(xgj,(Inst)0); }
	| gjparam setexpr	{ $$=$2; code2(xgj,(Inst)$1); }
	| gjtype gjparam setexpr { code2(xgj,(Inst)$2); }
	| gjtype MESGIN secondmsg { code3(xgj,(Inst)$3,(Inst)0); }
	| gjtype MESGIN secondmsg nodenum { code3(xgj,(Inst)(long)$3,(Inst)(long)$4); }
	| gjtype OPEN 		{ code2(xgj,(Inst)$2); }
	| gjtype CLOSE 		{ code2(xgj,(Inst)$2); }
	| gjtype REV 		{ code2(xgj,(Inst)$2); }
	;
gjparam:  AREA 
	| DIA 
	| VGAIN 
	| OFFSET 
	| TAUN 
	| GMAX 
	| GNV 
	;
chantype:  NA 			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  K  			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  KCa 			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  ClCa			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  GLU			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  GABA 		{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  GLY			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  CHRC			{ $$=(Inst *)code2(xchan,(Inst)$1);}
	|  chantype chanparm setexpr  { code2(xchan,(Inst)$2); }
	|  chantype noisetype	     { }
	;
catype:    CA			{ }
	|  psca			{ }
	;
psca:	   CGMP			{ }
	|  NMDA			{ }
	|  AMPA 		{ }
	|  SYN2 		{ }
	;
pscachan:  psca 			{ $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  pscachan chanparm setexpr	{ code2(xcachan,(Inst)$2);}
	|  pscachan cacomp      { }
	;
cachan:    catype 			{ $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  cachan chanparm setexpr	{ code2(xcachan,(Inst)$2);}
	|  cachan cacomparm	     { }
	|  cachan noisetype	     { }
	;
cacomp:    CACOMP	  	  { $$=(Inst*)code2(xcachan,(Inst)$1); }
	|  cacomp cacomparm	  { }
	;
cacomparm: caparm setexpr  	  { $$=(Inst*)code2(xcachan,(Inst)$1); }
	|  caexch
	|  capump
	|  cabuf
	|  cicr
	|  ip3
	|  cacomparm caparm setexpr  { code2(xcachan,(Inst)$2); }
	|  cacomparm caexch
	|  cacomparm capump
	|  cacomparm cabuf
	|  cacomparm cicr
	|  cacomparm ip3
	;
caparm:    CAO 
	|  CAI 
	|  TCAI 
	|  CBOUND 
	|  CASHELL 
	|  CAOSHELL 
	|  MG 
	;
caexch:	   CAEXCH		     { $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  caexch KEX setexpr 	     { code2(xcachan, (Inst)$2); }
	;
capump:	   CAPUMP		     { $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  capump VMAX setexpr 	     { code2(xcachan, (Inst)$2); }
	|  capump KM setexpr   	     { code2(xcachan, (Inst)$2); }
	;
cabuf:	   CABUF		     { $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  cabuf VMAX setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cabuf KD setexpr   	     { code2(xcachan, (Inst)$2); }
	|  cabuf BTOT setexpr        { code2(xcachan, (Inst)$2); }
	|  cabuf BTOTI setexpr       { code2(xcachan, (Inst)$2); }
	;
ip3:	   IP3			     { $$=(Inst *)code2(xcachan,(Inst)$1);}
 	|  ip3 CAS2 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  ip3 IP3I setexpr 	     { code2(xcachan, (Inst)$2); }
	|  ip3 BIP3 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  ip3 VIP3 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  ip3 OI setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 A3 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 D1 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 D2 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 D3 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 D4 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 A2 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 B2 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 K3 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 V3 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 V2 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 V4 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 MTYPE setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 MIP3 setexpr	     { code2(xcachan, (Inst)$2); }
	|  ip3 HIP3 setexpr	     { code2(xcachan, (Inst)$2); }
	;	
cicr:	   CICR			     { $$=(Inst *)code2(xcachan,(Inst)$1);}
	|  cicr CAS setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr VM2 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr VM3 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr KA setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr KF setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr KR setexpr 	     { code2(xcachan, (Inst)$2); }
        |  cicr K1 setexpr           { code2(xcachan, (Inst)$2); }
	|  cicr K2 setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr N setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr M setexpr 	     { code2(xcachan, (Inst)$2); }
	|  cicr P setexpr 	     { code2(xcachan, (Inst)$2); }
        |  cicr C1 setexpr           { code2(xcachan, (Inst)$2); }
	;
chanparm:  TYPE 
	|  VVREV 
	|  OFFSET 
	|  OFFSETM 
	|  OFFSETH 
	|  TAU 
	|  TAUM 
	|  TAUH 
	|  TAUN 
	|  TAUA 
	|  TAUB 
	|  TAUC 
	|  TAUD 
	|  TAUE 
	|  TAUF 
	|  TRCONC 
	|  MAXCOND 
	|  CAPERM 
	|  CAKD 
	|  CAHC 
	|  DENSITY 
	|  NDENSITY 
	|  K1 
	|  K2 
	|  D1 
	|  D2 
	|  N 
	|  UNIT 
	;
membrtype: chantype		{ }
	|  cachan		{ }
	|  memparm setexpr		{ $$=$2; code2(membtyp,(Inst)$1); }
	|  jnoisetype			{ }
	;
jnoisetype:JNOISE setexpr		{ $$=code2(membtyp,(Inst)$1); } 
	|  jnoisetype RSD setexpr	{ code2(membtyp,(Inst)$2); }
	;
memparm:   RM 
	|  CM 
	|  RI 
	|  RG 
	|  VVREV 
	|  VREST 
	;
connect:   AT begin nodenum  	{ $$=$2; code2(conn1,(Inst)(long)$3); } 
	|  AT begin nodenum LOC parglist     /* args backwards here: */ 
				{ $$=$2; code3(conn1l,(Inst)(long)$5,(Inst)(long)$3); } 

	|  AT begin nodenum ':' elemnum OFFSET setexpr PUT nodenum
				{ $$=$2; code3(conn1m,(Inst)(long)$9, (Inst)(long)$3); }

	|  CONNECT begin nodenum TO nodenum   /* args backwards here */
				{ $$=$2; code2(conn2d,(Inst)(long)$5);
					 code2(conn2s,(Inst)(long)$3); } 

	|  CONNECT begin nodenum LOC parglist TO nodenum
				{ $$=$2; code2(conn2d,(Inst)(long)$7);
					 code3(conn2sl,(Inst)(long)$5,(Inst)(long)$3); }

	|  CONNECT begin nodenum TO nodenum LOC parglist
				{ $$=$2; code3(conn2dl,(Inst)(long)$7,(Inst)(long)$5);
					 code2(conn2s,(Inst)(long)$3); }

					/* args backwards here: */	
	|  CONNECT begin nodenum LOC parglist TO nodenum LOC parglist
				{ $$=$2; code3(conn2dl,(Inst)(long)$9,(Inst)(long)$7);
					 code3(conn2sl,(Inst)(long)$5,(Inst)(long)$3); } 
	|  MODIFY expr		{ $$=$2; code2(xmod,(Inst)$1); }
	;
synaptype: SYNAPSE begin		{ $$=$2; }
	|  synaptype OPEN		{ code2(xsynapse,(Inst)$2); }
	|  synaptype CLOSE	 	{ code2(xsynapse,(Inst)$2); }
	|  synaptype SENS synapsens 	{ code2(xsynapse,(Inst)$3); }
	|  synaptype synparm setexpr 	{ code2(xsynapse,(Inst)$2); }
	|  synaptype timectype '[' timecspec ']'
					{ code3(xsynapse,(Inst)$2,(Inst)(long)$4); } 
	|  synaptype timectype setexpr 	{ code3(xsynapse,(Inst)$2,(Inst)1); }
	|  synaptype nfilt	 	{ }
	|  synaptype noisetype	 	{ }
	|  synaptype synapmesng		{ }
	|  synaptype synapresp		{ }
	;
synapsens: V 				{ }
	|  CA 				{ }
	;
synapmesng: MESGOUT ntrans 		{$$=code2(xsynapse,(Inst)$2); }
	;
secondmsg: CAMP				{ }
	|  CGMP				{ }
	;
ntrans:   GLU				{ }
	| AMPA				{ }
	| KAINATE			{ }
	| NMDA				{ }
	| CNQX				{ }
	| GABA				{ }
	| BIC				{ }
	| PTX				{ }
	| GLY				{ }
	| STRY				{ }
	| PH				{ }
	| ATP				{ }
	| secondmsg			{ }
	;
synapresp: RESP begin			{ $$=$2; code2(xsynapse,(Inst)$1); }
	| synapresp pscachan		{ }	
	| synapresp GLU			{ code2(xchan,(Inst)$2); }
	| synapresp GABA		{ code2(xchan,(Inst)$2); }
	| synapresp GLY			{ code2(xchan,(Inst)$2); }
	| synapresp chanparm setexpr	{ code2(xchan,(Inst)$2); }
	;	
synparm:   VVREV			/* list of params that take on values */
	|  THRESH 
	|  TFALL2
	|  TFALL3
	|  HGAIN 
	|  OFFSET 
	|  CGAIN 
	|  COFF 
	|  EGAIN 
	|  VGAIN 
	|  CKD 
	|  CHC 
	|  MAXSRATE 
	|  RRPOOL 
	|  RRPOOLG 
	|  MRRPOOL 
	|  MAXCOND 
	|  TRCONC 
	|  MESGCONC 
	|  KD 
	|  HCOF
	|  EXPON 
	|  LINEAR 
	|  DYAD 
	|  SPOST 
	|  N 
	|  UNIT 
	;
nfilt:     NFILT1 setexpr		{ $$=$2; code2(xsynapse,(Inst)$1); }
	|  NFILT1H setexpr		{ $$=$2; code2(xsynapse,(Inst)$1); }
	|  NFILT2 setexpr		{ $$=$2; code2(xsynapse,(Inst)$1); }
	|  NFILT3 setexpr		{ $$=$2; code2(xsynapse,(Inst)$1); }
	;
timectype: TIMEC1
	|  TIMEC1H 
	|  TIMEC2
	|  TIMEC3 
	;
timecspec: expr 			{ $$ = 1; }
	|  timecspec ',' expr		{ $$ = $1 + 1; } /* how many tau's? */
	;
noisetype: CCHNOISE setexpr 		{ $$=(Inst *)code2(noise,(Inst)$1); }
	|  VESNOISE setexpr 		{ $$=(Inst *)code2(noise,(Inst)$1); }
	|  noisetype N setexpr 		{ code2(noise,(Inst)$2); }
	|  noisetype VSIZE setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype VCOV setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype RSD setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype REFR setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype COV setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype UNIT setexpr 	{ code2(noise,(Inst)$2); }
	|  noisetype TAUF setexpr 	{ code2(noise,(Inst)$2); }
	;
transducer: ROD				{ }
	|   CONE			{ }
	|   CHR				{ }
	|   VTRANSDUCER			{ }
	|   ITRANSDUCER			{ }
	;
receptype: transducer begin parglist  { $$=$2;code3(xphotrec,(Inst)$1,(Inst)(long)$3);}
	|  receptype receparm setexpr	{ code2(xrecparm, (Inst)$2); }
	|  receptype SAVE 		{ code2(xrecparm, (Inst)$2); }
	|  receptype RESTORE 		{ code2(xrecparm, (Inst)$2); }
	;
receparm:  MAXCOND 	{}
	|  DIA 		{}
	|  PIGM 	{}
	|  PATHL 	{}
	|  ATTF 	{}
	|  FILT		{}
	|  TIMEC1	{}
	|  LOOPG	{}
	|  RSD		{}
	|  DARKNOISE	{}
	|  CHANNOISE	{}
	|  PHOTNOISE	{}
	|  UNIT 	{}
	|  LINIT	{}
	|  STIMCHAN	{}
	;
stimulus:  stimparm	{}
	|  stimulus stimparm {}
	;
stimparm:  stimtype     {}
	|  LOC    begin parglist { $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3); }
	|  CENTER begin parglist { $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3); }
	|  START  setexpr        { $$=$2; code2(xstim,(Inst)$1); }
	|  STOP  setexpr         { $$=$2; code2(xstim,(Inst)$1); }
	|  DUR setexpr 	         { $$=$2; code2(xstim,(Inst)$1); }
	|  TFREQ setexpr         { $$=$2; code2(xstim,(Inst)$1); }
	|  DRIFT setexpr         { $$=$2; code2(xstim,(Inst)$1); }
	|  ORIENT setexpr        { $$=$2; code2(xstim,(Inst)$1); }
	|  SPHASE setexpr        { $$=$2; code2(xstim,(Inst)$1); }
	|  XENV setexpr          { $$=$2; code2(xstim,(Inst)$1); }
	|  YENV setexpr          { $$=$2; code2(xstim,(Inst)$1); }
	|  SQ setexpr            { $$=$2; code2(xstim,(Inst)$1); }
	|  BLUR setexpr          { $$=$2; code2(xstim,(Inst)$1); }
	|  SSCALE setexpr        { $$=$2; code2(xstim,(Inst)$1); }
	|  SCATTER begin parglist { $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3); }
	|  stimwave	{}
	|  stiminten	{}
	|  stimfile	{}
	;
stimtype:  BAR expr 	{ $$=$2; code2(xstim,(Inst)$1); }
	|  SPOT expr	{ $$=$2; code2(xstim,(Inst)$1); }
	|  SINE expr	{ $$=$2; code2(xstim,(Inst)$1); }
	|  SINEANN expr	{ $$=$2; code2(xstim,(Inst)$1); }
	|  WINDMILL expr{ $$=$2; code2(xstim,(Inst)$1); }
	|  GABOR expr	{ $$=$2; code2(xstim,(Inst)$1); }
	|  RECT begin arglist 	{ $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3); }
	|  CHECKERBOARD begin arglist { $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3); }
	|  transducer expr 	{ $$=$2; code3(xstim,(Inst)$1,(Inst)1);}
	|  transducer begin dimlist { $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3);}
	|  NODE begin nodenum 	{ $$=$2; code3(xstim,(Inst)$1,(Inst)(long)$3);}
	|  SIMAGE filename 	{ $$=$2; code2(xstim,(Inst)$1);}
	|  SECTOR expr  { $$=$2; code2(xstim,(Inst)$1);}
	;
stimwave:  WAVEL setexpr    { $$=$2; code3(xstim,(Inst)$1,(Inst)0); }
	|  WAVEL SUN 	   { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$2); }
	|  WAVEL '=' SUN   { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$3); }
	|  WAVEL TUNGSTEN  { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$2); }
	|  WAVEL '=' TUNGSTEN  { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$3); }
	|  WAVEL XENON 	   { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$2); }
	|  WAVEL '=' XENON { $$=(Inst *)code3(xstim,(Inst)$1,(Inst)$3); }
	;
stiminten: INTEN setexpr   { $$=$2; code2(xstim,(Inst)$1); }
	|  BACKGR setexpr  { $$=$2; code2(xstim,(Inst)$1); }
	|  VCLAMP setexpr  { $$=$2; code2(xstim,(Inst)$1); }
	|  CCLAMP setexpr  { $$=$2; code2(xstim,(Inst)$1); }
	|  CONTRAST setexpr{ $$=$2; code2(xstim,(Inst)$1); }
	|  PUFF ntrans setexpr  { $$=$3; code3(xstim,(Inst)$1,(Inst)$2); }
	|  PUFF CA setexpr  { $$=$3; code3(xstim,(Inst)$1,(Inst)$2); }
	|  MASK setexpr	    { $$=$2; code2(xstim,(Inst)$1); }
	|  STIMCHAN setexpr { $$=$2; code2(xstim,(Inst)$1); }
	;
stimfile: SFILE filename {$$=$2; code2(xstim,(Inst)$1);}
	;
asgn:	  var '=' expr 	 { code(assign); }
	| var ADDEQ expr { code(addeq); }
	| var SUBEQ expr { code(subeq); }
	| var MULEQ expr { code(muleq); }
	| var DIVEQ expr { code(diveq); } 
	| var ANDEQ expr { code(andeq); } 
	| var OREQ  expr { code(oreq); } 
	;
within:				   { $$ = 0; }	/* nothing */
	| WITHIN2D expr NODE nodenum { $$ = (Inst*)(long) $4; }
	| WITHIN3D expr NODE nodenum { $$ = (Inst*)(long) -$4; }

stmt:     expr { code(expop); } 
	| RETURN { defnonly("return"); $$=code(procret); }
	| RETURN expr 
		{defnonly("return"); $$=$2; code(funcret); }
	| PROCEDURE begin parglist 
		{ $$=$2; code3(call, (Inst)$1, (Inst)(long)$3); }
	| VAR begin parglist 
		{ $$=$2; code3(call, (Inst)$1, (Inst)(long)$3); }
	| defn 	{ } 
	| EXIT { $$=code(xexit); }
	| DEBUGF { $$=code(xdebugf); }
	| FCLOSE '(' expr ')' { $$=$3; code(pfclose); }
	| PRINT prlist { $$ = $2; code(crlf); }
	| PRINTF  begin '(' prflist ')'{ $$=$2; code2(pprintf,(Inst)(long)$4);}
	| SPRINTF begin '(' sprflist ')'{ $$=$2; code2(psprintf,(Inst)(long)$4);}
	| FPRINTF begin '('expr','prflist')'{$$=$2; code2(pfprintf,(Inst)(long)$6);}
	| for '(' forexpr ';' forexpr ';' forexpr ')' stmt end {
		($1)[1] = (Inst)($5-$1-1);	/* cond of loop */
		($1)[2] = (Inst)($7-$1-1);	/* incr of loop */
		($1)[3] = (Inst)($9-$1-1);	/* body of loop */
		($1)[4] = (Inst)($10-$1-1); }	/* end expr */
	| BREAK { $$ = (Inst *)code(breakcode); }
	| CONTINUE { $$ = (Inst *)code(contcode); }
	| while cond stmt end {
		($1)[1] = (Inst)($3-$1-1);	/* body of loop */
		($1)[2] = (Inst)($4-$1-1); }	/* end, if cond fails */
	| if cond stmt end {	 	/* else-less if */
		($1)[1] = (Inst)($3-$1-1);	/* then part */
		($1)[3] = (Inst)($4-$1-1); }	/* end, if cond fails */
	| if cond stmt end ELSE stmt end { /* if with else */
		($1)[1] = (Inst)($3-$1-1);	/* thenpart */
		($1)[2] = (Inst)($6-$1-1);	/* elsepart */
		($1)[3] = (Inst)($7-$1-1); }	/* end, if cond fails */
	| foreach within stmt end { /* look through all nodes */
		($1)[2] = (Inst)0;
		($1)[3] = (Inst)0;
		($1)[4] = (Inst)0;
		($1)[5] = (Inst)0;
		($1)[6] = (Inst)0;
		($1)[7] = (Inst)$2;		/* within node */
		($1)[8] = (Inst)($3-$1-1);	/* body of loop */
		($1)[9] = (Inst)($4-$1-1); }	/* end, when done */
	| foreach NODE snode within end stmt end { /* look through all nodes */
		($1)[2] = (Inst)1;
		($1)[3] = (Inst)(long)$3;
		($1)[4] = (Inst)0;
		($1)[5] = (Inst)0;
		($1)[6] = (Inst)0;
		($1)[7] = (Inst)$4;		/* within node */
		($1)[8] = (Inst)($6-$1-1);	/* body of loop */
		($1)[9] = (Inst)($7-$1-1); }	/* end, when done */
	| foreach NODE snode snode within end stmt end {
		($1)[2] = (Inst)2;
		($1)[3] = (Inst)(long)$3;
		($1)[4] = (Inst)(long)$4;
		($1)[5] = (Inst)0;
		($1)[6] = (Inst)0;
		($1)[7] = (Inst)$5;		/* within node */
		($1)[8] = (Inst)($7-$1-1);	/* body of loop */
		($1)[9] = (Inst)($8-$1-1); }	/* end, when done */
	| foreach NODE snode snode snode within end stmt end {
		($1)[2] = (Inst)3;
		($1)[3] = (Inst)(long)$3;
		($1)[4] = (Inst)(long)$4;
		($1)[5] = (Inst)(long)$5;
		($1)[6] = (Inst)0;
		($1)[7] = (Inst)$6;		/* within node */
		($1)[8] = (Inst)($8-$1-1);	/* body of loop */
		($1)[9] = (Inst)($9-$1-1); }	/* end, when done */
	| foreach NODE snode snode snode snode within end stmt end {
		($1)[2] = (Inst)4;
		($1)[3] = (Inst)(long)$3;
		($1)[4] = (Inst)(long)$4;
		($1)[5] = (Inst)(long)$5;
		($1)[6] = (Inst)(long)$6;
		($1)[7] = (Inst)$7;		/* within node */
		($1)[8] = (Inst)($9-$1-1);	/* body of loop */
		($1)[9] = (Inst)($10-$1-1); }	/* end, when done */
	| stmtlist 		{ }
	| EDIT filename		{ $$=$2; code(edit); }
	| INCLUDE filename	{ $$=$2; code(pushfil); }
	| elemlist		{ } 
	| STIM stimulus		{ $$=$2;    code2(xstim,(Inst)$1); }
	| PLOT plotlist		{ $$=$2; }
	| DISPLAY displist	{ $$=$2;    code2(dispnod,(Inst)$1); }
	| gprim			{ }
	| GRAPH graphlist	{ $$=$2; }
	| FREAD begin '(' filename ',' arrayname ',' vararg ')'
		{ $$=$2; code3(xfread,(Inst)$6,(Inst)(long)$8); $6->type=ARRAY; }
	| FREADS begin '(' filename ',' arrayname ',' vararg ')'
		{ $$=$2; code3(xfreads,(Inst)$6,(Inst)(long)$8); $6->type=ARRAY; }
	| FWRITE begin '(' filename ',' ARRAY ')'
		{ $$=$2; code2(xfwrite,(Inst)$6); }
	| UNLINK begin '(' filename ')'
		{ $$=$2; code(xunlink); }
	| DIM begin arrayname dimlist  { $$=$2; 
	  code3(darray,(Inst)$3,(Inst)(long)$4); $3->type=ARRAY; 
			$3->arrp= (ncarray*)NULL; }
	| DIM begin arrayname dimlist { $$=$2; 
	  code3(darray,(Inst)$3,(Inst)(long)$4); $3->type=ARRAY; 
 		$3->arrp= (ncarray*)NULL; }
		'=' initlist 
		{$$=$2; code3(initarr,(Inst)$3,(Inst)(long)$7);}
	| DIM begin arrayname '[' ']' '=' initlist { $$=$2; 
	  code4(darray,(Inst)$3,(Inst)0,(Inst)(long)$7); $3->type=ARRAY; 
 		$3->arrp= (ncarray*)NULL; code3(initarr,(Inst)$3,(Inst)(long)$7);}
	| fftparm begin '(' ARRAY ')' { $$=$2; code3(dofft,(Inst)$1,(Inst)$4); }
	| LMFIT begin '(' FUNCTION ',' ARRAY ',' ARRAY ')' { 
			$$=$2; code4(do_lmfit,(Inst)$4,(Inst)$6,(Inst)$8); }
	| LMFIT2D begin '(' FUNCTION ',' ARRAY ',' ARRAY ')' { 
			$$=$2; code4(do_lmfit2d,(Inst)$4,(Inst)$6,(Inst)$8); }
	| ELIMIT limitlist 	{ $$=$2; }  
	| INITRAND rinitlist	{ $$=$2; code2(initrgen,(Inst)$1); }
	| ERASE ELEMENT elemnum	{ $$=$3; code(erelem); }
	| ERASE NODE begin nodenum { $$=$3; code2(eranode, (Inst)(long)$4); }
	| ERASE ARRAY		{ $$=(Inst *)code2(erasarr,(Inst)$2); }
	| ERASE VAR		{ $$=(Inst *)code2(erasarr,(Inst)$2); }
	| ERASE MODEL		{ $$=(Inst *)code(eramod); }
	| SAVE MODEL '(' filename ')'     { $$=$4; code(savemod); }
	| RESTORE MODEL '(' filename ')'  { $$=$4; code(restoremod); }
	| RUN 	 		{ $$ = (Inst *)code2(modrun,(Inst)$1); }
	| STEP expr  		{ $$ = $2;     code2(modrun,(Inst)$1); }
	;
gnnlist:  gnnterm		{  }
	| gnnlist ',' gnnterm	{  }
	;
gnnterm:  N setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| NND setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| NNSTD setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| DENSITY setexpr	{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| REG setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| RSD setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| GINFO setexpr		{ $$=$2; code2(xgausnn,(Inst)$1); } 
	| SIZE begin parglist	{ $$=$2; code3(xgausnn,(Inst)$1,(Inst)(long)$3); } 
	| CENTER begin parglist	{ $$=$2; code3(xgausnn,(Inst)$1,(Inst)(long)$3); } 
	;
fftparm:   FFT			{ }
	| IFFT			{ }
	| PFFT			{ }
	| ACOV			{ }
	;
plotype:  V 			{ }
	| VM			{ }
	| I			{ }
	| IM			{ }
	| L			{ }
	| ntrans		{ }
	;
plotypee: FA0			{ }
	| FA1			{ }
	| FA2			{ }
	| FA3			{ }
	| FA4			{ }
	| FA8			{ }
	| FA9			{ }
	| FH0			{ }
	| FH1			{ }
	| FH2			{ }
	| FH3			{ }
	| FH4			{ }
	| FB0			{ }
	| FB1			{ }
	| FB2			{ }
	| FB3			{ }
	| FB4			{ }
	| FC0			{ }
	| FC1			{ }
	| FC2			{ }
	| FC3			{ }
	| FC4			{ }
	| FC9			{ }
	| G0			{ }
	| G1			{ }
	| G2			{ }
	| G3			{ }
	;
ggplot:   VVREV			{ }
	| I			{ }
	| M			{ }
	| H			{ }
	| CA			{ }
	;
plotypeg: G 			{ $$=0; }
	| G '(' ggplot ')'	{ $$=$3->type; }
	;
plotypeh: G '(' expr ')'	{ $$=$3; }
	;
caplot:   VVREV			{ }
	| I			{ }
	| IP			{ }
	| IE			{ }
	| IPE			{ }
	| CAS			{ }
	| CICR			{ }
	;
plotypec: CA 			{ $$=0; }
	| CA '(' caplot	')'	{ $$=$3->type; }
	;
plotyped: CA '(' expr ')'	{ $$=$3; }
	;
plotypeb: CABUF  '(' expr ')'	{ $$=$3; }
	;
plotypef: CABUFB  '(' expr ')'	{ $$=$3; }
	;
plotlist:
	  plotype  begin nodenum { $$=$2; code3(vplot,(Inst)$1,(Inst)(long)$3); }
	| plotypec begin nodenum { $$=$2; 
				code4(vplot,(Inst)CA,(Inst)(long)$3,(Inst)(long)$1); }
	| plotyped nodenum { $$=$1; code4(vplot,(Inst)CA,(Inst)(long)$2,(Inst)1); }
	| plotypeb nodenum { $$=$1; code4(vplot,(Inst)CABUF,(Inst)(long)$2,(Inst)1); }
	| plotypef nodenum { $$=$1; code4(vplot,(Inst)CABUFB,(Inst)(long)$2,(Inst)1); }
	| plotypee begin elemnum { $$=$2; code3(vplot,(Inst)$1,(Inst)1); }
	| plotypeg begin elemnum { $$=$2; 
				code4(vplot,(Inst)G,(Inst)1,(Inst)(long)$1); }
	| plotypeh elemnum { $$=$1; code4(vplot,(Inst)G,(Inst)1,(Inst)1); }
	| V begin '@' CABLE elemnum ':' expr 
					 { $$=$2;code3(vplot,(Inst)$4,(Inst)2);}
	| S var 		{ $$=$2; code3(vplot,(Inst)$1,(Inst)1); }
	| procname 		{ $$=code3(vplot,(Inst)$1,(Inst)0);}
	| plotlist ',' plotype nodenum
				{ $$=$1; code3(vplot,(Inst)$3,(Inst)(long)$4);}
	| plotlist ',' plotypec nodenum
		{ $$=$1; code4(vplot,(Inst)CA,(Inst)(long)$4,(Inst)(long)$3);}
	| plotlist ',' plotyped nodenum
		{ $$=$1; code4(vplot,(Inst)CA,(Inst)(long)$4,(Inst)1);}
	| plotlist ',' plotypeb nodenum
		{$$=$1; code4(vplot,(Inst)CABUF,(Inst)(long)$4,(Inst)1);}
	| plotlist ',' plotypef nodenum
		{$$=$1; code4(vplot,(Inst)CABUFB,(Inst)(long)$4,(Inst)1);}
	| plotlist ',' plotypee elemnum
				{ $$=$1; code3(vplot,(Inst)$3,(Inst)1);}
	| plotlist ',' plotypeg elemnum
		{ $$=$1; code4(vplot,(Inst)G,(Inst)1,(Inst)(long)$3);}
	| plotlist ',' plotypeh elemnum
		{ $$=$1; code4(vplot,(Inst)G,(Inst)1,(Inst)1);}
	| plotlist ',' V begin '@' CABLE elemnum ':' expr 
				{ $$=$1; code3(vplot,(Inst)$6,(Inst)2);}
	| plotlist ',' S var 	{ $$=$1; code3(vplot,(Inst)$3,(Inst)1); }
	| plotlist ',' procname { $$=$1; code3(vplot,(Inst)$3,(Inst)0); }
	| plotlist maxmin	{ $$=$1; code(vplotm); } 
	| plotlist plltyp 	{ $$=$1; }
	;
gprim:    GMOVE begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GDRAW begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GRMOVE begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GRDRAW begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GCIRC begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GRECT begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GPEN  begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GVPEN begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GROT  begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GORIG begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GCWID begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GCROT begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GDASH begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GSIZ begin parglist	{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GWINDOW begin parglist{ $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GFRAME begin parglist { $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GPURGE begin parglist { $$=$2;    code3(gplot,(Inst)$1,(Inst)(long)$3); }
	| GTEXT begin '(' prflist ')'{ $$=$2; code2(txtf,(Inst)(long)$4);}
	| GLABEL begin '(' prflist ')'{ $$=$2; code2(gglabel,(Inst)(long)$4);}
	;
graphlist: graphtype		{ $$=$1; }
	| graphlist graphtype 	{ $$=$1; }
	;
graphtype: begin parglist 	{ $$=$1;      code3 (grph, (Inst)0L,(Inst)(long)$2);}
	| INIT			{ $$=(Inst *) code3 (grph,(Inst)$1, (Inst)0L); }
	| RESTART		{ $$=(Inst *) code3 (grph,(Inst)$1, 0L); }
	| PEN begin parglist	{ $$=$2;      code3 (grph,(Inst)$1,(Inst)(long)$3);}
	| X maxmin 		{ $$=$2;      code3 (grph,(Inst)$1, 0L); }
	| X maxmin grphlst	{ $$=$2;      code3 (grph,(Inst)$1, 0L); }
	| Y maxmin		{ $$=$2;      code3 (grph,(Inst)$1, 0L); }
	| Y maxmin grphlst	{ $$=$2;      code3 (grph,(Inst)$1, 0L); }
	;
maxmin:   MAX setexpr MIN setexpr { $$=$2; }
	| MIN setexpr MAX setexpr { $$=$2; }
	;
grphlst:  grltyp		{  }
	| grphlst grltyp	{  }
	;
grltyp:   CHAR setexpr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)-1); }
	| CCHAR setexpr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)-1); }
	| PEN  setexpr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)-1); }
	| VPEN procname		{ $$=    code3 (grph,(Inst)$1, (Inst)$2); }
	| SIZE setexpr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)-1); } 
	| FILT begin '[' timecspec ']'{$$=$2;code3 (grph,(Inst)$1,(Inst)(long)$4);} 
	| FILT expr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)1); } 
	| PLNAME setexpr	{ $$=$2; code3 (grph,(Inst)$1, (Inst)1); } 
	| PLNUM  setexpr	{ $$=$2; code3 (grph,(Inst)$1, (Inst)1); } 
	| PLSIZE setexpr	{ $$=$2; code3 (grph,(Inst)$1, (Inst)1); } 
	| PLVAL setexpr		{ $$=$2; code3 (grph,(Inst)$1, (Inst)1); } 
	;
plltyp:   CHAR setexpr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)-1); }
	| CCHAR setexpr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)-1); }
	| PEN  setexpr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)-1); }
	| VPEN procname		{ $$=    code3 (plparams,(Inst)$1, (Inst)$2); }
	| SIZE setexpr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)-1); } 
	| FILT begin '[' timecspec ']'{$$=$2;code3 (plparams,(Inst)$1,(Inst)(long)$4);} 
	| FILT expr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	| PLNAME setexpr	{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	| PLNUM  setexpr	{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	| PLSIZE setexpr	{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	| PLVAL setexpr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	| PLARR expr		{ $$=$2; code3 (plparams,(Inst)$1, (Inst)1); } 
	;
displist: disptype		{ $$=$1; }
	| displist disptype 	{ $$=$1; }
disptype: MATCHING begin nodenum
			  { $$=$2; code4(dispnod,(Inst)$1,(Inst)0,(Inst)(long)$3); }
	| CONNECT begin nodenum TO nodenum
			  { $$=$2; code4(dispnod,(Inst)$1,(Inst)(long)$5,(Inst)(long)$3); }
	| RANGE begin nodenum TO nodenum
			  { $$=$2; code4(dispnod,(Inst)$1,(Inst)(long)$5,(Inst)(long)$3); }
	| ELEMENT elemnum { $$=$2; code4(dispnod,(Inst)$1,(Inst)0,(Inst)1); }
	| elemtype        { $$=(Inst *)code2(dispnod,(Inst)$1); }
	| COMPS           { $$=(Inst *)code2(dispnod,(Inst)$1); }
	| EXCEPT 	  { $$=code2(dispnod,(Inst)$1); } 
	| rotatype setexpr{ $$=$2; code2(dispnod,(Inst)$1); }
	| SIZE setexpr	  { $$=$2; code2(dispnod,(Inst)$1); } 
	| WINDOW begin parglist  { $$=$2; code3(dispnod,(Inst)$1,(Inst)(long)$3); } 
	| DSCALE setexpr  { $$=$2; code2(dispnod,(Inst)$1); } 
	| COLOR setexpr	       { $$=$2; code3(dispnod,(Inst)$1,(Inst)0); } 
	| COLOR setexpr	maxmin { $$=$2; code3(dispnod,(Inst)$1,(Inst)1); } 
	| VPEN procname	   { $$= code3(dispnod,(Inst)$1,(Inst)$2); } 
	| CALIBLIN setexpr { $$=$2; code3(dispnod,(Inst)$1,(Inst)0); } 
	| CALIBLIN setexpr LOC parglist
			   { $$=$2; code3(dispnod,(Inst)$1,(Inst)(long)$4);} 
	| HIDE 	 	  { $$=(Inst *)code2(dispnod,(Inst)$1); } 
	| CENTER begin parglist { $$=$2; code3(dispnod,(Inst)$1,(Inst)(long)$3); } 
	| RMOVE  begin parglist { $$=$2; code3(dispnod,(Inst)$1,(Inst)(long)$3); } 
	| ONLY		  { $$=(Inst *)code2(dispnod,(Inst)$1); }
	| NODE	 	 { $$=(Inst *)code2(dispnod,(Inst)$1); }
	| STIM AT expr 	 { $$=$3; code3(dispnod,(Inst)$1,(Inst)1); }
	| STIM  	 { $$=0; code3(dispnod,(Inst)$1,(Inst)0); }
	| CMAP setexpr	 { $$=$2; code2(dispnod,(Inst)$1); }
	| NEWPAGE 	 { $$=(Inst *)code2(dispnod,(Inst)$1); }
	| Z maxmin  	 { $$=(Inst *)code3(dispnod,(Inst)$1,(Inst)1); }
	;
filename: expr 	     	{ $$ = $1; }
	;
cond:	'(' cexpr ')'	{code(STOPC); $$ = $2; }
	;
while:	  WHILE { $$ = (Inst *)code3(whilecode,STOPC,STOPC); }
	;
forexpr:  cexpr 	{code(STOPC); $$ = $1; }
	;
for:	  FOR { $$ = (Inst *)code(forcode); code4(STOPC,STOPC,STOPC,STOPC); }
	;
if:	  IF	{ $$ = (Inst *)code(ifcode); code3(STOPC,STOPC,STOPC); }
	;
foreach:  FOREACH 		  { $$ = (Inst *)code(foreacode);
				    code((Inst)0);
				    code4(STOPC,STOPC,STOPC,STOPC);
				    code4(STOPC,STOPC,STOPC,STOPC); }
	| FOREACH elemtypee '?' var { $$ = (Inst *)code(foreacode);
				    code((Inst)$2);
				    code4(STOPC,STOPC,STOPC,STOPC);
				    code4(STOPC,STOPC,STOPC,STOPC); }
	;
elemtype: CABLE 		{ }
	| SPHERE		{ }
	| SYNAPSE		{ }
	| CHAN 			{ }
	| transducer		{ }
	| GJ			{ }
	| PNX			{ }
	| LOAD 			{ }
	| RESISTOR 		{ }
	| DIODE 		{ }
	| CAP  			{ }
	| GNDCAP  		{ }
	| BATT 			{ }
	| GNDBATT 		{ }
	| ELECTRODE 		{ }
	| BUF			{ }
	| DYAD			{ }	/* not strictly an element type */
	| EXPON			{ }	/* not strictly an element type */
	| LINEAR		{ }	/* not strictly an element type */
	;
elemtypee: elemtype		{ }
	| ELEMENT		{ }
	;	
limitlist: limitype 		{ $$=$1; }
	|  limitlist limitype 	{ $$=$1; }
	;
limitype: X maxmin		{ $$=$2; code2(elimit,(Inst)$1); }
	| Y maxmin		{ $$=$2; code2(elimit,(Inst)$1); }
	| Z maxmin		{ $$=$2; code2(elimit,(Inst)$1); }
	| ELEMENT elemnum	{ $$=$2; code2(elimit,(Inst)$1); }
	;
rinitlist: expr			{ $$=$1; code2(initrgen,(Inst)0); } 
	| rinitlist RSD setexpr	{ $$=$1; code2(initrgen,(Inst)$2); }
	;
begin:	  /* nothing */		{ $$ = progp; }
	;
end:	  /* nothing */		{ code(STOPC); $$ = progp; }
	;
stmtlist: '{' stmtls '}'	{ $$ = $2; }
	| '{' begin localdef listerm stmtls '}'	
		{ $$ = $2; code2(locend,(Inst)(long)(-$3+1+argcount));
		   erasarg($3); argcount = $3 - 1; }
	;
stmtls:   /* nothing */		{ $$ = progp; }
	| stmtls stmt listerm 
	;
expr:	  NUMBER     { $$ = (Inst *)code2(constpush, (Inst)$1); }
	| CONST	     { $$ = 
		(Inst *)code4(varpush,(Inst)(long)$1->type,(Inst)$1,evalvar);}
	| initlist  { $$ = 
		(Inst *)code4(varpush,(Inst)ARRAYCONST,(Inst)(long)$1,evalvar);} 
	| DIMS '(' var ')' {  code(xdims); }
	| LITCHAR    { $$ = (Inst *)code2(constpush,(Inst)$1); }
	| STRING     { $$ = (Inst *)code2(constpush,(Inst)$1); }
	| var  	     { $$ = $1; code(evalvar); }
	| var INCROP { $$ = $1; code(postinc); } 
	| var DECROP { $$ = $1; code(postdec); }
	| INCROP var { $$ = $2; code(preinc); }
	| DECROP var { $$ = $2; code(predec); }
	| asgn
	| FUNCTION begin parglist
		     { $$ = $2; code3(call,(Inst)$1,(Inst)(long)$3); }
	| READ '(' var ')' { $$=$3; code(varread); }
	| NTYPE '(' elemtype ')' { $$=code2(gettype,(Inst)$3); }
	| BLTIN begin parglist { 
			 $$=$2; code3(bltin,(Inst)(long)$3,(Inst)$1->ptr);}
	| GAUSNN begin '(' arrayname ',' gnnlist ')'
			 { $$=$2; code3(xgausnn,(Inst)$1,(Inst)$4);
				$4->type = ARRAY; }
	| NOTINIT '(' var ')'	 { $$=$3; code(notinitx); }
	| VARNUM  '(' var ')' 	{ $$=$3; code(varnum); }
	| VARSTR  '(' var ')' 	{ $$=$3; code(varstr); }
	| VARCHR  '(' var ')' 	{ $$=$3; code(varchr); }
	| FOPEN  '(' filename ',' expr ')'   { $$=$3; code(pfopen); } 
	| FGETS  '(' var ',' expr','expr ')' { $$=$3; code(pfgets); }
	| FGETC  '(' var ',' expr ')'        { $$=$3; code(pfgetc); }
	| FPUTC  '(' var ',' expr ')'        { $$=$3; code(pfputc); }
	| GETFLDS '(' expr ')' { $$=$3; code(getflds); }
	| SCANF  begin '(' scnflist ')'      { $$=$2; code2(pscanf,(Inst)(long)$4);}
	| FSCANF begin '(' expr ','scnflist')'{$$=$2; code2(pfscanf,(Inst)(long)$6);}
	| SSCANF begin'(' expr ',' scnflist ')'{$$=$2; code2(psscanf,(Inst)(long)$6);}
	| comndline		{ } 
	| plotype begin nodenumex
				 { $$=$2; code3(xrecord, (Inst)$1, (Inst)(long)$3); }
	| plotypec begin nodenumex
			{ $$=$2; code4(xrecord, (Inst)CA, (Inst)(long)$3, (Inst)(long)$1); }
	| plotyped nodenumex
			{ $$=$1; code4(xrecord, (Inst)CA, (Inst)(long)$2, (Inst)1); }
	| plotypeb nodenumex
			{ $$=$1; code4(xrecord, (Inst)CABUF,(Inst)(long)$2,(Inst)1);}
	| plotypef nodenumex
			{ $$=$1; code4(xrecord, (Inst)CABUFB,(Inst)(long)$2,(Inst)1);}
	| plotypee begin elemnum
				 { $$=$2; code3(xrecord, (Inst)$1, (Inst)1); }
	| plotypeg begin elemnum
			{ $$=$2; code4(xrecord, (Inst)G, (Inst)1, (Inst)(long)$1); }
	| plotypeh elemnum
			{ $$=$1; code4(xrecord, (Inst)G, (Inst)1, (Inst)1); }
	| V begin '@' CABLE elemnum ':' expr 
				 { $$=$2;code3(xrecord,(Inst)$4,(Inst)2);}
	| EEDIST begin '(' elemnum ',' elemnum ')'
				 { $$=$2; code(eedist); }
	| E3DIST begin '(' nodenum ',' elemnum ')'
				 { $$=$2; code2(e3dist, (Inst)(long)$4); }
	| E2DIST begin '(' nodenum ',' elemnum ')'
				 { $$=$2; code2(e2dist, (Inst)(long)$4); }
	| EZDIST begin '(' nodenum ',' elemnum ')'
				 { $$=$2; code2(ezdist, (Inst)(long)$4); }
	| EFRAC begin '(' nodenum ',' elemnum ')'
				 { $$=$2; code2(efrac, (Inst)(long)$4); }
	| N3DIST begin '(' nodenum ',' nodenum ')'  /* args are backward here: */
				 { $$=$2; code3(n3dist, (Inst)(long)$6, (Inst)(long)$4); }
	| N2DIST begin '(' nodenum ',' nodenum ')'  /* args are backward here: */
				 { $$=$2; code3(n2dist, (Inst)(long)$6, (Inst)(long)$4); }
	| NZDIST begin '(' nodenum ',' nodenum ')'  /* args are backward here: */
				 { $$=$2; code3(nzdist, (Inst)(long)$6, (Inst)(long)$4); }
	| ELEMENT elemnum PFIELD elemfield { $$=$2; code2(efield, (Inst)$4); }
	| CHAN    elemnum PFIELD chanfield { $$=$2; code2(cfield, (Inst)$4); }
	| NODE begin nodenum PFIELD nodefield 
				{ $$=$2; code3(nfield,(Inst)$5, (Inst)(long)$3);}
	| '(' cexpr ')'	{ $$ = $2; }
	| expr '+' expr { code (add); }
	| expr '-' expr { code (sub); }
	| expr '*' expr { code (mul); }
	| expr '/' expr { code (xdiv); }
	| expr '%' expr { code (mod); }
	| expr '^' expr { code (power); }
	| '-' expr %prec UNARYMINUS	{ $$=$2; code(negate); }
	| '+' expr %prec UNARYPLUS	{ $$=$2; }
	| expr GT expr  { code (gt); }
	| expr GE expr  { code (ge); }
	| expr LT expr  { code (lt); }
	| expr LE expr  { code (le); }
	| expr EQ expr  { code (eq); }
	| expr NE expr  { code (ne); }
	| expr AND expr { code (xand); }
	| expr OR expr  { code (orx); }
	| expr BITAND expr { code (bitand_x); }
	| expr BITOR expr { code (bitor_x); }
	| expr XOR expr { code (xxor); }
	| NOT expr 	{ $$ = $2; code (xnot); }
	| SIZEOF '(' var  ')' { $$=$3; code (xsizeof); }
	;
comndline: COMNDLINE		{ $$=code2(xsystem, (Inst)$1); } 
	| SYSTEM '(' expr ')'	{ $$=$3; code2(xsystem, (Inst)0); } 
	;
commaexpr: stmt  ',' expr	{ } 
	| commaexpr ',' expr	{ code(popone); } 
	;
cexpr:	  expr			{ }
	| commaexpr		{ }
	;
setexpr:  expr		{ $$ = $1; }
	| "=" expr	{ $$ = $2; }
	;
var:	  VAR	   { $$=(Inst*)code3(varpush,(Inst)(long)$1->type,(Inst)$1); }
	| ARRAY begin dimlist{ $$ = $2;
		code4(varpush,(Inst)ARRAYVAL,(Inst)$1,(Inst)(long)$3);}

	| ARRAY { $$ = code4(varpush,(Inst)(long)$1->type,(Inst)$1,(Inst)0);}

	| LOCALARR { $$ = code4(varpush,(Inst)(long)$1->type,(Inst)(long)$1->argnum,(Inst)0);}

	| LOCALARR begin dimlist{ $$ = $2; 
		code4(varpush,(Inst)LOCALARRVAL,(Inst)(long)$1->argnum,(Inst)(long)$3);}

	| ARG	   {$$=(Inst *)code3(varpush,(Inst)(long)$1->type,(Inst)(long)$1->argnum);}
	| ARG begin dimlist {$$ = $2;
		code4(varpush,(Inst)ARGARRVAL,(Inst)(long)$1->argnum,(Inst)(long)$3);}
	| LOCALVAR {$$=(Inst *)code3(varpush,(Inst)(long)$1->type,(Inst)(long)$1->argnum);}
	;
prlist:	  expr			{ code(prexpr); }
	| prlist ',' expr	{ code(prexpr); }
	;
prflist:  expr 			{ $$ =  0; }
        | expr ',' prfargs 	{ $$ = $3; }
	;
sprflist: var ',' expr	 	{ $$ =  0; }
        | var ',' expr ',' prfargs { $$ = $5; }
	;
prfargs:  expr			{ $$ = 1; }
	| prfargs ',' expr	{ $$ = $1 + 1; }	/* how many args? */
	;
scnflist:  expr 		{ $$ =  0; }
        | expr ',' scnfargs 	{ $$ = $3; }
	;
scnfargs:  var			{ $$ = 1; }
	| scnfargs ',' var	{ $$ = $1 + 1; }	/* how many args? */
	;
defn:	  FUNC procname 
		{ $2->type=FUNCTION; indef=1; formal=1; starg = argcount; }
	 '(' formarg ')' { formal=0; }
          stmtlist {code(procret); define($2,(int)$5, $8, starg);
		erasarg(starg+1); argcount=starg; indef=0; }
	| PROC procname 
		{ $2->type=PROCEDURE; indef=1;formal=1; starg = argcount; }
	 '(' formarg ')' { formal=0; }
	 stmtlist {code(procret); define($2,(int)$5, $8, starg);
		erasarg(starg+1); argcount=starg; indef=0;}
	| PROC procname '=' PROCEDURE 
			{$2->type=PROCEDURE; define($2,-1,(Inst*)$4,starg);}
	| FUNC procname '=' FUNCTION  
			{$2->type=FUNCTION; define($2,-1,(Inst*)$4,starg);} 
	| FUNC procname '=' BLTIN  
			{$2->type=BLTIN; define($2,-2,(Inst*)$4,starg);} 
	;
/*
dimlistm: '[' expr ']'		{ $$ = 3; }
	| '[' ':' ']'		{ $$ = 2; }
	| dimlist '[' expr ']' 	{ $$ = $1<<1 + 1; }
	| dimlist '[' ':'  ']' 	{ $$ = $1<<1; }
	;
*/
dimlist:  '[' expr ']'		{ $$ = 1; }
	| dimlist '[' expr ']' 	{ $$ = $1 + 1; }
	;
nodenum: expr			{ $$ = 1; }
	| dimlist		{ $$ = $1; }
	;
nodenumex: dimlist		{ $$ = $1; }
	;
elemnum: var			{ $$ = $1; code(evalvar); }
	| NUMBER		{ $$ = (Inst *)code2(constpush, (Inst)$1); }
	| '(' expr ')'		{ $$ = $2; }
	;
snode: 	  '?' var		{ $$ = 1; }
	| '[' expr ']'		{ $$ = 0; }
	;
elemfield: TYPE			{ }
	|  NTYPE		{ }
	|  REGION		{ }
	|  ELABL		{ }
	|  LENGTH		{ }
	|  N3DIST		{ }
	|  N2DIST		{ }
	|  DIA			{ }
	|  DIA2			{ }
	|  RM			{ }
	|  RI			{ }
	|  CM			{ }
	|  CPLAM		{ }
	|  NODE1A		{ }
	|  NODE1B		{ }
	|  NODE1C		{ }
	|  NODE1D		{ }
	|  NODE2A		{ }
	|  NODE2B		{ }
	|  NODE2C		{ }
	|  NODE2D		{ }
	|  DYAD			{ }
	|  SDYAD		{ }
	|  SCURVE		{ }
	|  NUMDYAD		{ }
	|  NUMSPRE		{ }
	|  MODIFY		{ }
	;
chanfield: TYPE			{ }
	|  NTYPE		{ }
	|  STYPE		{ }
	|  NSTATE		{ }
	;
nodefield: NUMCONN		{ }
	| NUMSYN		{ }
	| EXIST			{ }
	| XLOC			{ }
	| YLOC			{ }
	| ZLOC			{ }
	| CACOMP		{ }
	| expr			{ $$ =  (Inst *)0; }	/* relative elem # */
	| SYNAPSE expr		{ }			/* relative syn # */
	;
procname: VAR
	| FUNCTION
	| PROCEDURE
	| BLTIN
	;
rotatype: XROT			{ }
	| YROT			{ }
	| ZROT			{ }
	;
undefarrname:  VAR
	;
arrayname:  undefarrname	{ }
	|   ARRAY		{ }
	;
arglist:   /* nothing */	{ $$ = 0; }
	| expr 			{ $$ = 1; }
	| arglist ',' expr	{ $$ = $1 + 1; }
	;
parglist: '(' arglist ')'	{ $$ = $2; }
	;
initlist: '{''{' arglist '}''}'	{ $$ = $3; }
	| '{' '{' expr ':' expr '}''}' { $$ = -2; }
	| '{' '{' expr ':' expr ':' expr '}' '}'{ $$ = -3; }
	;
formarg:  /* nothing */		{ $$ = 0; }
	| ARG			{ $$ = 1;      $1->argnum = ++argcount; } 
	| formarg ',' ARG	{ $$ = $1 + 1; $3->argnum = ++argcount; }
	;
localarrdef: DIM LOCALVAR {$$=$2; savlocal=inlocal; inlocal=0;}
	;
localvar: LOCALVAR  		{ $$ = 1; $1->argnum = ++argcount; } 
	| localarrdef dimlist 
			{ inlocal=savlocal; $$ = 1; 
			code3(dlocarray,(Inst)(long)++argcount,(Inst)(long)$2); 
			 $1->argnum=argcount; arrcount++;
			 $1->type=LOCALARR; } 

	| localarrdef dimlist 
		{ code3(dlocarray,(Inst)(long)++argcount,(Inst)(long)$2); 
		 $1->argnum=argcount; arrcount++;
		 $1->type=LOCALARR; } 
		'=' initlist 
		{ $$ = 1; 
		inlocal=savlocal; code3(initarr,(Inst)0,(Inst)(long)$5);}

        | localarrdef '[' ']' '=' initlist
		{ $$=1; code4(dlocarray,(Inst)(long)++argcount,(Inst)0,(Inst)(long)$5);
		 code2(stackmove,(Inst)(long)$5); 
		$1->argnum=argcount; arrcount++; $1->type=LOCALARR; 
		inlocal=savlocal; code3(initarr,(Inst)0,(Inst)(long)$5); }

	| localvar ',' LOCALVAR { $$ = $1 + 1; $3->argnum = ++argcount; }
	| localvar ',' localarrdef dimlist 
			{ inlocal=savlocal; $$ = $1 + 1; 
			code3(dlocarray,(Inst)(long)++argcount,(Inst)(long)$4); 
			  $3->argnum=argcount; arrcount++;
			  $3->type=LOCALARR; } 
	;
	| localvar ',' localarrdef dimlist 
		{ code3(dlocarray,(Inst)(long)++argcount,(Inst)(long)$4); 
		  $3->argnum=argcount; arrcount++; $3->type=LOCALARR; } 
			'=' initlist 
		{ $$ = $1 + 1; 
		inlocal=savlocal; code3(initarr,(Inst)0,(Inst)(long)$7);}

        | localvar ',' localarrdef '[' ']' '=' initlist
		{ $$=$1+1; code4(dlocarray,(Inst)(long)++argcount,(Inst)0,(Inst)(long)$7);
		 code2(stackmove,(Inst)(long)$7); 
		$3->argnum=argcount; arrcount++; $3->type=LOCALARR; 
		inlocal=savlocal; code3(initarr,(Inst)0,(Inst)(long)$7); }

	;
localdef: LOCAL { inlocal=argcount+1; arrcount=0; } localvar 
	    { $$=argcount-$3+1; inlocal=0; 
		(Inst *)code3(local,(Inst)(long)$3,(Inst)(long)arrcount); }
        | localdef listerm LOCAL { inlocal=argcount+1; arrcount=0; } localvar 
    	    { $$=$1; inlocal=0; (Inst *)code3(local,(Inst)(long)$5,(Inst)(long)arrcount); }
	;
vararg:  var			{ $$ = 1; }
	| vararg ',' var	{ $$ = $1 + 1; }
	;
%%
	/* end of grammar */
