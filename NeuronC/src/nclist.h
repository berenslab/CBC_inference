/* Header file nclist.h for program "nc". */
/* Contains list pointers */


int cumelem=0;
int cumattr=0;
int cumnattr=0;
int cumcattr=0;
int cumnode=0;
int cumcomp=0;
int cumconn=0;
int cumchan=0;
int cumsynap=0;
int cumload=0;
int cumphotrec=0;
int cumrecst=0;
int cumcacomp=0;
int cumcycacomp=0;
int cumcycgcomp=0;
int cumclst=0;
int cumdlst=0;
int cumrand=0;
int cumlpfilt=0;
int cumntcomp=0;
int cumgj=0;
int cumdbuf=0;
int cumnbuf=0;
int reccum = 0;
int delelmt = 0;
int ndomains = 0;

comp *compnt=0;
comp *compend=0;
conn *connpnt=0;
conn *connend=0;
synap *synpnt=0;
synap *synend=0;
chan *chanpnt=0;
chan *chanend=0;
load *loadpnt=0;
load *loadend=0;
photrec *recpnt=0;
photrec *recend=0;
recstim *recspnt=0;
recstim *recsend=0;

#include "stim.h"
extern recnod *reclist;
extern recnod *reclend;

