/* Segment ncslv in Program nc */

/* Runs slave task for nc. */
/* Receives compartments from master nc task, */
/*   and runs compartment computations in parallel. */

#include "ncomp.h"
#include "ncvirt.h"

virtobj *virtcomp_list;

updatobj *recvlist;
updatobj *sendlist;

int cumvirt;
int cumrecv;
int cumsend;

#define NPROC 5

/*------------------------------------------------*/

main(int argc,char **argv)

{
   int mytid, master, task_ids[NPROC], retfl;

  mytid = pvm_mytid();
  master = pvm_parent();

  pvm_recv (master, MINIT);
  pvm_upkdouble (&timinc,1,1);
  pvm_upkbyte(&mcomnd,1,1);
  for (retfl=0; !retfl; ) {
    pvm_upkbyte(&mcomnd,1,1);
    switch {mcomnd) {

     case MCOMP:  recvcomp();
                  break;

     case MRUN:   runiter();
                  break;

     case MEXIT:  retfl = 1;
                  break;
    }
  }
}

/*------------------------------------------------*/
 
