/* Segment ncvirt in Program nc */

/* Runs communication between nc parallel tasks. */
/* Moves compartments between tasks, and sets up */
/*   send and receive lists. */

#include "ncomp.h"
#include "ncvirt.h"
#include "ncio.h"

virtobj *virtcomp_list;

updatobj *recvlist;
updatobj *sendlist;

int cumvirt;
int cumrecv;
int cumsend;

/*------------------------------------------------*/

void zerocpnt ( xxx connpnt, comp *cpnt)

/* Zero the comp pointer in a connection */
/*  Also set the comp number */

{
 if (!connpnt || !cpnt) return;
 if (connpnt->comp1 == cpnt) {
    connpnt->comp1 = NULL;
    connpnt->num1 = cpnt->num;
 }
 if (connpnt->comp2 == cpnt) {
    connpnt->comp2 = NULL;
    connpnt->num2 = cpnt->num;
 }
}

/*------------------------------------------------*/

void zeroclst (conlst *pnt, comp *cpnt)

/* Zero the pointers that point to a comp in its list of
connections (conlst). Then delete the list.   */

{
   conlst *lpnt,*tlpnt;

   for (lpnt=pnt; lpnt; ) {
        tlpnt = lpnt;
        lpnt = lpnt->next;
        zerocpnt (tlpnt->conpnt,cpnt);
        efree (tlpnt);
    }
}

/*------------------------------------------------*/

void zeronlst (conlst *pnt, comp *cpnt)

/* Zero the pointers that point to a comp in its list of
nodes (nodlst). Then delete the list.   */

{
   conlst *lpnt,*tlpnt;

   for (lpnt=pnt; lpnt; ) {
        tlpnt = lpnt;
        lpnt = lpnt->next;
        if (tlpnt->conpnt) {
          ((node *)tlpnt->conpnt)->comptr = NULL;
          ((node *)tlpnt->conpnt)->compnum = cpnt->num;
	}
        efree (tlpnt);
    }
}

/*------------------------------------------------*/

void virtcomp (comp *cpnt, int tid)

/* Set a compartment to virtual, and save its address in 
   master list.
*/

{
    virtobj *vpnt;

 cpnt->virt = 1;			/* set type to virtual */
 if ((vpnt=(virtobj *)emalloc(sizeof(virtobj))) == NULL) {
    ncfprintf (stderr,"no space left for virtobj %d\n",cumvirt);
    return (NULL);
 }
 vpnt->next = virtcomp_list;
 vpnt->last = NULL;
 if (virtcomp_list) virtcomp_list->last = vpnt;
 virtcomp_list = vpnt;
 vpnt->structype = VCOMP;
 vpnt->num = cpnt->num; 
 vpnt->tid = tid;
 cumvirt++;
}

/*------------------------------------------------*/

sendcomp (comp* cpnt; int tid)

/* Send a compartment to a pvm task */

{
    updatobj *upnt;

 if ((upnt=(updatobj *)emalloc(sizeof(updatobj))) == NULL) {
    ncfprintf (stderr,"no space left for updatobj %d\n",cumsend);
    return (NULL);
 }
 upnt->next = sendlist;
 upnt->last = NULL;
 if (sendlist) sendlist->last = upnt;
 sendlist = upnt;
 upnt->structype = VCOMP;
 upnt->field = VVOLT;
 upnt->compnum = cpnt->num; 
 upnt->tid = tid;
 cumsend++;
 pvm_initsend(PvmDataDefault);
 pk_comp(cpnt);

}

/*------------------------------------------------*/

movecomp (comp* cpnt; int tid)

/* Move a compartment from one pvm task to another */

{

  sendcomp (cpnt,tid);
/*  zeroclst (cpnt->clst);	/* zero pointers in connection list */
/*  zeronlst (cpnt->nodlst);	/* zero pointers in node list */
  virtcomp (cpnt,tid);		/* set the compartment to virtual */
}

/*------------------------------------------------*/
 
