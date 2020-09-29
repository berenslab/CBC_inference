/* cone function */

extern double conerest;
extern double conerm;

int mcone (double xpos, double ypos, int n) 

{ 
   double dia;
   photorec *p;
   cable *c;
   sphere *s;

 p = make_cone (nd(n), dia=1.674); p->xpos=xpos; p->ypos=ypos; p->attf=1.5; 

 /*c = make_cable(nd(n), nd(n,1) c->dia=10;c->length=10;c->Ri=5000;c->Rm=1e4;/* 5 x rod os area*/ 

 s = make_sphere(nd(n), dia=16); s->Rm=500000; s->vrest=conerest; 

 /*c = make_cable(nd(n)), nd(n+1));c->dia=.1;c->length=.2;c->Rm=conerm;/* connecting cilium */
 
 s = make_sphere(nd(n), dia=3); s->Rm=conerm; s->vrest=conerest;

 c = make_cable(nd(n),nd(n,1));c->dia=1;c->length=50;c->Rm=conerm; c->vrest=conerest;
 s = make_sphere(nd(n,1), dia=5); s->Rm=conerm; s->vrest=conerest; 

 return (1);     /* return minor node num for cone pedicle */
};

int mconep (double xpos, double ypos, int n, int pigm, double linit) 

{ 
   double dia;
   photorec *p;
   cable *c;
   sphere *s;

 p = make_cone (nd(n), pigm, dia=1.674); 
 p->xpos=xpos; p->ypos=ypos; p->attf=1.5; p->linit=linit; 
 if (pigm<19 || pigm>21) {
    p->timec1=.5; 		/* speed up cone response 2x */
    p->loopg = 1.0;		/* gain of ca loop, set this from 0.4 to 1.0 */
 }

 /*c = make_cable(nd(n), nd(n,1) c->dia=10;c->length=10;c->Ri=5000;c->Rm=1e4;/* 5 x rod os area*/ 

 s = make_sphere(nd(n), dia=16); s->Rm=500000; s->vrest=conerest; 

 /*c = make_cable(nd(n)), nd(n+1));c->dia=.1;c->length=.2;c->Rm=conerm;/* connecting cilium */
 
 s = make_sphere(nd(n), dia=3); s->Rm=conerm; s->vrest=conerest;

 c = make_cable(nd(n),nd(n,1));c->dia=1;c->length=50;c->Rm=conerm; c->vrest=conerest;
 s = make_sphere(nd(n,1), dia=5); s->Rm=conerm; s->vrest=conerest; 

 return (1);     /* return minor node num for cone pedicle */
};

