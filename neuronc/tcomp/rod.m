/* rod function */

func mrod (xpos,ypos,n) { 

/* rodrm = 5000; */
at   [n] rod (xpos,ypos); 

/* conn [n] to [n,1]   cable dia  2 length 20 ri 2000 rm 1e5; /* */

at   [n]          sphere dia  7       rm 3000;	/* = area of 2 x 25 um rod os */
/*conn [n] to [n][1] cable  dia .1 length .2 rm rodrm;  /* connecting cilium */
at   [n]          sphere dia  3           rm rodrm;
conn [n] to [n][1] cable  dia .2 length 50 rm rodrm;
at   [n][1]        sphere dia  3           rm rodrm;

/* conn [] to[] synapse expon 6 igain 7.5 thresh -.049 maxcon 1e-10 delay 1.5 */

return (1);		/* return minor node number for rod spherule */
};

