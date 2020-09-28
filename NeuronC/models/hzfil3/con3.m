/* cone function */

func mcone (xpos,ypos,n) { 

/* conerm = 3000; */
at   [n] cone (xpos,ypos) attf 1.5; 

/*conn [n] to [n][1] cable dia 10 length 10 ri 5000 rm 1e4;/* 5 x rod os area*/ 

at   [n]            sphere dia  16           rm 500000 vrest conrest; 
/*conn [n] to [n+1] cable  dia  .1 length .2 rm conerm;  /* connecting cilium */
at   [n]            sphere dia   3           rm conerm vrest conrest;
conn [n] to [n][1]  cable  dia   1 length 50 rm conerm vrest conrest;
at   [n][1]         sphere dia   5           rm conerm vrest conrestt;

return (1);	/* return minor node num for cone pedicle */
};

