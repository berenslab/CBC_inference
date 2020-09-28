
func mid(siz)

/* function to return center element of array */

{
  if (int(siz/2)*2==siz) {
/*    print "even"; */
    m = (siz+1) * siz / 2;
  }
  else {
/*    print "odd"; */
    m = (siz*siz-1) / 2;
  };
return m;
};

func midrow(siz)

/* function to return start of middle row of array */

{
  return (siz*int(siz/2));
};

