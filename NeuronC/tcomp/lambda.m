;
func lambda(xdia,xrm) {

 rrm = xrm * 1e8 / (3.14159 * xdia);
 rri = 200 * 1e4 / (3.14159 * xdia * xdia / 4); 
 return (sqrt(rrm/rri));
};


