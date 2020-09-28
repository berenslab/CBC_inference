;
RI = 200;

/* function lambda returns space constant in microns,
    given diameter of dendrite in microns (xdia)
    and Rm (xrm) in Ohm-cm2.
*/

func lambda(xdia,xrm) {

 rrm = xrm * 1e8 / (3.14159 * xdia);
 rri = RI * 1e4 / (3.14159 * xdia * xdia / 4); 
 return (sqrt(rrm/rri));
};

/* function rri returns ri (ohms per micron),
    given diameter of dendrite in microns (xdia)
*/

func xrri(xdia) {
 xri = RI * 1e4 / (3.14159 * xdia * xdia / 4); 
 return (xri);
};

func xrrg(gjsize,rg) {
  xrg = 1 / (gjsize / rg );
  return (xrg);
};


