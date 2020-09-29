/* ncray */

/* Rudimentary interface to generate symbolic output file for 
   "POV" ray-tracer.  Requires "nc.pov" include file to set up
   lights and default surface.  Camera location defines 
   the rotations from the "nc" file.
*/


extern "C" {

#include <stdio.h>
#ifdef CPML
#include <cpml.h>
#endif
#include <string.h>
#include <math.h>
#include <ctype.h>
}

#include "nc.h"
#include "ndef.h"
#include "ncsub.h"
#include "ncelem.h"
#include "colors.h"
#include "y.tab.h"
#include "ncio.h"

extern double ncmatc[4][4];

void rayconn (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	int ctype, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4]);
void transf (double x,double y,double z,
	      double *tx,double *ty,double *tz, double (*mat)[4]); 
void setray (double xrot, double yrot, double zrot, 
	double xcent, double ycent, double zcent,
	double xtran, double ytran, double ztran,
	double scale);
void rprintf (char *fmt,double x, double y, double z, double size, int color);

/*--------------------------------------*/

/* These colors duplicate the colors in "colors.h" as defined in
   mprintx.c and mprintc.c.
*/

const char *color_strings[]= {

/* standard VGA colors */

        "Black","Blue","Green","Cyan","Red","Magenta","Yellow*0.6","White",
        "Gray60","LightBlue","Aquamarine","CadetBlue","Firebrick",
        "VioletRed","Yellow", "VLightGrey",

/* for vcolor, lcolor, cacolor: blue - green - yellow - red */

	"color red .024 green 0    blue 0.2",
	"color red .082 green 0    blue .57",
	"color red .082 green 0    blue  1.",
	"color red 0    green .22  blue  1.",
	"color red 0    green .64  blue .9",
	"color red 0    green .93  blue .47",
	"color red 0    green  1.  blue  0",
	"color red .35  green  1.  blue  0",
	"color red .67  green  1.  blue  0",
	"color red .98  green .98  blue  0",
	"color red  1   green .75  blue  0",
	"color red  1   green .53  blue  0",
	"color red  1   green .3   blue  0",
	"color red  1   green .12  blue  0",
	"color red  1   green .016 blue .0078",
	"color red  1   green .0   blue .043",

        "color red  0    green  0  blue  0.63",                /* blue - purple - red */
        "color red 0.031 green  0  blue  0.86",
        "color red 0.047 green  0  blue  0.98",
        "color red 0.078 green  0  blue  0.94",
        "color red 0.16  green  0  blue  0.86",
        "color red 0.24  green  0  blue  0.78",
        "color red 0.31  green  0  blue  0.71",
        "color red 0.39  green  0  blue  0.63",
        "color red 0.47  green  0  blue  0.55",
        "color red 0.55  green  0  blue  0.47",
        "color red 0.63  green  0  blue  0.39",
        "color red 0.71  green  0  blue  0.31",
        "color red 0.78  green  0  blue  0.24",
        "color red 0.86  green  0  blue  0.16",
        "color red 0.94  green  0  blue  0.078",
        "color red 1     green  0  blue  0",


        "color red 0      green 0    blue   0",          /* black - red */
        "color red 0.031  green 0    blue  0",
        "color red 0.047  green 0    blue 0",
        "color red 0.078  green 0    blue 0",
        "color red 0.16   green 0    blue 0",
        "color red 0.24   green 0    blue 0",
        "color red 0.31   green 0    blue 0",
        "color red 0.39   green 0    blue 0",
        "color red 0.47   green 0    blue 0",
        "color red 0.55   green 0    blue 0",
        "color red 0.63   green 0    blue 0",
        "color red 0.71   green 0    blue 0",
        "color red 0.78   green 0    blue 0",
        "color red 0.86   green 0    blue 0",
        "color red 0.94   green 0    blue 0",
        "color red 1      green 0    blue 0",

        "color red 0    green 0      blue 0",          /* black - green */
        "color red 0    green 0.031  blue 0",
        "color red 0    green 0.047  blue 0",
        "color red 0    green 0.078  blue 0",
        "color red 0    green 0.16   blue 0",
        "color red 0    green 0.24   blue 0",
        "color red 0    green 0.31   blue 0",
        "color red 0    green 0.39   blue 0",
        "color red 0    green 0.47   blue 0",
        "color red 0    green 0.55   blue 0",
        "color red 0    green 0.63   blue 0",
        "color red 0    green 0.71   blue 0",
        "color red 0    green 0.78   blue 0",
        "color red 0    green 0.86   blue 0",
        "color red 0    green 0.94   blue 0",
        "color red 0    green 1      blue 0",

        "color red 0    green 0     blue 0",            /* black - blue */
        "color red 0    green 0     blue 0.031",
        "color red 0    green 0     blue 0.047",
        "color red 0    green 0     blue 0.078",
        "color red 0    green 0     blue 0.16",
        "color red 0    green 0     blue 0.24",
        "color red 0    green 0     blue 0.31",
        "color red 0    green 0     blue 0.39",
        "color red 0    green 0     blue 0.47",
        "color red 0    green 0     blue 0.55",
        "color red 0    green 0     blue 0.63",
        "color red 0    green 0     blue 0.71",
        "color red 0    green 0     blue 0.78",
        "color red 0    green 0     blue 0.86",
        "color red 0    green 0     blue 0.94",
        "color red 0    green 0     blue 1",

/* for B/W grayscale */

	"White*.00",
	"White*.01",
	"White*.02",
	"White*.03",
	"White*.04",
	"White*.05",
	"White*.06",
	"White*.07",
	"White*.08",
	"White*.09",
	"White*.10",
	"White*.11",
	"White*.12",
	"White*.13",
	"White*.14",
	"White*.15",
	"White*.16",
	"White*.17",
	"White*.18",
	"White*.19",
	"White*.20",
	"White*.21",
	"White*.22",
	"White*.23",
	"White*.24",
	"White*.25",
	"White*.26",
	"White*.27",
	"White*.28",
	"White*.29",
	"White*.30",
	"White*.31",
	"White*.32",
	"White*.33",
	"White*.34",
	"White*.35",
	"White*.36",
	"White*.37",
	"White*.38",
	"White*.39",
	"White*.40",
	"White*.41",
	"White*.42",
	"White*.43",
	"White*.44",
	"White*.45",
	"White*.46",
	"White*.47",
	"White*.48",
	"White*.49",
	"White*.50",
	"White*.51",
	"White*.52",
	"White*.53",
	"White*.54",
	"White*.55",
	"White*.56",
	"White*.57",
	"White*.58",
	"White*.59",
	"White*.60",
	"White*.61",
	"White*.62",
	"White*.63",
	"White*.64",
	"White*.65",
	"White*.66",
	"White*.67",
	"White*.68",
	"White*.69",
	"White*.70",
	"White*.71",
	"White*.72",
	"White*.73",
	"White*.74",
	"White*.75",
	"White*.76",
	"White*.77",
	"White*.78",
	"White*.79",
	"White*.80",
	"White*.81",
	"White*.82",
	"White*.83",
	"White*.84",
	"White*.85",
	"White*.86",
	"White*.87",
	"White*.88",
	"White*.89",
	"White*.90",
	"White*.91",
	"White*.92",
	"White*.93",
	"White*.94",
	"White*.95",
	"White*.96",
	"White*.97",
	"White*.98",
	"White*.99",

/* for a variety of colors */

/* taken from ccols.txt, converted by "convcol2", to make "ccols.x2" */
/* from /usr/lib/X11/rgb.txt */

	"color red 0.980 green 0.941 blue 0.902", /* linen */ 
	"color red 0.980 green 0.922 blue 0.843", /* AntiqueWhite */ 
	"color red 1.000 green 0.937 blue 0.835", /* PapayaWhip */ 
	"color red 1.000 green 0.922 blue 0.804", /* BlanchedAlmond */ 
	"color red 1.000 green 0.855 blue 0.725", /* PeachPuff */ 
	"color red 1.000 green 0.894 blue 0.710", /* moccasin */ 
	"color red 1.000 green 0.980 blue 0.804", /* LemonChiffon */ 
	"color red 1.000 green 0.961 blue 0.933", /* seashell */ 
	"color red 0.941 green 1.000 blue 0.941", /* honeydew */ 
	"color red 0.961 green 1.000 blue 0.980", /* MintCream */ 
	"color red 0.941 green 1.000 blue 1.000", /* azure */ 
	"color red 0.902 green 0.902 blue 0.980", /* lavender */ 
	"color red 1.000 green 0.941 blue 0.961", /* LavenderBlush */ 
	"color red 1.000 green 0.894 blue 0.882", /* MistyRose */ 
	"color red 0.098 green 0.098 blue 0.439", /* MidnightBlue */ 
	"color red 0.000 green 0.000 blue 0.502", /* NavyBlue */ 
	"color red 0.392 green 0.584 blue 0.929", /* CornflowerBlue */ 
	"color red 0.282 green 0.239 blue 0.545", /* DarkSlateBlue */ 
	"color red 0.416 green 0.353 blue 0.804", /* SlateBlue */ 
	"color red 0.518 green 0.439 blue 1.000", /* LightSlateBlue */ 
	"color red 0.000 green 0.000 blue 0.804", /* MediumBlue */ 
	"color red 0.255 green 0.412 blue 0.882", /* RoyalBlue */ 
	"color red 0.118 green 0.565 blue 1.000", /* DodgerBlue */ 
	"color red 0.000 green 0.749 blue 1.000", /* DeepSkyBlue */ 
	"color red 0.529 green 0.808 blue 0.922", /* SkyBlue */ 
	"color red 0.275 green 0.510 blue 0.706", /* SteelBlue */ 
	"color red 0.690 green 0.769 blue 0.871", /* LightSteelBlue */ 
	"color red 0.678 green 0.847 blue 0.902", /* LightBlue */ 
	"color red 0.686 green 0.933 blue 0.933", /* PaleTurquoise */ 
	"color red 0.000 green 0.808 blue 0.820", /* DarkTurquoise */ 
	"color red 0.282 green 0.820 blue 0.800", /* MediumTurquoise */ 
	"color red 0.251 green 0.878 blue 0.816", /* turquoise */ 
	"color red 0.373 green 0.620 blue 0.627", /* CadetBlue */ 
	"color red 0.400 green 0.804 blue 0.667", /* MediumAquamarine */ 
	"color red 0.498 green 1.000 blue 0.831", /* aquamarine */ 
	"color red 0.333 green 0.420 blue 0.184", /* DarkOliveGreen */ 
	"color red 0.561 green 0.737 blue 0.561", /* DarkSeaGreen */ 
	"color red 0.180 green 0.545 blue 0.341", /* SeaGreen */ 
	"color red 0.235 green 0.702 blue 0.443", /* MediumSeaGreen */ 
	"color red 0.125 green 0.698 blue 0.667", /* LightSeaGreen */ 
	"color red 0.596 green 0.984 blue 0.596", /* PaleGreen */ 
	"color red 0.000 green 1.000 blue 0.498", /* SpringGreen */ 
	"color red 0.486 green 0.988 blue 0.000", /* LawnGreen */ 
	"color red 0.498 green 1.000 blue 0.000", /* chartreuse */ 
	"color red 0.000 green 0.980 blue 0.604", /* MediumSpringGreen */ 
	"color red 0.678 green 1.000 blue 0.184", /* GreenYellow */ 
	"color red 0.196 green 0.804 blue 0.196", /* LimeGreen */ 
	"color red 0.604 green 0.804 blue 0.196", /* YellowGreen */ 
	"color red 0.133 green 0.545 blue 0.133", /* ForestGreen */ 
	"color red 0.420 green 0.557 blue 0.137", /* OliveDrab */ 
	"color red 0.741 green 0.718 blue 0.420", /* DarkKhaki */ 
	"color red 0.941 green 0.902 blue 0.549", /* khaki */ 
	"color red 0.933 green 0.910 blue 0.667", /* PaleGoldenrod */ 
	"color red 0.980 green 0.980 blue 0.824", /* LightGoldenrodYellow */ 
	"color red 1.000 green 1.000 blue 0.878", /* LightYellow */ 
	"color red 1.000 green 0.843 blue 0.000", /* gold */ 
	"color red 0.933 green 0.867 blue 0.510", /* LightGoldenrod */ 
	"color red 0.855 green 0.647 blue 0.125", /* goldenrod */ 
	"color red 0.722 green 0.525 blue 0.043", /* DarkGoldenrod */ 
	"color red 0.737 green 0.561 blue 0.561", /* RosyBrown */ 
	"color red 0.804 green 0.361 blue 0.361", /* IndianRed */ 
	"color red 0.545 green 0.271 blue 0.075", /* SaddleBrown */ 
	"color red 0.627 green 0.322 blue 0.176", /* sienna */ 
	"color red 0.804 green 0.522 blue 0.247", /* peru */ 
	"color red 0.871 green 0.722 blue 0.529", /* burlywood */ 
	"color red 0.961 green 0.961 blue 0.863", /* beige */ 
	"color red 0.961 green 0.871 blue 0.702", /* wheat */ 
	"color red 0.957 green 0.643 blue 0.376", /* SandyBrown */ 
	"color red 0.824 green 0.706 blue 0.549", /* tan */ 
	"color red 0.824 green 0.412 blue 0.118", /* chocolate */ 
	"color red 0.698 green 0.133 blue 0.133", /* firebrick */ 
	"color red 0.647 green 0.165 blue 0.165", /* brown */ 
	"color red 0.914 green 0.588 blue 0.478", /* DarkSalmon */ 
	"color red 0.980 green 0.502 blue 0.447", /* salmon */ 
	"color red 1.000 green 0.627 blue 0.478", /* LightSalmon */ 
	"color red 1.000 green 0.647 blue 0.000", /* orange */ 
	"color red 1.000 green 0.549 blue 0.000", /* DarkOrange */ 
	"color red 1.000 green 0.498 blue 0.314", /* coral */ 
	"color red 0.941 green 0.502 blue 0.502", /* LightCoral */ 
	"color red 1.000 green 0.388 blue 0.278", /* tomato */ 
	"color red 1.000 green 0.271 blue 0.000", /* OrangeRed */ 
	"color red 1.000 green 0.000 blue 0.000", /* red */ 
	"color red 1.000 green 0.412 blue 0.706", /* HotPink */ 
	"color red 1.000 green 0.078 blue 0.576", /* DeepPink */ 
	"color red 1.000 green 0.753 blue 0.796", /* pink */ 
	"color red 1.000 green 0.714 blue 0.757", /* LightPink */ 
	"color red 0.859 green 0.439 blue 0.576", /* PaleVioletRed */ 
	"color red 0.690 green 0.188 blue 0.376", /* maroon */ 
	"color red 0.780 green 0.082 blue 0.522", /* MediumVioletRed */ 
	"color red 0.816 green 0.125 blue 0.565", /* VioletRed */ 
	"color red 1.000 green 0.000 blue 1.000", /* magenta */ 
	"color red 0.933 green 0.510 blue 0.933", /* violet */ 
	"color red 0.867 green 0.627 blue 0.867", /* plum */ 
	"color red 0.855 green 0.439 blue 0.839", /* orchid */ 
	"color red 0.729 green 0.333 blue 0.827", /* MediumOrchid */ 
	"color red 0.600 green 0.196 blue 0.800", /* DarkOrchid */ 
	"color red 0.580 green 0.000 blue 0.827", /* DarkViolet */ 
	"color red 0.541 green 0.169 blue 0.886", /* BlueViolet */ 
	"color red 0.627 green 0.125 blue 0.941", /* purple */ 
	"color red 0.576 green 0.439 blue 0.859", /* MediumPurple */ 
};

const char *getcol (int color)

/* Return the char string color.  */

{
 return (color_strings[color]);
}

double xr,yr,zr;

/*--------------------------------------*/

void initray (double xrot, double yrot, double zrot, 
	double xcent, double ycent, double zcent, 
	double rxcent, double rycent, double rzcent, double scal)
{

  double xcam,ycam,zcam,aspect;
  double xcamr,ycamr,zcamr;
  double xzoom,yzoom,zzoom;
  static int runyet=0;

  if (!runyet) runyet = 1;
  else return; 
    
  aspect = 1.0;
  xcam = 0;			/* original camera position */
  ycam = 0;
  zcam = 10*scal;
  xzoom = 0; 			/* frame size */ 
  yzoom = 9.5;  
  zzoom = 0;  
  ncfprintf (stdout,"/* Scene code for POV ray-tracer */\n");
  ncfprintf (stdout,"/* generated by nc. */\n");
  ncfprintf (stdout,"\n");
  ncfprintf (stdout,"#version 3.7;\n\n");
  ncfprintf (stdout,"#include \"nc.pov\"\n\n");

  setray (xrot,yrot,zrot,xcent,ycent,zcent,rxcent,rycent,rzcent,1.0);

/*  transf (xcam,ycam,zcam, &xcamr,&ycamr,&zcamr,ncmatc); */

xr = xrot;
yr = yrot;
zr = zrot;

  ncfprintf (stdout,"  camera {\n");
  ncfprintf (stdout,"    location <%g,  %g,  %g>  /* viewpoint location */\n",
			xcam,ycam,zcam);
  ncfprintf (stdout,"     up      <%g, %g, %g>	/* aspect ratio */\n",
			0.0,1.0,0.0);
  ncfprintf (stdout,"     right   <%g, %g, %g>	/* aspect ratio */\n",
			1.0*aspect,0.0,0.0);
  ncfprintf (stdout,"     direction <%g, %g, %g>	/* length of viewing lens */\n",
			xzoom,yzoom,zzoom);
  ncfprintf (stdout,"     look_at <%g, %g, %g>	/* viewpoint direction */\n",
			0.0,0.0,0.0);
  ncfprintf (stdout,"  }\n");
  
}

/*--------------------------------------*/

void raycomp(double x, double y, double z, double dia)

/* draw a compartment */

{
  int color;

  if (z==NULLVAL) z = 0.0;
  color = BLUE;
  ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		x,y,z,dia/2.0,getcol(color));
}

/*--------------------------------------*/

void raynode(node *npnt, double dscale, int color, double (*mat)[4])

/* Display node numbers in 3D location. */
/* See "drnode() and dr_node()" in ncdisp.cc. */

{
    double x,y,z,tx,ty,tz,xoffs,charsiz;
    static char numbuf[20];
    double dispsize = 40.0;
    int d;

//  if (!color) color = WHITE; 
//  color = BLACK; 
  if (npnt==NULL) return;
  x = npnt->xloc;
  y = npnt->yloc;
  z = npnt->zloc;
  transf (x,y,z,&tx,&ty,&tz,mat);

  if (dscale<0) {
      dscale = -dscale;
      d = int(dscale);
      dscale = (dscale - d) * 10;
     if (dscale==0) dscale = 1.0;
     if (npnt->labeltext != NULL) {
         sprintf (numbuf,"%s",npnt->labeltext);
	 color = npnt->label;
     } else {
       switch (d) {
         case 1: sprintf (numbuf,"%d",npnt->nodenm1); break;
         case 2: sprintf (numbuf,"%d",npnt->nodenm2); break;
         case 3: sprintf (numbuf,"%d",npnt->nodenm3); break;
         case 4: sprintf (numbuf,"%d",npnt->nodenm4); break;
         case 5: sprintf (numbuf,"%d %d",npnt->nodenm1,npnt->nodenm2); break;
         case 6: sprintf (numbuf,"%d %d",npnt->nodenm2,npnt->nodenm3); break;
         case 7: sprintf (numbuf,"%d %d %d",npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
               break;
       }
       color = BLACK;
     }
   }
   else if (npnt->nodenm3 != NULLVAL)
      sprintf (numbuf,"%d %d %d",npnt->nodenm1,npnt->nodenm2,npnt->nodenm3);
   else if (npnt->nodenm2 != NULLVAL)
      sprintf (numbuf,"%d %d",npnt->nodenm1,npnt->nodenm2);
   else
      sprintf (numbuf,"%d",npnt->nodenm1);

  charsiz = .015 * dscale;
  // gcwidth (charsiz);
  xoffs = strlen(numbuf) / 2.0;		/* find center of number */
  xoffs *= charsiz * .7;
  //gmove (tx-xoffs, ty - charsiz*1.1);	/* label node numbers */
  //gtext (numbuf);
  //rprintf(numbuf,(tx-xoffs)*dispsize*0.5,(ty-charsiz*1.1)*dispsize*0.5,1,0,color);
  rprintf(numbuf,tx,ty,tz,dscale,color);
}

/*--------------------------------------*/

void makcyl(double x1, double y1, double z1,
       double x2, double y2, double z2, double dia, int color)

{
     static char end1[40], end2[40];

 sprintf (end1,"<%.5g,%.5g,%.5g>", x1,y1,z1);
 sprintf (end2,"<%.5g,%.5g,%.5g>", x2,y2,z2);
 if (strcmp(end1,end2) != 0) 
    ncfprintf (stdout,"cylinder { %s %s,%g pigment {%s} }\n\n",
	end1,end2,dia/2.0,getcol(color));
}

/*--------------------------------------*/

void makcone(double x1, double y1, double z1,
       double x2, double y2, double z2, double dia, double dia2, int color)

{
     static char end1[40], end2[40];

 sprintf (end1,"<%.5g,%.5g,%.5g>,%g", x1,y1,z1,dia/2.0);
 sprintf (end2,"<%.5g,%.5g,%.5g>,%g", x2,y2,z2,dia2/2.0);
 if (strcmp(end1,end2) != 0) 
    ncfprintf (stdout,"cone { %s %s pigment {%s} }\n\n",
	end1,end2,getcol(color));
}

/*--------------------------------------*/

void raycable (double x1, double y1, double z1, 
		double x2, double y2, double z2, 
		double dia, double dia2, 
		double dscale, double n1dia, double n2dia, 
		int color, int hide, double (*mat)[4])
                
/*
   Describe a cable element in 3D for ray-tracer.
   Include file sets up camera angle and default surfaces.
*/

{
  double tx1,ty1,tz1;
  double tx2,ty2,tz2;

  if (z1==NULLVAL) z1 = 0.0;
  if (z2==NULLVAL) z2 = 0.0;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  if (color==NULLVAL) color = GREEN;
  if (dscale<0) {
     dscale = -dscale;
     if (dia < dscale)   dia = dscale;
     if (dia2 < dscale) dia2 = dscale;
  }
  else { 
    dia *= dscale;
    dia2 *= dscale;
  }
  if (dia == dia2) {
   if (x1!=x2 || y1!=y2 || z1!=z2) 
    makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia,color);
  }
  else {
   if (x1!=x2 || y1!=y2 || z1!=z2)
    makcone(tx1,ty1,tz1,tx2,ty2,tz2,dia,dia2,color);
  }
}

/*--------------------------------------*/

void raysphere (double x, double y, double z, double dia, double dscale, 
	int color, int hide, double (*mat)[4],int fill)
{
    double tx,ty,tz;

  if (z==NULLVAL) z = 0.0;
  if (color==NULLVAL) color = BLUE;
  if (dscale<0) {
     dscale = -dscale;
     dia *= dscale;
     if (dia < 0.5) dia = 0.5;
  }
  else {
     dia *= dscale;
  }
  transf (x,y,z,&tx,&ty,&tz,mat);
  ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		tx,ty,tz,dia/2.0,getcol(color));
}

/*--------------------------------------*/

void raysynap (double x1, double y1, double z1, 
	double x2, double y2, double z2, synapse* epnt,
	double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4])
{
    double tx1,ty1,tz1, tx2,ty2,tz2;
    double dia;
    double cx,cy,cz;

  if (color==NULLVAL) color=RED;
  if (z1==NULLVAL) z1 = 0.0;
  if (z2==NULLVAL) z2 = 0.0;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  dia = 0.5;
  cx = (tx1+tx2)*0.5;
  cy = (ty1+ty2)*0.5;
  cz = (tz1+tz2)*0.5;
  if (x1!=x2 || y1!=y2 || z1!=z2) {
    makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia,color);
  }
  ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		cx,cy,cz,dia/2.0,getcol(color));
}

/*--------------------------------------*/

void raycconn (double x1, double y1, double z1, double x2, double y2, double z2)
{
   rayconn (x1,y1,z1,x2,y2,z2, AXIALRES, 1.0, 0.0, 0.0, YELLOW, 0, 0);
}

/*--------------------------------------*/

void rayconn (double x1, double y1, double z1, 
	double x2, double y2, double z2, 
	int ctype, double dscale, double n1dia, double n2dia, 
	int color, int hide, double (*mat)[4])
{
  double dia,cx,cy,cz;
  double tx1,ty1,tz1, tx2,ty2,tz2;

  if (z1==NULLVAL) z1 = 0.0;
  if (z2==NULLVAL) z2 = 0.0;
  if (color==NULLVAL) color = YELLOW;
  transf (x1,y1,z1,&tx1,&ty1,&tz1,mat);
  transf (x2,y2,z2,&tx2,&ty2,&tz2,mat);
  cx = (tx1+tx2)*0.5;
  cy = (ty1+ty2)*0.5;
  cz = (tz1+tz2)*0.5;

  switch (ctype){

  case GJ: 
    if (dscale>0) {
      dia = dscale/2;
      ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		cx,cy,cz,dia/2.0,getcol(color));
      makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia/4,color);
    }  
    else {
      dia = 1.0 * -dscale;
      ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		cx,cy,cz,dia/2.0,getcol(color));
    }
    break;

  case RESISTOR:
  case AXIALRES:
    dia = 2.0;
    makcyl(tx1,ty1,tz1,tx2,ty2,tz2,dia,color);
    break;
  } 
}

/*--------------------------------------*/

void rayphotrec (double x, double y, double z, 
		     photorec *epnt, double dscale,
		     double n1dia,double n2dia,
		     int color, int hide, double (*mat)[4])
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double dia,len;
    int ctype,pigm;
    static int color_trans[TOTREC] = {MAGENTA,RED,GREEN,BLUE};

  ctype = epnt->ctype;
  dia = epnt->dia*dscale;
  if (color==NULLVAL) {
      if ((pigm=(int)(epnt->pigm)) == NULLVAL)
         pigm = (ctype==CONE ? 1 : 0);
      pigm = limit(pigm,0,TOTREC);
      color= color_trans[pigm];
  }
  len = 1.0;				/* length of photorec in um */
  if (dia<=0) dia = 2.0;
  if (z==NULLVAL) z = 0.0;
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  switch (ctype) {

  case VTRANSDUCER:
  case ITRANSDUCER:
  case CHR:
  case CONE:
  ncfprintf (stdout,"object { nc_cone scale <%g,%g,1> rotate <%g,%g,%g> translate <%g,%g,%g> pigment {%s} }\n\n",
	 dia,dia,xr,yr,zr,tx1,ty1,tz1,getcol(color));
  break;

  case ROD:
  ncfprintf (stdout,"object { nc_rod scale <%g,%g,1> rotate <%g,%g,%g> translate <%g,%g,%g> pigment {%s} }\n\n",
	 dia,dia,xr,yr,zr,tx1,ty1,tz1,getcol(color));
  break;

 }  /* switch */
}

/*--------------------------------------*/

void rayelec (double x, double y, double z, 
		     electrode *epnt, double dscale,
		     int color, int hide, double (*mat)[4])
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double dia,len;
    int ctype;

  ctype = epnt->ctype;
  dia = epnt->dia*dscale;
  len = epnt->length;
  if (color==NULLVAL) {
      color= MAGENTA;
  }
  if (dia<=0) dia = 2.0;
  if (z==NULLVAL) z = 0.0;
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  switch (ctype) {

  case ELECTRODE:
  ncfprintf (stdout,"object { nc_electrod scale <%g,%g,1> rotate <%g,%g,%g> translate <%g,%g,%g> pigment {%s} }\n\n",
	 dia,dia,xr,yr,zr,tx1,ty1,tz1,getcol(color));
  break;

 }  /* switch */
}

/*--------------------------------------*/

void rayload (double x, double y, double z, int ctype, 
	double n1dia, double n2dia, int color, int hide, double (*mat)[4])
                              
/* see "drphotorec()" above for explanation of
   frame, origin() and rotation. */
 
{
    double tx1,ty1,tz1,tx2,ty2,tz2;
    double dia,len;

  len = 10; 
  transf (x,y,z,&tx1,&ty1,&tz1,mat);
  transf (x,y,z+len,&tx2,&ty2,&tz2,mat);
  if (color==NULLVAL) color = BLUE;
  dia = 0.5;
  ncfprintf (stdout,"sphere { <%g,%g,%g>, %g pigment {%s} }\n\n",
		tx1,ty1,tz1,dia/2.0,getcol(color));
}

/*--------------------------------------*/

void rprintf (char *fmt,double x, double y, double z, double size, int color)

/* print a row of chars into 3D ray-traced picture */

{
   int i,slen,ch;

 slen = strlen(fmt);
/* for (i=0; i<slen; i++) {
   ch = toupper(*(fmt+i));
   if (ch== ' ') continue;
   ncfprintf (stdout,"object {char_%c scale %g translate <%g,%g,%g> pigment {%s} }\n\n",
      ch, size, x+size*5*i,y,z,getcol(color));
 }
*/
   ncfprintf(stdout,"text { ttf \"timrom.ttf\" \"%s\" 1, 0 scale %g translate <%g,%g,%g> pigment {%s} }\n\n",
		   // fmt,1.0,x+size*8,y,z,getcol(color));
		   fmt,size,x,y,z,getcol(color));


}
/*--------------------------------------*/

void raycalib(double x, double y, double len, double size, int color)
{
    static char numbuf[20];
    double xb,yb,zb,z,width,ls;

  if (color==WHITE) color = BLACK;
  if (color==NULLVAL) color = CYAN;
  z = 1;
  x = (x-0.5) * size * .95;
  y = (y-0.5) * size * .95;
  xb = x - len;
  ls = 0.02 * size;			/* letter scale */
  width = len * 0.01;
  yb = y + width;
  zb = z + width;
  sprintf (numbuf,"%g um",len);  
  ncfprintf (stdout,"box { <%g,%g,%g> <%g,%g,%g> pigment {%s} }\n\n",
			x,y,z,xb,yb,zb,getcol(color));
  rprintf (numbuf, x-len*0.7, y-ls,z,ls,color);
}

/*--------------------------------------*/

