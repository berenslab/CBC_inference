/* nc.inc */

/* Generates scene and objects for nc */
/*    output file to "POV" ray-tracer. */

/* Read "povray.doc" to understand the format */

#include "shapes.inc"
#include "colors.inc"
#include "textures.inc"
#include "chars.inc"

global_settings { assumed_gamma 1.0 }

/*--------------------------------------*/

   /* two light sources, one point, one diffuse */

  light_source {<-500, -100, 300> color Gray80}
  light_source {<-500, 5000, 100> color Gray80}

  light_source {
   <300, -100, 300> color White
  
     area_light <200, 0, 0>, <0, 200, 0>, 4, 4
     adaptive 1
   }

//  camera {
//     location <0,  -20  10>		/* viewpoint location */
//     direction <0, 2, 0>		/* length of zoom lens */
//     up    <0, 1, 0>			/* aspect ratio */
//     right <1, 0, 0>			/* aspect ratio */
//     look_at <0, 0, 0>			/* viewpoint direction */
//  }
  
  default {				/* default finish */
     finish {
        phong 1 phong_size 100
        ambient .15
        diffuse .7
	reflection .08
     }
  }

//  plane {				/* add ground plane */
//    <0, 0, 1>, -10
//    pigment { Gray40 }
//  }
// 			/* put in some fog to make bkgnd light */ 
//  fog {	
//   color White
//   distance  100000.0
//  }
 			/* put in some fog to make bkgnd light */ 
  fog {	
   color Gray80
   distance  5000000.0
  }


/*--------------------------------------*/

/* nc-generated objects: */

#declare
  nc_sphere =
    sphere {
      <0,0,0>,1
      pigment {CadetBlue}
  }

#declare
  nc_cable =
    cylinder {
        <0,0,0>,       // Center of one end
        <1,0,0>,       // Center of other end
        1              // Radius
        pigment {Aquamarine  }
  }

#declare
  nc_cone =
    cone {
      <0,0,10>,0.1
      <0,0,1>,0.5
      pigment {CadetBlue}
  }

#declare
  nc_rod =
    cone {
      <0,0,10>,0.5
      <0,0,1>,0.5
      pigment {CadetBlue}
  }

#declare
  nc_electrod =
    cone {
        <0,0,0>, 0.2       // Center, radius of one end
        <0.5,20,0>, 2     // Center, radius of other end
        open
        pigment { Yellow filter 0.1 }
//        finish { refraction 1 ior 1.5 }
  }

/* Use predefined objects like this: */

//object {
//  nc_cone
//  pigment {Cadetblue}
//  scale 1.0
//  rotate <0,0,0>
//  translate <0,0,0>
//}

