#include "colors.inc" 
#include "textures.inc" 


camera {
  location <2,8,-14>
  look_at <2,0,2>
  angle 30
}

light_source {
  <1,4,-7>
  color Grey
}

light_source {
  <0,20,0>
  color White
}
#declare HMIN=0;
#declare HMAX=1;

#declare Tex=texture{
  pigment {
    gradient y
    scale y*(HMAX-HMIN)*0.1
    color_map {
      [0.00 color Light_Purple] 
      [0.95 color Light_Purple] 
      [1.00 color White]    
} } }

