 
#include "colors.inc"

camera {
  location <-5,.15,-2>
  look_at <1,0.5,1>
  angle 45
}

light_source { <500,500,-1000> White }

plane { y,0
  texture {
    pigment { color rgbt <0,0,1,0.8> }
    finish { reflection .35 specular 1 }
    normal { ripples .35 turbulence .5 scale .25 }
  }
}

plane { y, -10
  texture {
    pigment { Blue }
  }
}


text { ttf "timrom.ttf" "deal.II  3.3" .25, 0
  pigment { BrightGold }
  finish { specular 1 }
  translate -3*x
  translate 0.02*y
}

text { ttf "timrom.ttf" "deal.II  3.3" .15, 0
  pigment { BrightGold }
  finish { specular 1 }
  rotate 90*x
  rotate 90*y
  translate -2.25*x
  translate 2.5*z
  translate -0.25*y
}


#include "skies.inc"

sky_sphere { S_Cloud2 }

