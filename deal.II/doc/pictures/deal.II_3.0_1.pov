 
#include "colors.inc"

camera {
  location <-5,.15,-2>
  look_at <1,.5,1>
  angle 45
}

light_source { <500,500,-1000> White }

plane { y,0
  texture {
    pigment { color red 0 green 0 blue 1 }
    finish { reflection .35 specular 1 }
    normal { ripples .35 turbulence .5 scale .25 }
  }
}

text { ttf "timrom.ttf" "deal.II 3.3" .25, 0
  pigment { BrightGold }
  finish { specular 1 }
  translate -3*x
  translate 0.02*y
}

#include "skies.inc"

sky_sphere { S_Cloud2 }

