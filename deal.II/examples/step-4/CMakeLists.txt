CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.0 REQUIRED
  HINTS
    "../../"
    "../../../../../"
    $ENV{DEAL_DIR}
    $ENV{DEAL_II_DIR}
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(step-4)

ADD_EXECUTABLE(step-4 step-4.cc)
DEAL_II_SETUP_TARGET(step-4)

ADD_CUSTOM_TARGET(run
  COMMAND step-4
  )
ADD_CUSTOM_TARGET(distclean
  COMMAND make clean
  COMMAND rm -f *gmv *gnuplot *gpl *eps *pov *vtk *ucd *.d2
  )
