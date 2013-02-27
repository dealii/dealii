#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
## Author: Matthias Maier <matthias.maier@iwr.uni-heidelberg.de>
##
#####

#
# Configuration for the ARPACK library:
#


MACRO(FEATURE_ARPACK_CONFIGURE_EXTERNAL var)
  LIST(APPEND DEAL_II_EXTERNAL_LIBRARIES ${ARPACK_LIBRARIES})
  ADD_FLAGS(CMAKE_SHARED_LINKER_FLAGS "${ARPACK_LINKER_FLAGS}")

  SET(${var} TRUE)
ENDMACRO()


CONFIGURE_FEATURE(ARPACK)
