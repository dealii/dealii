#####
##
## Copyright (C) 2013 by the deal.II authors
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
# This macro toggles the preference for static/shared libraries if
# DEAL_II_PREFER_STATIC_LIBS=TRUE but the final executable will still be
# dynamically linked, i.e. DEAL_II_STATIC_EXECUTABLE=OFF
#
# Usage:
#     SWITCH_LIBRARY_PREFERENCE()
#

MACRO(SWITCH_LIBRARY_PREFERENCE)
  IF(DEAL_II_PREFER_STATIC_LIBS AND NOT DEAL_II_STATIC_EXECUTABLE)
    #
    # Invert the search order for libraries when DEAL_II_PREFER_STATIC_LIBS
    # is set. This will prefer static archives instead of shared libraries:
    LIST(REVERSE CMAKE_FIND_LIBRARY_SUFFIXES)
  ENDIF()
ENDMACRO()
