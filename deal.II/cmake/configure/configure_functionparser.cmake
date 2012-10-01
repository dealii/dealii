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
# Configuration for functionparser
#

OPTION(DEAL_II_WITH_FUNCTIONPARSER
  "Build deal.II with support for functionparser."
  OFF)


MACRO(FEATURE_FUNCTIONPARSER_FIND_EXTERNAL var)
  MESSAGE(STATUS
    "No module available for finding functionparser externally."
    )
ENDMACRO()


MACRO(FEATURE_FUNCTIONPARSER_CONFIGURE_EXTERNAL var)
ENDMACRO()


MACRO(FEATURE_FUNCTIONPARSER_CONFIGURE_BUNDLED var)
  INCLUDE_DIRECTORIES(${FUNCTIONPARSER_FOLDER})
  SET(HAVE_FUNCTIONPARSER TRUE)
  SET(${var} TRUE)
ENDMACRO()


SET(FEATURE_FUNCTIONPARSER_CUSTOM_ERROR_MESSAGE TRUE)


MACRO(FEATURE_FUNCTIONPARSER_ERROR_MESSAGE)
  MESSAGE(SEND_ERROR "\n"
    "No module available for finding functionparser externally.\n"
    "Disable DEAL_II_WITH_FUNCTIONPARSER, or enable DEAL_II_ALLOW_BUNDLED.\n\n"
    )
ENDMACRO()


CONFIGURE_FEATURE(FUNCTIONPARSER)

