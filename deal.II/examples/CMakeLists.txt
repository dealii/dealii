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

IF(DEAL_II_COMPONENT_EXAMPLES)
  INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${DEAL_II_EXAMPLES_RELDIR}
    COMPONENT examples
    FILES_MATCHING
    #
    # Exclude folder structures: doc, doxygen, CMakeFiles,...
    #
    PATTERN "CMakeFiles*" EXCLUDE
    PATTERN "doc*" EXCLUDE
    PATTERN "doxygen*" EXCLUDE
    PATTERN ".svn*" EXCLUDE
    #
    # Glob Includes:
    #
    PATTERN "*.cc"
    PATTERN "*.prm"
    PATTERN "*.inp"
    PATTERN "step*/CMakeLists.txt"
    #
    # Special files:
    #
    PATTERN "output.reference.dat" # step-39
    PATTERN "postprocess.pl"       # step-39
    PATTERN "obstacle_file.pbm"    # step-42
    PATTERN "untitled.geo"         # step-49
    PATTERN "untitled.msh"         # step-49
    )
ENDIF()
