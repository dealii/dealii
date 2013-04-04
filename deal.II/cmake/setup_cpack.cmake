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
#
###########################################################################
#                                                                         #
#                              Setup cpack:                               #
#                                                                         #
###########################################################################

#
# General setup:
#
SET(CPACK_PACKAGE_NAME "${DEAL_II_PACKAGE_NAME}")
SET(CPACK_PACKAGE_VENDOR "${DEAL_II_PACKAGE_VENDOR}")
SET(CPACK_PACKAGE_VERSION "${DEAL_II_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_VERSION_MAJOR "${DEAL_II_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${DEAL_II_VERSION_MINOR}")

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${DEAL_II_PACKAGE_DESCRIPTION}")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")

SET(CPACK_SOURCE_GENERATOR TGZ ZIP) # TBZ2
SET(CPACK_SOURCE_IGNORE_FILES "/\\\\.svn/;\\\\.swp$;.*~$")

SET(CPACK_GENERATOR TGZ ZIP) # TBZ2

#
# Generator specific configuration:
#

SET(CPACK_ARCHIVE_COMPONENT_INSTALL TRUE)

# TODO

#INCLUDE(InstallRequiredSystemLibraries)
INCLUDE(CPack)
