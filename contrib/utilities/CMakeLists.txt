## ---------------------------------------------------------------------
##
## Copyright (C) 2018 - 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# At the present time we only install one script from this directory.
#

SET(_scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/indent.py
  )
FILE(COPY ${_scripts}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts
  )
INSTALL(FILES ${_scripts}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/scripts
  COMPONENT library
  )
