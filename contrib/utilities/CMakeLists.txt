## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# At the present time we only install one script from this directory.
#

set(_scripts
  ${CMAKE_CURRENT_SOURCE_DIR}/indent.py
  )
file(COPY ${_scripts}
  DESTINATION ${CMAKE_BINARY_DIR}/${DEAL_II_SHARE_RELDIR}/scripts
  )
install(FILES ${_scripts}
  DESTINATION ${DEAL_II_SHARE_RELDIR}/scripts
  COMPONENT library
  )
