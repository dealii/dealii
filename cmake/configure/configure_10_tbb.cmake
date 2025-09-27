## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2020 - 2022 by the deal.II authors
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
# Configuration for thread support in deal.II with the help of the tbb
# library:
#

macro(feature_tbb_find_external var)
  find_package(DEAL_II_TBB)

  if(TBB_FOUND)
    set(${var} TRUE)
  endif()

  set(DEAL_II_TBB_WITH_ONEAPI ${TBB_WITH_ONEAPI})

  #
  # TBB currently uses the version numbering scheme
  #
  #     YYYY.X
  #
  # (e.g., 2018.0) where YYYY is the year of the release and X is the yearly
  # release number. Older versions use
  #
  #     X.Y.Z
  #
  # (e.g., 4.2.1). Since we are compatible with all versions that use the new
  # numbering scheme we only check for very old versions here.
  #
  # TBB versions before 4.2 are missing some explicit calls to std::atomic::load
  # in ternary expressions; these cause compilation errors in some compilers
  # (such as GCC 8.1 and newer). To fix this we simply disallow all older
  # versions:
  #
  if(TBB_VERSION VERSION_LESS "4.2")
    # Clear the previously determined version numbers to avoid confusion
    set(TBB_VERSION "bundled")
    set(TBB_VERSION_MAJOR "")
    set(TBB_VERSION_MINOR "")

    message(STATUS
      "The externally provided TBB library is older than version 4.2.0, which "
      "cannot be used with deal.II."
      )
    set(TBB_ADDITIONAL_ERROR_STRING
      "The externally provided TBB library is older than version\n"
      "4.2.0, which is the oldest version compatible with deal.II and its\n"
      "supported compilers."
      )
    set(${var} FALSE)
  endif()
endmacro()


configure_feature(TBB)
