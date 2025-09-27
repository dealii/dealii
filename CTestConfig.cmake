## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 by the deal.II authors
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
# Dashboard configuration:
#

set(CTEST_PROJECT_NAME "deal.II")

set(CTEST_DROP_METHOD "https")
set(CTEST_DROP_SITE "cdash.dealii.org")
set(CTEST_DROP_LOCATION "/submit.php?project=deal.II")
set(CTEST_DROP_SITE_CDASH TRUE)

#
# We use the CTEST_UPDATE routine to query information about the repository
# but we don't want to actually perform an update. Thus, replace the update
# by a noop:
#
set(CTEST_GIT_COMMAND "git")
set(CTEST_UPDATE_VERSION_ONLY true)

set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   100)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 300)

# number of lines to submit before an error:
set(CTEST_CUSTOM_ERROR_PRE_CONTEXT            5)
# number of lines to submit after an error:
set(CTEST_CUSTOM_ERROR_POST_CONTEXT          20)

#
# Coverage options:
#

set(CTEST_EXTRA_COVERAGE_GLOB
  # These files should have executable lines and therefore coverage:
  # source/**/*.cc
  )

set(CTEST_CUSTOM_COVERAGE_EXCLUDE
  "/bundled"
  "/cmake/scripts/"
  "/contrib"
  "/examples"
  "/tests"
  )
