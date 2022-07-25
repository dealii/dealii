## ---------------------------------------------------------------------
##
## Copyright (C) 2013 by the deal.II authors
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
# Dashboard configuration:
#

SET(CTEST_PROJECT_NAME "deal.II")

SET(CTEST_DROP_METHOD "https")
SET(CTEST_DROP_SITE "cdash.dealii.org")
SET(CTEST_DROP_LOCATION "/submit.php?project=deal.II")
SET(CTEST_DROP_SITE_CDASH TRUE)

#
# We use the CTEST_UPDATE routine to query information about the repository
# but we don't want to actually perform an update. Thus, replace the update
# by a noop:
#
SET(CTEST_GIT_COMMAND "git")
SET(CTEST_UPDATE_VERSION_ONLY true)

SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS   100)
SET(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 300)

# number of lines to submit before an error:
SET(CTEST_CUSTOM_ERROR_PRE_CONTEXT            5)
# number of lines to submit after an error:
SET(CTEST_CUSTOM_ERROR_POST_CONTEXT          20)

#
# Coverage options:
#

SET(CTEST_EXTRA_COVERAGE_GLOB
  # These files should have executable lines and therefore coverage:
  # source/**/*.cc
  )

SET(CTEST_CUSTOM_COVERAGE_EXCLUDE
  "/bundled"
  "/cmake/scripts/"
  "/contrib"
  "/examples"
  "/tests"
  )
