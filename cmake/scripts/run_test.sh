## ---------------------------------------------------------------------
##
## Copyright (C) 2015 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Helper script used in testsuite targets to run ("run") and compare
# ("diff") individual tests.
#
# Usage:
#   run_test.sh run TEST_FULL RUN_COMMAND
#
#   run_test.sh diff TEST_FULL DIFF_EXECUTABLE DIFF_EXECUTABLE COMPARISON_FILE
#

set -u

STAGE="$1"
TEST_FULL="$2"
RUN_COMMAND="$3"
NUMDIFF_EXECUTABLE="$4"
DIFF_EXECUTABLE="$5"
COMPARISON_FILE="$6"

# Ensure uniform sorting for pathname expansion
export LC_ALL=C

#
# Add a top level target to run and compare the test:
#

run(){
  rm -f failing_output
  rm -f output
  rm -f stdout

  ${RUN_COMMAND} > stdout 2>&1
  RETURN_VALUE=$?

  [ -f output ] || mv stdout output

  if [ $RETURN_VALUE -ne 0 ]; then
    mv output failing_output
    echo "${TEST_FULL}: BUILD successful."
    echo "${TEST_FULL}: RUN failed. ------ Return code $RETURN_VALUE"
    echo "${TEST_FULL}: RUN failed. ------ Result: `pwd`/failing_output"
    echo "${TEST_FULL}: RUN failed. ------ Partial output:"
    cat failing_output
    if [ -f stdout ]; then
      echo ""
      echo "${TEST_FULL}: RUN failed. ------ Additional output on stdout/stderr:"
      echo ""
      cat stdout
    fi
    exit 1
  fi
}

diff() {
  rm -f failing_diff*
  rm -f diff*
  touch diff

  test_successful=false

  #
  # Pick up main comparison file and all variants. A valid variant name is
  # of the form [...].output.[STRING]
  #
  for file in "${COMPARISON_FILE}"*; do
    # determine variant name (empty string for main comparison file):
    variant=${file#*.output}

    #
    # Run diff or numdiff (if available) to determine whether files are the
    # same. Create a diff file "diff${variant}" for each variant file that
    # is found (including the main comparison file).
    #
    case ${NUMDIFF_EXECUTABLE} in
      *numdiff*)
        ${NUMDIFF_EXECUTABLE} -a 1e-6 -r 1e-8 -s ' \t\n:<>=,;' \
                              "${file}" output > diff${variant}
        ;;
      *)
        "${DIFF_EXECUTABLE}" "${file}" output > diff${variant}
        ;;
    esac

    if [ $? -eq 0 ]; then
      #
      # Ensure that only a single diff file with no contents remains
      # (numdiff has the bad habit of being very verbose...):
      #
      rm -f diff*
      touch diff

      if [ -n "${variant}" ]; then
        #
        # In case of a successful comparison against a variant, store the
        # fact that we compared against a variant in the diff file.
        #
        echo "${TEST_FULL}: DIFF successful. - Variant: ${file}" > diff
      fi

      test_successful=true
      break
    fi
  done

  #
  # If none of the diffs succeeded, use the diff against the main comparison
  # file. Output the first few lines of the output of numdiff, followed by
  # the results of regular diff since the latter is just more readable.
  #
  if [ $test_successful = false ] ; then
    for file in diff*; do
      mv "$file" failing_"$file"
    done
    echo "${TEST_FULL}: BUILD successful."
    echo "${TEST_FULL}: RUN successful."
    echo "${TEST_FULL}: DIFF failed. ------ Source: ${COMPARISON_FILE}"
    echo "${TEST_FULL}: DIFF failed. ------ Result: `pwd`/output"
    echo "Check `pwd`/output ${COMPARISON_FILE}"
    echo "${TEST_FULL}: DIFF failed. ------ Diff:   `pwd`/failing_diff"
    echo "${TEST_FULL}: DIFF failed. ------ First 8 lines of numdiff/diff output:"
    cat failing_diff | head -n 8
    echo "${TEST_FULL}: DIFF failed. ------ First 50 lines diff output:"
    "${DIFF_EXECUTABLE}" -c "${COMPARISON_FILE}" output | head -n 50
    exit 1
  fi
  exit 0
}

case $STAGE in
  run)
    run;;
  diff)
    diff;;
  *)
    exit 1;;
esac
