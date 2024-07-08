## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2015 - 2024 by the deal.II authors
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
# Helper script used in testsuite targets to run ("run") and compare
# ("diff") individual tests.
#
# Usage:
#   run_test.sh run TEST_FULL [COMMAND and ARGS]
#
#   run_test.sh diff TEST_FULL NUMDIFF_EXECUTABLE COMPARISON_FILE
#

set -u

STAGE="$1"
TEST_FULL="$2"
shift 2

# Ensure uniform sorting for pathname expansion
export LC_ALL=C

case $STAGE in
  run)
    ##
    # run stage:
    #   - run specified command and parameters given by $@
    #   - creates file "output" on success containing stdout and stderr of
    #     the test
    #   - if test exits with non-zero return value, output is renamed to
    #     failing_output
    ##

    #
    # If TEST_N_THREADS is not equal to zero then:
    #  - Export the environment variable DEAL_II_NUM_THREADS set to
    #    $TEST_N_THREADS. This will enforce an upper bound of
    #    DEAL_II_NUM_THREADS during thread initialization of the threading
    #    pool in deal.II.
    #  - Export TEST_N_THREADS which is internally used in the deal.II
    #    testsuite to explicitly set the number of threads in the header
    #    file tests.h.
    #

    # Read in TEST_N_THREADS from command line:
    if [[ $1 == TEST_N_THREADS=* ]]; then
      TEST_N_THREADS="${1#TEST_N_THREADS=}"
      shift
    fi

    if [ "${TEST_N_THREADS:-0}" -ne 0 ]; then
      export DEAL_II_NUM_THREADS="${TEST_N_THREADS}"
      export TEST_N_THREADS
    fi

    # Limit the OpenMP pool to two threads. Set both variables in case the
    # caller has either one set.
    #
    # These variables do different things and are interpreted by, e.g., GOMP and
    # openBLAS in slightly different ways so it is best to set both to avoid
    # inconsistencies. For reference:
    # 1. OMP_NUM_THREADS permits nesting, e.g., if we were using OpenMP we could
    #    set it to 4,2,1 to set parallelization levels inside parallelized blocks
    # 2. OMP_THREAD_LIMIT limits the total number of threads, independent of
    #    nesting
    #
    # 3. In addition we set the OMP_PROC_BIND variable to false to allow
    #    free movement of the two worker threads and silence a KOKKOS
    #    warning (which might insist on this variable to be defined).
    export OMP_NUM_THREADS="2"
    export OMP_THREAD_LIMIT="2"
    export OMP_PROC_BIND="false"

    #
    # OpenMPI parameters:
    #  - Allow to oversubsribe the system, meaning that we can issue tests
    #    with more mpi ranks than available cores.
    #  - Ensure that we do not bind mpi tests to specific
    #    cores/processors/sockets. Otherwise we run the risk that multiple
    #    mpi tests (for example with two ranks) are all pinned to the same
    #    processor core, bringing everything to a grinding halt.
    #
    #    If the test runs exclusively, however, do not override MPI binding
    #    policies. This is important to ensure that performance tests are
    #    scheduled properly.
    #
    # for OpenMPI 4 and older:
    export OMPI_MCA_rmaps_base_oversubscribe=1
    # for OpenMPI 5:
    export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
    if [ -z "${TEST_IS_EXCLUSIVE+x}" ]; then
      # for OpenMPI 4 and older:
      export OMPI_MCA_hwloc_base_binding_policy=none
      # for OpenMPI 5:
      export PRTE_MCA_hwloc_default_binding_policy=none
    fi

    #
    # Kokkos parameters:
    #  - Disable Kokkos runtime warnings so that we can oversubscribe
    #    threads without Kokkos complaining. This should also help with
    #    some spurious warnings that we get in tests depending on who
    #    initialized openmp first, kokkos or another external dependency.
    export KOKKOS_DISABLE_WARNINGS=1

    rm -f failing_output
    rm -f output
    rm -f stdout

    "$@" > stdout 2>&1
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

    ;;

  diff)

    ##
    # diff stage:
    #   - compares the file "output" against various comparison files.
    ##

    NUMDIFF_EXECUTABLE="$1"
    COMPARISON_FILE="$2"

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
      variant="${file#*.output}"

      #
      # Configure numdiff with
      #   - absolute precision set to 1e-6, everything below is
      #     regarded as being equal to 0
      #   - relative differences of 1e-8
      #   - [space][tab][newline]=,:;<>[](){}^ as separators between
      #     numbers
      #
      "${NUMDIFF_EXECUTABLE}" -a 1e-6 -r 1e-8 -s ' \t\r\n=,:;<>[](){}^' \
                              "${file}" output > diff${variant}

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
    # file. Output the first few lines of the output of numdiff.
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
      echo "${TEST_FULL}: DIFF failed. ------ First 20 lines of numdiff output:"
      head -n 20 failing_diff
      exit 1
    fi
    exit 0

    ;;
  *)

    exit 1
    ;;
esac
