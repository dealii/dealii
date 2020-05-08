#!/bin/sh
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
# This is a little script that can be used to count the lines of source
# code in both the regular source/include/example files as well as the
# tests directories. This is occasionally useful to assess the growth
# of the library over time using a crude metric of its size.
#
# The script assumes that it is called in the top-level directory
# of a deal.II git repository. It counts the number of code lines
# for every 20th commit going back from the current 'master' branch
# to the beginning of the repository history.
#
# It takes, on average, about 1/5th of a second to count lines for one
# commit, averaged over the lifetime of the deal.II project. At the time
# when this script was written, the repository had about 41,000 commits,
# so counting lines for every 20th commit takes approximately 400 seconds,
# or 6m30s.
#
# NOTE: The script has existed in this repository only for a finite time.
#       So, if you go back in history, as this script does, it will
#       eventually disappear. As a consequence, you cannot execute it
#       in its current place -- copy it to the top level of the working
#       copy you are working in, where it will be an untracked file that
#       git leaves alone as it cycles through repository commits.
#
# NOTE: Separately, since the script requires git to check out many
#       different versions, all of which have files that may not exist
#       in the current tip of the git repository you have on your hard
#       drive, or have files that have changed over time, it is going
#       to lead to certain heartbreak if you run this script in a directory
#       that has anything other than a pristine clone of the current git
#       repository. Do not run it in directory in which you do or have
#       done development work.
#
# The output of this file consists of three numbers per line, showing
#   date source-lines test-lines
# displaying the number of lines of code in .h and .cc files (and .cu and
# .cuh files for CUDA), but excluding files in contributed libraries. The
# resulting output can be piped into a data file and then be visualized
# by importing into a spreadsheet, or using the following GNUPLOT script:
#
#    set xdata time
#    set timefmt "%Y-%m-%d"
#    set xrange ["1997-11-01":"2018-08-01"]
#    set format x "%Y"
#    set style data lines
#    set key top left
#    
#    set terminal png
#    set output "line-count.png"
#    
#    plot "< cat line-count.dat | sort" using 1:2 title "Lines of code in source files", \
#         "" using 1:3 title "Lines of code in tests"
#
#

current_branch=$(git rev-parse --abbrev-ref HEAD)
git checkout -q master

commits=$(git log | \
          grep -E '^commit ' | \
          perl -p -e 's/^commit //g;' | \
          perl -e '$i=0; while (<>) { ++$i; if ($i % 200 == 0) { print; } }')

for commit in $commits ; do
  git checkout -q "$commit"
  
  date=$(git log --date=short --format="%ad" -n 1)

  files_source=$(find . -name '*.h' -or -name '*.cc' -or -name '*.cu' -or -name '*.cuh' -type f | \
                 grep -E -i -v '(tests|boost|umfpack|bundled)/')
  lines_source=$(cat $files_source | wc -l)

  files_tests=$(find . -name '*.h' -or -name '*.cc' -type f | \
                grep -E -i -v '(boost|umfpack|bundled)/' | \
                grep tests/)
  lines_tests=$(cat $files_tests | wc -l)
  
  echo "$date" "$lines_source" "$lines_tests"
done

git checkout -q ${current_branch}

