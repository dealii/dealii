#!/bin/sh
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2018 - 2024 by the deal.II authors
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
# or 6m30s. This can be very substantially accelerated by cloning the
# repository in a directory such as /tmp.
#
# NOTE: The script has existed in this repository only for a finite time.
#       So, if you go back in history, as this script does, it will
#       eventually disappear. As a consequence, you cannot execute it
#       in its current place -- copy it to the top level of the working
#       copy you are working in, where it will be an untracked file that
#       git leaves alone as it cycles through repository commits.
#
#       You will then want to run the script using a command such as
#         ./count_lines.sh | tee ../line-count.dat
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
# displaying the number of lines of code in .h and .cc files, but excluding
# files in contributed libraries. The resulting output can be piped into a data
# file and then be visualized by importing into a spreadsheet, or using the
# following GNUPLOT script:
#
#  set xdata time
#  set timefmt "%Y-%m-%d"
#  set xrange ["1997-11-01":"2023-12-31"]
#  set yrange [0:1000000]
#  set format x "%Y"
#  set style data lines
#  set key top left
#  
#  unset arrow
#  unset label
#  
#  set arrow from "2023-07-07",500000 to "2023-07-07",800000 head
#  set label "9.5" at "2023-07-07",475000
#  
#  set arrow from "2022-06-24",450000 to "2022-06-24",750000 head
#  set label "9.4" at "2022-06-24",425000
#  
#  set arrow from "2021-06-18",400000 to "2021-06-18",700000 head
#  set label "9.3" at "2021-06-18",375000
#  
#  set arrow from "2020-05-20",350000 to "2020-05-20",650000 head
#  set label "9.2" at "2020-05-20",325000
#  
#  set arrow from "2019-05-21",300000 to "2019-05-21",600000 head
#  set label "9.1" at "2019-05-21",275000
#  
#  set arrow from "2018-05-11",250000 to "2018-05-11",550000 head lw 2
#  set label "9.0" at "2018-05-11",225000
#  
#  set arrow from "2013-07-24",150000 to "2013-07-24",450000 head lw 2
#  set label "8.0" at "2013-07-24",125000
#  
#  set arrow from "2011-01-09",100000 to "2011-01-09",400000 head lw 2
#  set label "7.0" at "2011-01-09",75000
#  
#  
#  set terminal png
#  set output "line-count.png"
#      
#  plot "< cat line-count.dat | sort" using 1:2 title "Lines of code in source files", \
#       "" using 1:3 title "Lines of code in tests", \
#       "< echo 2018-01-01 700000 ; echo 2023-12-31 1000000" using 1:2 title "50,000 lines per year"
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

  files_source=$(find . -name '*.h' -or -name '*.cc' -type f | \
                 grep -E -i -v '(tests|boost|umfpack|bundled)/')
  lines_source=$(cat $files_source | wc -l)

  files_tests=$(find . -name '*.h' -or -name '*.cc' -type f | \
                grep -E -i -v '(boost|umfpack|bundled)/' | \
                grep tests/)
  lines_tests=$(cat $files_tests | wc -l)
  
  echo "$date" "$lines_source" "$lines_tests"
done

git checkout -q ${current_branch}
