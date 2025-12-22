#!/bin/sh
## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
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
# This is a little script that checks if the feature branch is linear,
# i.e., no merges of the parent branch into the current branch are present.
#
# The parent branch to check against can be supplied as a parameter to
# this script, or otherwise it will be set to 'master' by default.
#
# NOTE: This script will do nothing on the parent branch since the log
#       will always be empty.
#

branch="${1:-master}"

get_merge_commits_since_parent_branch () {
  echo "$(git log --merges --pretty=format:"%h" $branch..)"
}

get_commit_parents () {
  result=$(git rev-list --parents -n1 $1)
  # first result is the commit itself, which we omit
  echo "$(cut --delimiter=' '  --fields=1 --complement <<< $result)"
}

commit_is_in_parent_branch () {
  return $(git merge-base --is-ancestor $1 $branch)
}

readarray -t merge_hash_array < <(get_merge_commits_since_parent_branch)

if [[ -z "${merge_hash_array}" ]] ; then
    echo "No merge commits present at all, everything is good!"
    exit 0
fi

echo "Merge commits found, checking if they originate from $branch."

for hash in "${merge_hash_array[@]}"; do
  readarray -d' ' merge_parents < <(get_commit_parents $hash)

  for parent in "${merge_parents[@]}"; do
    if (commit_is_in_parent_branch $parent) ; then
      echo "There is a merge commit coming from $branch!"
      echo "The commit hash is $parent"
      exit 1
    fi
  done
done

echo "None do, everything is good!"
exit 0
