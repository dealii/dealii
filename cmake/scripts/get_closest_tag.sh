#!/bin/sh
#
# Find the TAG that
#   - has a common ancestry with current HEAD
#   - with shortest (positive) distance of the common ancestor to HEAD
#

if ! git --version | grep -q "git version 2"; then
  # This script requires version 2 or newer
  exit 1
fi

head="$(git rev-parse HEAD)"

tags="$(git tag --sort=-creatordate)"

min_distance=0
min_tag=""
for tag in $tags; do
  tag_commit="$(git rev-parse $tag^{commit})"

  if [ "$tag_commit" = "$head" ]; then
    echo "$tag"
    exit 0;
  fi

  ancestor="$(git merge-base HEAD $tag)"

  if [ "$ancestor" != "$head" ]; then
    distance="$(git rev-list HEAD ^${ancestor} --count)"

    if [ "$min_distance" = "0" ]; then
      min_distance="$distance"
      min_tag="$tag"
    fi

    if [ "$distance" -lt "$min_distance" ]; then
      min_distance="$distance"
      min_tag="$tag"
    fi
  fi

done

if [ "$min_distance" != "0" ]; then
  echo $min_tag
  exit 0
else
  exit 1
fi
