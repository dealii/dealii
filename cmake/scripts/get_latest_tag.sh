#!/bin/sh
#
# Find the latest TAG that has a common ancestry (of positive distance) with
# current HEAD
#

if ! git --version | grep -q "git version 2"; then
  # This script requires version 2 or newer
  exit 1
fi

head="$(git rev-parse HEAD)"

tags="$(git tag --sort=-creatordate)"

for tag in $tags; do
  tag_commit="$(git rev-parse $tag^{commit})"

  if [ "$tag_commit" = "$head" ]; then
    echo "$tag"
    exit 0;
  fi

  ancestor="$(git merge-base HEAD $tag)"

  if [ "$ancestor" != "$head" ]; then
    echo "$tag"
    exit 0;
  fi
done

exit 1
