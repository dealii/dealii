#!/bin/sh

######################################################################
# Find equal reference files
#
# Called from one of the top level test directories, finds all doubled
# reference files in */cmp and outputs them for possible removal.
#
######################################################################

# go through all files
for dir in *; do
  # if file is a directory,
  if ! test -d $dir; then continue; fi
  # change to the subdirectory with reference files
  cd $dir/cmp

  # now diff all possible combinations and
  # write a line if files are equal.
  for file in *; do
    for cmp in *; do
      if test "$file" = "$cmp"; then break; fi
      if diff -q $file $cmp &> /dev/null ; then
        echo "$dir $file $cmp equal"
      fi
    done
   done
   cd ../..
 done
