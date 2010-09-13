#!/bin/bash

cd $1/cmp
echo Files: *
for name in * ; do
    echo Comparing file $name
    if diff -w ../output $name > ../DIFF_$name ; then
	exit;
    else
	less ../DIFF_$name
	rm ../DIFF_$name
	echo "Replace compare file? (y/n)"
	read
	if test "x$REPLY" == "xy" ; then
	    cp ../output $name
	fi
    fi
done
