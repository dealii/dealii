touch ../results/$2/$1/.cvsignore
foreach f (*.output)
    if (-e $f:r.OK) then
	set old=`ls -l ../compare | perl -p -e 's/.*results\///;'`
	pushd ../results/$2/$1
	ln -s ../../$old/$1/$f $f
	echo $f >> .cvsignore
	popd
    else
	pushd ../results/$2/$1    
	cp =1/$f .
	cvs add $f
	popd
    endif
end
