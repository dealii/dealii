# $Id$
# Copyright W. Bangerth, University of Heidelberg, 1998


#Make a dependency file tree
#usage: make_dep library_base_name -Iinc_path1 -Iinc_path2 ... files

#This program makes for each of the given files a makefile dependency
#list, also considering nested includes. It only considers included
#files which are located in the given include pathes (you can give any
#number of pathes). The output looks like this:
#
#lib_base_name.a(.o-file): file included_files
#
#secondly, a similar entry is created for the .g.a library with the
#.go-file.

#Author: Wolfgang Bangerth, March 1998




# list of include pathes (from command line)
@include_path = ();
# list of files to be checked (from command line)
@input_files  = ();
# associate list of files together with the files they include
%include_files = ();


$library_base_name = shift;


# fill include pathes
while ($ARGV[0] =~ /^-I/) {
    $_ = shift;
    /^-I(.*)/;
    $_ = $1;
    if (m![^/]$!) {
	@include_path = (@include_path, $_ . "/");
    } else {
	@include_path = (@include_path, $_);
    }
}

#fill list of files to be processed
while ($ARGV[0]) {
    @input_files = (@input_files, shift);
}


foreach $file (@input_files) {
    make_include_tree ($file);
};


foreach $file (keys %include_files) {
    # complete list of included files by nested include files
    foreach $include (split(' ',$include_files{$file})) {
	complete ($file, $include);
    }
}

# print dependency list
foreach $file (@input_files) {
    $file =~ /(.*)\.(cc)/;

    print "$library_base_name.a($1.o):";
    print "\\\n    $file";
    foreach $f (split (' ', $include_files{$file})) {
	print "\\\n    $f";
    }
    print "\n";

    print "$library_base_name.g.a($1.go):";
    print "\\\n    $file";
    foreach $f (split (' ', $include_files{$file})) {
	print "\\\n    $f";
    }
    print "\n";
}








# complete the list of included files by the files
# included by a file included by the original one
sub complete {
    local ($file, $include) = ($_[0], $_[1]);
    foreach $second_include (split(' ',$include_files{$include})) {
	if (! ($include_files{$file} =~ $second_include)) {
	    #second_include not yet in list of included files
	    $include_files{$file} =
		join(' ', $second_include, $include_files{$file});
	    complete ($file, $second_include);
	}
    }
}




# make the include tree for a file
sub make_include_tree {
    local ($filename) = $_[0];

    open (FILE, $filename);
    while (<FILE>) {
	# look out for include statements
	if (/^#\s*include\s*(["<])([^">]*)[">]/) {
	    local($include) = $2;
	    if ($1 =~ /</) {
		# include by <...>. Try to find real path
		for $include_dir (@include_path) {
		    $if = $include_dir . $include;
		    if (-r $if) {
			# file found; delete ./ at beginning
			if ($if =~ m!^\./(.*)!) {
			    $if = $1;
			}
			$include_files{$filename} =
			    join (' ', $if, $include_files{$filename});
		    }
		}
	    } else {
		# included by "..."
		# only append file if it exists
		if (-r $include) {
		    $include_files{$filename} =
			join (' ', $include, $include_files{$filename});
		}
	    }
	}
    }

    # for each file included here: make up include tree itself
    for $if (split(/ /,$include_files{$filename})) {
	# if include file list not yet made up
	if (! defined ($include_files{$if})) {
	    make_include_tree($if);
	}
    }
}
