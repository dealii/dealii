## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2013 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Author: Wolfgang Bangerth, Guido Kanschat 2000, 2004
#
# Check whether references in HTML files are valid or
# point to non-existing files/links/etc
#


# set this to 1 if you want verbose output
$debug = 0;

$startdir = `pwd`;
chop $startdir;

foreach $filename (@ARGV)
{
    chdir $startdir || die "Could not change dir to $startdir\n";
    open IN, $filename
        or die "---Can't open file `$filename'\n";

    print "File: $filename\n" if $debug;
    if ($filename =~ m!(.+)/([^/]+)!)
    {
	chdir $1;
	$filename = $2;
    }

    while (<IN>) {
	# save the entire line for simpler grepping when an error
	# occurs
	$this_line = $_;

	# if line ends with an = character, the concatenate it with the next
	# one
	while ( /=\s*$/ ) {
	    $newline = <IN>;
	    $newline =~ s/^\s*//g;
	    $_ = $thisline . $newline;
	    $this_line = $_;
        }

        # first find all hrefs
        while ( /<\s*a\s+href=\"?(.*?)[\s\"]/gi ) {
	    # then decide whether they are relevant for
            # our purpose
	    $link = $1;

	    if ( $link =~ /^mailto|http(s)?:\/\//i ) {
	        # link is external. don't check it
	        print "external link: $link\n" if $debug;
	        next;
	    }
	    elsif ( $link =~ m/^#(.*)/ )
		{
			# this is a reference within this file. try to
	        # find its anchor
	        $internal_ref = $1;
	        print "internal reference: $link\n" if $debug;

	        open IN2, $filename;
	        $found = 0;
	        while ( <IN2> ) {
		    while ( /<a[^>]* (name=|class=\"anchor\" id=)\"?(.*?)[\s\"]/gi ) {
		        if ( $2 eq $internal_ref)
			{
			    print "                    found.\n" if $debug;
			    $found = 1;
			    last;
			}
		    }
		}

		die "---Internal reference `$internal_ref' not found in file $filename\n This line is: $this_line.\n"
		    unless $found
                            # work around a bug in doxygen 1.6.3:
			    #   https://bugzilla.gnome.org/show_bug.cgi?id=620372
                            || ($internal_ref =~ /^index_[:_~]/) ;
		next;
	    }
	    elsif ( $link =~ /^(.*?)#(.*)/ )
	    {
		# this is a reference within another file. try to
		# find its anchor
		$external_file = $1;
		$external_ref = $2;

		# if the file name was prepended with http: (but is a local file,
		# so no double-slash), then split off http:
		$external_file =~ s/^http(s)?://g;

		print "external reference: $link\n" if $debug;

		open IN2, $external_file;
		$found = 0;
		while ( <IN2> ) {
		    while ( /<a[^>]* (name=|class=\"anchor\" id=)\"?(.*?)[\s\"]/gi ) {
			if ( $2 eq $external_ref)
			{
			    print "                    found.\n" if $debug;
			    $found = 1;
			    last;
			}
		    }
		}

		die "---External reference `$external_file#$external_ref' not found in file $filename\n This line is: $this_line.\n"
		    unless $found
                            # work around a bug in doxygen 1.6.3:
			    #   https://bugzilla.gnome.org/show_bug.cgi?id=620372
                            || ($external_ref =~ /^index_[:_~]/) ;
		next;
	    }
	    else {
		# this must now be a regular file which is
		# referenced. the file must be local

		# if the file name was prepended with http: (but is a local file,
		# so no double-slash), then split off http:
		$link =~ s/^http(s)?://g;

		die "---Local file `$link' not found in file `$filename'\n This line is: $this_line.\n"
		    unless ((-r $link) && (-f $link));
	    }
	}

	# check whether references to images are valid
	while ( /img\s+src=\"?(.*?)[\s\"]/gi ) {
	    # check whether the file for the image is present
	    $link = $1;

	    # ignore online links
	    if ($link =~ /^http/)
	    {
		next;
	    }

	    die "---Local image `$link' not found in file `$filename'\n This line is: $this_line.\n"
		unless ((-r $link) && (-f $link));
	}
   }
}

