# $Id$
# Check whether references in HTML files are valid or 
# point to non-existing files/links/etc
#
# Author: Wolfgang Bangerth, 2000


# set this to 1 if you want verbose output
$debug = 0;

foreach $filename (@ARGV) {
    open IN, $filename
        or die "---Can't open file `$filename'\n";

    print "File: $filename\n" if $debug;
    while (<IN>) {
        # first find all hrefs
        while ( /href=\"?(.*?)[\s\"]/gi ) {
	    # then decide whether they are relevant for 
            # our purpose
	    $link = $1;
    
	    if ( $link =~ /^mailto|http/i ) {
	        # link is external. don't check it
	        print "external link: $link\n" if $debug;
	        next;
	    }
	    elsif ( $link =~ /^#(.*)/ )
	    {
		# this is a reference within this file. try to 
	        # find its anchor
	        $internal_ref = $1;
	        print "internal reference: $link\n" if $debug;
	    
	        open IN2, $filename;
	        $found = 0;
	        while ( <IN2> ) {
		    if ( /<a name=\"?(.*?)[\s\"]/i ) {
		        if ( $1 eq $internal_ref) {
			    print "                    found.\n" if $debug;
			    $found = 1;
			    last;
			}
		    }
		}
		
		die "---Internal reference `$internal_ref' not found in file $filename\n"
		    unless $found;
		next;
	    }
	    elsif ( $link =~ /^(.*?)#(.*)/ )
	    {
		# this is a reference within another file. try to 
		# find its anchor
		$external_file = $1;
		$external_ref = $2;
		print "external reference: $link\n" if $debug;
		
		open IN2, $external_file;
		$found = 0;
		while ( <IN2> ) {
		    if ( /<a name=\"?(.*?)[\s\"]/i ) {
			if ( $1 eq $external_ref) {
			    print "                    found.\n" if $debug;
			    $found = 1;
			    last;
			}
		    }
		}
		
		die "---External reference `$internal_ref' not found in file $filename\n"
		    unless $found;
		next;
	    }
	    else {
		# this must now be a regular file which is
		# referenced. the file must be local
		die "---Local file `$link' not found in file `$filename'\n"
		    unless ((-r $link) && (-f $link));
	    }
		}
		}
}

