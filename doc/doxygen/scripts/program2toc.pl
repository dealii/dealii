## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

print "    <ul>\n";
use Getopt::Long;
GetOptions("prefix:s" => \$prefix);

$level = 3;
while (<>) {
    if ( /\@sect(.)\{(.*)\}/ ) {
	$newlevel = $1;
	$text = $2;
	
	# only allow header levels 3 through 6, since higher ones are
	# reserved for top-level document headers
	if (! ($newlevel =~ /[3456]/)) {
	    print STDERR "Only header levels 3 through 6 are allowed.\n";
	    print STDERR "You had $newlevel.\n";
	    die;
	}

	if ($newlevel > $level) {
	    for ($i=$level; $i<$newlevel; ++$i) {
	        print "      <ul>\n";
            }
	} elsif ($newlevel < $level) {
	    for ($i=$newlevel; $i<$level; ++$i) {
	        print "      </ul>\n";
            }
	}

	$reftext = $text;

	# for the anchor, use the name of the section but discard
	# everything except for letters, numbers, and underscores
	$reftext =~ s/[^a-zA-Z0-9_]//g;

	# replace quotation marks by the appropriate HTML quotation marks
	$text =~ s!``!&#8220;!g;
	$text =~ s!''!&#8221;!g;

        # replace double dashes in comments by &mdash;
	$text =~ s!--!&mdash;!g;

	print "        <li><a href=\"#$prefix-$reftext\">$text</a>\n";

	$level = $newlevel;
    } 
}

for (; $level>=3; --$level) {
    print "      </ul>\n";
}
