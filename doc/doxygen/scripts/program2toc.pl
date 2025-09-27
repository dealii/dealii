## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2006 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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

print "<br>\n";
