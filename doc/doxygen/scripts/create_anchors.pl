## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2007 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


# If we find a heading in a .dox file, create an HTML anchor for it.
use Getopt::Long;
GetOptions("prefix:s" => \$prefix);
while (<>) {
    if ( /<h.>(.*)<\/h.>\s*/ ) {
	$reftext = $1;

	# for the anchor, use the name of the section but discard
	# everything except for letters, numbers, and underscores
	$reftext =~ s/[^a-zA-Z0-9_]//g;

	print "<a name=\"$prefix-$reftext\"></a>$_\n";
    } else {
        print;
    }
}
