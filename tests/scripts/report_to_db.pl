#!/bin/perl

use strict;

my $path = '.';
my $revision;
my $id;
my $date;

while(<>)
{
    if (m/^Date:\s*(\d+)\s+(\d+)/)
    {
	my $name = sprintf("%s/%04d-%03d", $path, $1, $2);
#	print $name;
    }
    $revision = $1 if (m/^Revision:\s*(\d+)/);
    $id = $1 . $2 if (m/^Id:\s*(\S+)\s+(\S+)?\n/);
    if (m/^\d\d\d\d-\d\d-\d\d \d\d:\d\d/)
    {
	chop;
	s/\s+/ /g;
	s/ \+ / 4 /;
	print "$revision $_ $id\n"
    }
}
