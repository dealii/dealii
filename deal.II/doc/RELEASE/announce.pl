#!/usr/bin/perl

######################################################################
# How to use this script:
######################################################################
# 1. Produce a file announce-$version, where
#   $version = $Major.$minor[.$patchlevel]
#
# 2. Run this script like
#
#         perl announce.pl 3.3
#
#   This will send email to Guido, Ralf, Wolfgang.
#
# 3. Check if this email was correctly delivered.
#
# 4. Run the script again with additional 'ok':
#
#         perl announce.pl 3.3 ok
#
# 5. Remove addresses with delivery problems from the list!!!
#
######################################################################

use strict;

my $os = `uname -s`;

# die "Use on Linux machines only !" unless ($os =~ m'Linux');
die "Call as deal user" unless (`who am i` =~ m/deal/);

my @recipients;

my $test = 1;

die "Usage: perl announce.pl <version> [ok]" if ($#ARGV<0);

my $version = $ARGV[0];

my $file = "announce-$version";
die "Announcement file $file does not exist" unless (-r $file);


$test = 0 if ($ARGV[1] eq 'ok');

print "=====Announcing version $version=====\n";

######################################################################
# Undeliverable addresses are commented out.
######################################################################

if ($test)
{
    @recipients = ('wolf', 'guido.kanschat@gmx.net', 'hartmann');
} else {
    @recipients = (
		  'deal@iwr.uni-heidelberg.de',
#		  'kc@isc.tamu.edu',
		  'tveldhui@extreme.indiana.edu',
		  'sullivan@mathcom.com',
		  'Ian_MacPhedran@engr.usask.ca',
		  'roger@maths.grace.cri.nz',
		  'oon-digest@oonumerics.org',
		  'scicomp@uni-erlangen.de',
		  'na.digest@na-net.ornl.gov',
		  'num.info@hermes.iwr.uni-heidelberg.de',
		  'kollegiaten@iwr.uni-heidelberg.de');
}
my $r;

foreach $r (@recipients)
{
    print "$r\n";

    if ($test)
    {
	system ("mailx -s 'deal.II Version $version released' -r deal\@iwr.uni-heidelberg.de $r < $file");
    } else {
	system ("mailx -s 'deal.II Version $version released' -r deal\@iwr.uni-heidelberg.de $r < $file");
    }
}
