# $Id$

# Convert the contents of the OK file to a key readable by the
# database system for regression tests

# used by tests/Makefile.rules


my $contents = <>;

$contents =~ s/\s+//g;

if ($contents eq 'diffok')
{
    print '  + ';
}
elsif ($contents eq 'Compiling')
{
    print ' 0  ';
}
elsif ($contents eq 'Linking')
{
    print ' 1  ';
}
elsif ($contents eq 'Running')
{
    print ' 2  ';
}
elsif ($contents eq 'Checking')
{
    print ' 3  ';
}
else
{
    print ' ?? ';
}
