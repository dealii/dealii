

# Modify these to enter the current data automatically
my $year = 2013;
my $version = "8.1.pre";

if (m'</head>')
{
    print '<link href="$relpath$deal.css" rel="stylesheet" type="text/css"></link>', "\n";
    print '<link rel="SHORTCUT ICON" href="http://www.dealii.org/deal.ico"></link>', "\n";
    print '<meta name="author" content="The deal.II Authors <authors@dealii.org>"></meta>', "\n";
    print '<meta name="copyright" content="Copyright (C) 1998 - ', $year, ' by the deal.II authors"></meta>', "\n";
    print '<meta name="deal.II-version" content="', $version, '"></meta>', "\n";
}
