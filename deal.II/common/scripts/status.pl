$/ = 'File:';

my $compress = 0;
my $short = 0;
my $cvs_args = '';
my $process_cvs_args = 0;
my $debug = 0;

foreach (@ARGV)
{
    if ($process_cvs_args)
    {
	$cvs_args .= " $_";
    }
    $compress = 1 if (m/-c/);
    $short = 1 if (m/-s/);
    $debug = 1 if (m/-d/);
    $process_cvs_args = 1 if (m/^--$/);
}

my $format = '(\S+)\s*Status: (.+)\s*Working revision:\s*(\S+)'
    . '\s*Repository revision:\s*(\S+)\s*(\S+)'
    . '\s*Sticky Tag:\s*(.+)'
    . '\s*Sticky Date:\s*(.+)'
    . '\s*Sticky Options:\s*(.+)'
    ;

print STDERR "cvs status $cvs_args |" if ($debug);
open CVS, "cvs status $cvs_args |";

while(<CVS>)
{
    next if (m/^\?/);
    next if (m/^===================================================================/);
    if (m/$format/)
    {
	my $file = $1;
	my $status = $2;
	my $work = $3;
	my $rep = $4;
	my $rep_file = $5;
	my $tag = $6;
	my $date = $7;
	my $opt = $8;

	$rep_file =~ s!/home/people/cvs/deal/!!;
	$rep_file =~ s!,v!!;

	if ($status eq 'Up-to-date')
	{
	    $status = 'OK ';
	} elsif ($status eq 'Needs Patch')
	{
	    $status = "P $work<-$rep ";
	} elsif ($status eq 'Locally Modified')
	{
	    $status = "M $work->$rep ";
	} else {
	    $status = "$status $work $rep ";
	}

	$tag = '' if ($tag eq '(none)');
	$date = '' if ($date eq '(none)');
	$opt = '' if ($opt eq '(none)');
	my $stick = $tag . $date . $opt;
	$stick = "Stick $stick" unless ($stick eq '');
	$status .= $stick;
	printf "%-30s %-40s %s\n", $file, $status,
	(($short) ? '' : $rep_file)
	    unless ($compress && ($status eq 'OK '));
    }
    else
    {
	print STDERR "Ignored:\n$_\n",
	"===================================================================\n"
	    if ($debug);
    }
}
