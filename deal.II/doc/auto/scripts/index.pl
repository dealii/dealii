@entries = ();

while(<>)
{
    next if (m/DeclException/);
    chop;
    if (m/(.+)=(.+)/)
    {
	@n = split "::", $1;

	$entry = pop @n;
	while ($#n>=0)
	{
	    $entry .= " in " . pop @n;
	}
	$entry .= '=' . $2 . '=' . $library;
	push @entries, $entry;
    }
    else
    {
	$library = $_;
    }
}

print << 'EOT'
<HTML><HEAD><TITLE>DEAL Name Index</TITLE></HEAD>
<BODY><H1>Name Index of DEAL</H1><hr>
Remark: Lowercase is sorted behind uppercase and destructors are isolated at
the end.
<hr>
<UL>
EOT
    ;
foreach $entry (sort @entries)
{
    @l = split "=", $entry;
    
    @f = split "#", @l[1];

    print "<LI> <A HREF = \"$l[2]/@f[0]\">@l[0]</A>\n";
}
print << 'EOT'
</UL>
</BODY>
</HTML>
EOT
    ;
