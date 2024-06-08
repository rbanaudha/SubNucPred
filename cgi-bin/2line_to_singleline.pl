#!/usr/bin/perl

open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    #print "$line\n";
    if($line=~ m/^>/)
    {
	print "\n$line:";
    }
    if($line !~ m/^>/)
    {
	print "$line";
    }
}
print "\n";
