#!/usr/bin/perl

open(FH,"$ARGV[0]") or die "$!";
while($line=<FH>)
{
    chomp($line);
    #print "$line\n";
    @array=split(/:/,$line);
    print "$array[0]\n$array[1]\n";
}
