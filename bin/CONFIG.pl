#!/usr/bin/perl
sub getDetails
{
	my ($configFile) = @_;
  	open(BUFF,"$configFile");
	while($line = <BUFF>)
	{
		chomp($line);
		my @arr = split('=',$line);
		chomp($arr[0]);
		#$arr[0]=~ s/_/./g;
		chomp($arr[1]);
		$config{$arr[0]} = "$arr[1]";
		#print "$arr[0] $arr[1]\n";			
	}
	return \%config;
}
return \%config;
#getDetails("kk.txt","geneEXP");
