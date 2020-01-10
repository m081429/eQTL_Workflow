#$input = "193sgenome";
$dir = $ARGV[0];
chomp($dir);
open(BUFF,"$dir/input_summary/Tped_Summary") or die "no file $dir/Tped_Summary\n";
$line=<BUFF>;
while(<BUFF>)
{
	@array=split("\t",$_);
	$hash{$array[0]}=$array[4];
}
close(BUFF);
open(BUFF,"$dir/processing/processed_inputdata.tped") or die "no $dir/processed_inputdata.tped file found\n";
#system("head -1 GE.txt > SNP.txt");
#system("sed 's/gene/rs/g' SNP.txt > SNP1.txt");
#system("mv SNP1.txt SNP.txt");
open(BUFF1,"$dir/processing/processed_inputdata.tfam") or die "no $dir/processed_inputdata.tfam file found\n";


open(WRBUFF,">$dir/processing/SNP.txt");
open(WRBUFF1,">$dir/processing/snploc.txt");
print WRBUFF1 "SNP\tchrm_snp\tpos\n";
print WRBUFF "rsid";
while(<BUFF1>)
{
	chomp($_);
	@a = split(" ",$_);
	print WRBUFF "\t$a[1]";
}
print WRBUFF "\n";
while(<BUFF>)
{
	chomp($_);
	@array = split(" ",$_);
	$chr = shift(@array);
	if($chr ne "0")
	{
	$snp = shift(@array);
	shift(@array);
	$pos = shift(@array);
	print WRBUFF1 "$snp\t$chr\t$pos\n";
	@uniq = uniq(@array);
	@uniq = sort(@uniq);
	if($uniq[0] eq "0")
	{
		shift(@uniq);
	}
	$gnt1 = $uniq[0];
	$gnt2 = $uniq[1];
	if(exists($hash{$snp}) && $gnt1 ne $hash{$snp})
	{
		$tmp = $gnt1;
		$gnt1 = $gnt2;
		$gnt2 = $tmp;	
	}
	print WRBUFF $snp;
	for($i=0;$i<@array;$i++)
	{
		$snp1 = $array[$i];
		$i++;
		$snp2 = $array[$i];
		if($snp1 eq $gnt1 && $snp2 eq $gnt1)
		{
			print WRBUFF "\t2";
		}
		elsif($snp1 eq $gnt2 && $snp2 eq $gnt2)
                {
                        print WRBUFF "\t0";
                }
		elsif($snp1 eq "0" && $snp2 eq "0")
                {
                        print WRBUFF "\tNA";
                }
		else
		{
			print WRBUFF "\t1";
		}	
	}
	print WRBUFF "\n";
	#print "$snp\t@uniq\n"; 
#die;	
	}
}
sub uniq {
    return keys %{{ map { $_ => 1 } @_ }};
}
close(BUFF);
close(WRBUFF);
close(WRBUFF1);
close(BUFF1);
