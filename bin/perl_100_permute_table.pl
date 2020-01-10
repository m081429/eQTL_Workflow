$value_sort =3;#pvaluea
$end_pval=$ARGV[2];
chomp($end_pval);
$y =0;
for($i=$end_pval;$i> 1e-74;$i = $i/2)
{
	$pval[$y] = $i;
	$y++;
}
@pval = reverse(@pval);
$file = $ARGV[0];
chomp($file);
=head
if(!(-e $file))
{
	die "no file $file exists\n";
}
elsif(-e "$file.gz")
{
	$file="$file.gz";
}
else
{
	system("gzip $file");
	$file="$file.gz";
}
=cut
open(BUFF,"gunzip -c $file.gz|") or die "no input $file file found\n";
$outfile = $ARGV[1];
chomp($outfile);
open(WRBUFF,">$outfile") or die " no $outfile exists\n";
$line = <BUFF>;
@cis=();
@trans=();
for($y=0;$y<@pval;$y++)
{
	$cis[$y] = 0;
	$trans[$y] = 0;
}
while($line=<BUFF>)
{
	@array = split("\t",$line);
	$val = $array[$value_sort];
	$type = $array[5];
	for($y=0;$y<@pval;$y++)
	{
		if($val < $pval[$y])
		{
			if($type == 1)
			{
				$cis[$y]++;	
			}
			else
			{
				$trans[$y]++;	
			}
			last;
		}
	}
}
$cis = 0;
$trans = 0;
for($y=0;$y<@pval;$y++)
{
	$cis = $cis+$cis[$y];
	$trans = $trans+$trans[$y];
	$cis_trans1[$y] = $cis;
	$cis_trans2[$y] = $trans;
}
$line_cis = join("\t",@cis_trans1);
$line_trans = join("\t",@cis_trans2);
print WRBUFF $line_cis."\t".$file."\n";
print WRBUFF $line_trans."\t".$file."\n";
close(CIS);
close(TRANS);
close(BUFF);
close(WRBUFF);
