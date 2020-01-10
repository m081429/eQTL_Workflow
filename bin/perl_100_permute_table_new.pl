$value_sort =3;#pvaluea
$end_pval=$ARGV[3];
$cisfile = $ARGV[1];
$outfile = $ARGV[2];
$transfile = $ARGV[0];
chomp($end_pval);
chomp($cisfile);
chomp($outfile);
chomp($transfile);
$y =0;
for($i=$end_pval;$i> 1e-74;$i = $i/2)
{
	$pval[$y] = $i;
	$y++;
}
@pval = reverse(@pval);

open(CBUFF,"gunzip -c $cisfile.gz|") or die "no input $cisfile file found\n";
open(TBUFF,"gunzip -c $transfile.gz|") or die "no input $transfile file found\n";
open(WRBUFF,">$outfile") or die " no $outfile exists\n";
$line = <CBUFF>;
$line = <TBUFF>;
@cis=();
@trans=();
for($y=0;$y<@pval;$y++)
{
	$cis[$y] = 0;
	$trans[$y] = 0;
}
while($line=<CBUFF>)
{
	@array = split("\t",$line);
	$val = $array[$value_sort];
	for($y=0;$y<@pval;$y++)
	{
		if($val < $pval[$y])
		{
				$cis[$y]++;	
			last;
		}
	}
}
while($line=<TBUFF>)
{
        @array = split("\t",$line);
        $val = $array[$value_sort];
        for($y=0;$y<@pval;$y++)
        {
                if($val < $pval[$y])
                {
                        $trans[$y]++;
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
print WRBUFF $line_cis."\t".$cisfile."\n";
print WRBUFF $line_trans."\t".$transfile."\n";
close(CIS);
close(TRANS);
close(BUFF);
close(WRBUFF);
