$table = $ARGV[0];
chomp($table);
open(BUFF,$table) or die "no file found $table\n";
$line=<BUFF>;
chomp($line);
@array=split("\t",$line);
my $cis_ori,$trans_ori,$cis_fdr,$trans_fdr;
for($i=0;$i<@array;$i++)
{
	if($array[$i] eq "CIS_ORI")
	{	
		$cis_ori=$i;
	}
	if($array[$i] eq "TRANS_ORI")
        {
                $trans_ori=$i;
        }
	if($array[$i] eq "CIS_FDR")
        {
                $cis_fdr=$i;
        }
	if($array[$i] eq "TRANS_FDR")
        {
                $trans_fdr=$i;
        }			
}
#print "$cis_ori,$trans_ori,$cis_fdr,$trans_fdr\n";
$line=<BUFF>;
my @cis,@trans;
while($line=<BUFF>)
{
	chomp($line);
	@array=split("\t",$line);
	for($i=@cis;$i<$array[$cis_ori];$i++)
	{
		push(@cis,$array[$cis_fdr]/100);
	}
	for($i=@trans;$i<$array[$trans_ori];$i++)
        {
                push(@trans,$array[$trans_fdr]/100);
        }
}
undef(@array);
#print @cis." ".@trans." "."\n";
$infile = $ARGV[1];
chomp($infile);
open(BUFF,"gunzip -c $infile |") or die "no file found $infile\n";
$line=<BUFF>;
chomp($line);
print "$line\tperm_fdr\n";
$cis=0;
$trans=0;
while($line=<BUFF>)
{
        chomp($line);
	@array=split("\t",$line);
	if($array[5] ==1)
	{
		print $line."\t$cis[$cis]\n";
		$cis++;
	}
	if($array[5] ==0)
        {
                print $line."\t$trans[$trans]\n";
		$trans++;
        }		
}
