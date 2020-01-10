if(@ARGV != 12)
{
        die "not enough arguments\n";
}
$file_tped = $ARGV[0];
$file_expr = $ARGV[1];
$file_infile = $ARGV[2];
$file_outfile = $ARGV[3];
$in_GNT_COUNT = $ARGV[4];
$in_GNT_NA = $ARGV[5];
$in_GNT_MAF = $ARGV[6];
$in_GENE_MEAN_MAX=$ARGV[7];
$in_GENE_MEAN_MIN=$ARGV[8];
$in_GENE_SD_MIN=$ARGV[9];
$in_GENE_SD_MAX=$ARGV[10];
$in_GENE_NA=$ARGV[11];
print "Input Parameters\n";
print "file_tped: $file_tped\n";
print "file_expr: $file_expr\n";
print "file_infile: $file_infile\n";
print "file_outfile: $file_outfile\n";
print "in_GNT_COUNT: $in_GNT_COUNT\n";
print "in_GNT_NA: $in_GNT_NA\n";
print "in_GNT_MAF: $in_GNT_MAF\n";
print "in_GENE_MEAN_MAX: $in_GENE_MEAN_MAX\n";
print "in_GENE_MEAN_MIN: $in_GENE_MEAN_MIN\n";
print "in_GENE_SD_MIN: $in_GENE_SD_MIN\n";
print "in_GENE_SD_MAX: $in_GENE_SD_MAX\n";
print "in_GENE_NA: $in_GENE_NA\n";


chomp($file_tped);
open(GNT,"$file_tped") or die " no file found $file_tped\n";
open(GEXP,"$file_expr") or die " no file found $file_expr\n";
open(INFILE,"gunzip -c $file_infile |") or die " no file found $file_infile\n";
open(OUTFILE,"|gzip > $file_outfile") or die " no file found $file_outfile\n";
$line = <GNT>;
$line = <GEXP>;
$line= <INFILE>;
print OUTFILE $line;
while(<GNT>)
{
	chomp($_);
	@a=split(/\t/,$_);
	#print "$a[7] $in_GNT_NA $a[5] $in_GNT_MAF\n";
	if($a[7] < $in_GNT_NA && $a[5] > $in_GNT_MAF)
	{
		#print "$a[7] $in_GNT_NA $a[5] $in_GNT_MAF\n";
		#die;
		@a1=split('_',$a[6]);
		#@a2=split('_',$in_GNT_COUNT);
		#print "@a1\t@a2\n";
		#if($a1[0] > $a2[0] && $a1[1] > $a2[1] && $a1[2] > $a2[2] )
		$a2[0]=$in_GNT_COUNT;
		$a2[1]=$in_GNT_COUNT;
		$a2[2]=$in_GNT_COUNT;
		#die "@a1\t@a2\n";
		if(($a1[0] > $a2[0] && $a1[1] > $a2[1]) || ($a1[1] > $a2[1] && $a1[2] > $a2[2]) || ($a1[2] > $a2[2] && $a1[0] > $a2[0]))
		{
			$gnt{$a[0]}=1;
			#print "$a[0]\n";	
			#die $a[0]."\n";
		}
	}
}
undef(@a);
while(<GEXP>)
{
	chomp($_);
	@a=split(/\t/,$_);
	#die "$_\n";
	#print "$a[2]> $in_GENE_SD_MIN\n";
	#print "$a[3]<$in_GENE_NA && $a[1] > $in_GENE_MEAN_MIN && $a[1] < $in_GENE_MEAN_MAX && $a[2]> $in_GENE_SD_MIN && $a[2] < $in_GENE_SD_MAX\n";
	if($a[3]<$in_GENE_NA && $a[1] > $in_GENE_MEAN_MIN && $a[1] < $in_GENE_MEAN_MAX && $a[2]> $in_GENE_SD_MIN && $a[2] < $in_GENE_SD_MAX)
	{
		$gene{$a[0]}=1;
		#print $a[0]."\n";
	}
}
undef(@a);
while(<INFILE>)
{
	chomp($_);
	@a=split("\t",$_);
	if(exists($gnt{$a[0]}) && exists($gene{$a[1]}))
	{
		print OUTFILE $_."\n";
	}
}
