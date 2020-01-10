#!/usr/bin/perl
#
# This script is used to run the eQTL workflow 


#get current directory
use Cwd 'abs_path';
$line = abs_path($0);
chomp $line;
@DR_array = split('/',$line);
pop(@DR_array);
$dir = join("/",@DR_array);

use Getopt::Long;
&Getopt::Long::GetOptions(
'run_config=s'      => \$config,
'tool_config=s'      => \$toolconfig
);
if($config eq "" || $toolconfig  eq "")
{
	print "config file missing.Check for this program manual and create a config file \n";
	
	die "Retry with : perl  perl_workflow_eqtl.pl  -run_config <PATH TO THE RUN CONFIG FILE> -tool_config <PATH TO THE TOOL CONFIG FILE>\n";
}
require "$dir/bin/CONFIG.pl";
#reading config info parameters
getDetails($config);
my $expr= $config{"expr.data"};
my $tped = $config{"tped"};
my $tfam = $config{"tfam"};
my $outdir = $config{"output.dir"};
my $dirtemp = $config{"temp.folder"};
my $summ_ind = $config{"only.summary"};
my $expr_na = $config{"expr.na.symbol"};
my $expr_loc = $config{"expr.loc"};
my $impute_queue = $config{"cluster.queue"};
my $impute_mem = $config{"cluster.max.mem"};
my $email = $config{"email"};
my $permutation=$config{"num.permutation"};
my $cis_pval_in=$config{"cis.pval"};
my $trans_pval_in=$config{"trans.pval"};
my $rounded = $config{"run.dir"};
my $covariates = $config{"covariate.data"};
my $sge=$config{"sge"};
my $filter= $config{"filtering"};
my $filter_geno_count= $config{"fil.min.num.gnt"};
my $filter_geno_na= $config{"fil.max.na.gnt"};
my $filter_geno_maf= $config{"fil.maf.cutoff.gnt"};
my $filter_gene_mean_max= $config{"fil.max.mean.geneexp"};
my $filter_gene_mean_min= $config{"fil.min.mean.geneexp"};
my $filter_gene_sd_max= $config{"fil.max.std.geneexp"};
my $filter_gene_sd_min= $config{"fil.min.std.geneexp"};
my $filter_gene_na= $config{"fil.max.na.geneexp"};
my $cisDist= $config{"cisDist"};
my $geneannot_attach=$config{"geneannot.attach.genesumm"};
#reading tool info parameters
getDetails($toolconfig);
my $plink= $config{"path_plink"};
my $perl=$config{"path_perl"};
my $sh=$config{"path_sh"};
my $qsub=$config{"path_qsub"};
my $rscript=$config{"path_rscript"};
my $rlib=$config{"path_rlib"};
#tool info
$rscript=~ s/\s|\t|\r|\n//g;
$perl=~ s/\s|\t|\r|\n//g;
$sh=~ s/\s|\t|\r|\n//g;
$qsub=~ s/\s|\t|\r|\n//g;
$rscript=~ s/\s|\t|\r|\n//g;
$rlib =~ s/\s|\t|\r|\n//g;

#runinfo
$expr =~ s/\s|\t|\r|\n//g;
$tped =~ s/\s|\t|\r|\n//g;
$tfam =~ s/\s|\t|\r|\n//g;
$outdir =~ s/\s|\t|\r|\n//g;
$dirtemp =~ s/\s|\t|\r|\n//g;
$summ_ind =~ s/\s|\t|\r|\n//g;
$expr_na =~ s/\s|\t|\r|\n//g;
$expr_loc =~ s/\s|\t|\r|\n//g;
$impute_queue =~ s/\s|\t|\r|\n//g;
$impute_mem =~ s/\s|\t|\r|\n//g;
$email =~ s/\s|\t|\r|\n//g;
$filter_geno_count =~ s/\s|\t|\r|\n//g;
$filter_geno_na =~ s/\s|\t|\r|\n//g;
$filter_geno_maf =~ s/\s|\t|\r|\n//g;
$filter_gene_mean_max =~ s/\s|\t|\r|\n//g;
$filter_gene_mean_min =~ s/\s|\t|\r|\n//g;
$filter_gene_sd_max =~ s/\s|\t|\r|\n//g;
$filter_gene_sd_min =~ s/\s|\t|\r|\n//g;
$filter_gene_na =~ s/\s|\t|\r|\n//g;
$permutation  =~ s/\s|\t|\r|\n//g;
$trans_pval_in =~ s/\s|\t|\r|\n//g;
$cis_pval_in =~ s/\s|\t|\r|\n//g;
$sge =~ s/\s|\t|\r|\n//g;
$rounded =~ s/\s|\t|\r|\n//g;
$covariates =~ s/\s|\t|\r|\n//g;
$filter =~ s/\s|\t|\r|\n//g;
$cisDist =~ s/\s|\t|\r|\n//g;
$geneannot_attach =~ s/\s|\t|\r|\n//g;
print "***********RUNINFO CONFIG ARGUMENTS***********\n";
print "expr_data: $expr\n";
print "tped: $tped\n";
print "tfam: $tfam\n";
print "output.dir: $outdir\n";
print "temp.folder: $dirtemp\n";
print "only.summary: $summ_ind\n";
print "expr.na.symbol: $expr_na\n";
print "expr.loc: $expr_loc\n";
print "cluster.queue: $impute_queue\n";
print "cluster.max.mem: $impute_mem\n";
print "email: $email\n";
print "num.permutation : $permutation\n";
print "cis.pval: $cis_pval_in\n";
print "trans.pval: $trans_pval_in\n";
print "sge: $sge\n";
print "run.dir: $rounded\n";
print "filtering: $filter\n";
print "cisDist : $cisDist\n";
print "geneannot.attach.genesumm : $geneannot_attach\n";
if($filter eq "YES")
{
	print "fil.min.num.gnt: $filter_geno_count\n";
	print "fil.max.na.gnt: $filter_geno_na\n";
	print "fil.maf.cutoff.gnt: $filter_geno_maf\n";
	print "fil.max.mean.geneexp: $filter_gene_mean_max\n";
	print "fil.min.mean.geneexp: $filter_gene_mean_min\n";
	print "fil.max.std.geneexp: $filter_gene_sd_max\n";
	print "fil.min.std.geneexp: $filter_gene_sd_min\n";
	print "fil.max.na.geneexp: $filter_gene_na\n";
}
print "***********TOOLINFO CONFIG ARGUMENTS***********\n";
print "path_rscript: $rscript\n";
print "path_perl: $perl\n";
print "path_sh: $sh\n";
print "path_qsub: $qsub\n";
print "path_plink: $plink\n";
print "path_rlib : $rlib\n";


#print "MODEL : $model\n";
if($covariates ne "")
{
	print "COVARIATES: $covariates\n";
}
if( $geneannot_attach eq "" |  $cisDist eq "" | $filter eq "" | $rounded eq "" | $sge eq "" | $trans_pval_in eq "" | $cis_pval_in eq "" |$permutation eq "" |  $email eq "" | $impute_mem eq "" |$impute_queue eq "" |$expr_loc eq "" | $expr_na eq "" | $dirtemp eq "" | $tped eq "" | $tfam eq "" | $outdir eq "" | $expr eq "" | $summ_ind eq "")
{
	die "some of the input run info parameters are Empty\n";
}

if($filter eq "YES")
{
	if($filter_geno_count eq "" | $filter_geno_na eq "" | $filter_geno_maf eq "" | $filter_gene_mean_max eq "" | $filter_gene_mean_min eq "" | $filter_gene_sd_max eq "" | $filter_gene_na eq "")
	{
		die "Filter : YES.so filter options cannot be empty\n";
	}	
}
if($sge eq "YES" & ($qsub eq "" || !(-e $qsub)))
{
	die "some of the input tool info parameter qsub $qsub is Empty or invalid\n";
}
if(!(-e $sh) | !(-e $perl) | !(-e $plink) | !(-e $rlib) | !(-d $rlib))
{
	die "some of the input tool info parameters are Empty or invalid\n";
}
if(!(-e $tped))
{
	die "input tped ile does not exist\n";
}
if(!(-e $tfam))
{
	die "input tfam file does not exist\n";
}
if(!(-e $expr))
{
	die "input expr file does not exist\n";
}
unless(-d $dirtemp)
{
    system("mkdir -p $dirtemp");
}
#if($model ne "modelLINEAR" && $model ne "modelANOVA")
#{
#	die "Model should be either modelLINEAR or  modelANOVA\n";
#}
if(uc($rounded) eq "NA")
{
	$round = sprintf("%.0f", rand()*time());
	$rounded = "temp".$round;
}
#$rounded = "temp_eqtl";
$dirtemp ="$dirtemp/$rounded";

#creating the temporary directory
system("mkdir -p $dirtemp");
system("mkdir $dirtemp/input_summary");
system("mkdir $dirtemp/sgelog");
system("mkdir $dirtemp/processing");
system("mkdir $dirtemp/eqtl_run");

#die "$dirtemp\n";
$outdir=$outdir."/$rounded";
unless(-d $outdir)
{
    system("mkdir -p $outdir");
}
#summary begins
#system("cp $tfam $dirtemp/processed_inputdata.tfam");
#reading tped file and creating the summary file
open(TFAM,$tfam) or die "no $tfam file exists\n";
open(WRTFAM,">$dirtemp/processing/processed_inputdata.tfam") or die "not able to write the tfam file\n";
$i = 1;
while(<TFAM>)
{	
	chomp($_);
	@tfam=split(" ",$_);
	$sampleorder{$tfam[1]} = $i++;
	print WRTFAM $_."\n";
}

open(TPED,$tped) or die "no $tped file exists\n";
open(WRBUFF,">$dirtemp/input_summary/Tped_Summary") or die " not able to write the file tped summary\n";
open(WRTPED,">$dirtemp/processing/processed_inputdata.tped") or die "not able to write the tped file\n";
print WRBUFF "SNPid\tChr\tPos\tMajorAllele\tMinorAllele\tMinorAlleleFreq\tgenotype_distribution\tnum_na\n";
open(EXPR,"$expr") or die "not able to read the expression file $expr\n";
open(EXPR_LOC,"$expr_loc") or die "not able to read the expression file $expr_loc\n";
open(WREXPR,">$dirtemp/processing/expression.txt") or die "not able to write the file $dirtemp/processing/expression.txt\n";
open(WRGSUMM,">$dirtemp/input_summary/Gene_Summary") or die " not able to write the file gene summary\n";
print WRGSUMM "SegmentID\tmean\tstd\tnum_na\tChr\tStart\tStop\n";

$line=<EXPR>;
chomp($line);
@expr=split(/\t/,$line);
$num_geneexp=@expr-1;
$i--;
if($num_geneexp != $i)
{
	die "number of samples in the plink files are $i and number of samples in the expression are $num_geneexp\n";
}
@order=();
for($i=1;$i<@expr;$i++)
{	
	if($i > 0)
	{
		if(!exists($sampleorder{$expr[$i]}))
		{
			die "sample exist in genotype but not in the expression\n";
		}
		push(@order,$sampleorder{$expr[$i]});	
	}	
}
#print @order."\n";

while(<TPED>)
{
	chomp($_);
	@tped=split(" ",$_);
	$chr = shift(@tped);
	$rsid=shift(@tped);
	$dist = shift(@tped);
	$pos=shift(@tped);
	$_=join(" ",@tped);
	@allele=();
	if(uc($_) =~ m/A/)
	{
		unshift(@allele,"A");
	}
	if(uc($_) =~ m/T/)
	{
		unshift(@allele,"T");
	}
	if(uc($_) =~ m/G/)
	{
		unshift(@allele,"G");
	}
	if(uc($_) =~ m/C/)
	{
		unshift(@allele,"C");
	}
	$num = @tped;
	if(@allele > 1 && @allele < 3)
	{
		$count1 = 0;
		$count2 = 0;
		$count3 = 0;
		$missing = 0;
		for($i=0;$i<$num;$i++)
		{
			$a1 = $tped[$i++];
			$a2 = $tped[$i];	
			if($a1 eq $allele[0] && $a2 eq $allele[0])
			{
				$count1++;
			}
			elsif($a1 eq $allele[1] && $a2 eq $allele[1])
			{
				$count2++;
			}
			elsif(($a1 eq $allele[0] && $a2 eq $allele[1] ) || ($a2 eq $allele[0] && $a1 eq $allele[1] ))
			{
				$count3++;
			}	
			else
			{
				$missing++;
			}
		}
		print WRTPED "$chr $rsid $dist $pos $_\n";
		print WRBUFF "$rsid\t$chr\t$pos";
		if($count1 >= $count2)
		{
			$maf = ($count2*2+$count3*1)/(($count1+$count2+$count3)*2);
			$count = $count1.'_'.$count3.'_'.$count2;
			$missing_p = ($missing*2)/$num;
			print WRBUFF "\t$allele[0]\t$allele[1]\t$maf\t$count\t$missing_p\n";
		}
		else
		{
			$maf = ($count1*2+$count3*1)/(($count1+$count2+$count3)*2);
			$count = $count2.'_'.$count3.'_'.$count1;
			$missing_p = ($missing*2)/$num;
			print WRBUFF "\t$allele[1]\t$allele[0]\t$maf\t$count\t$missing_p\n";
		}	
	}
} 

#Gene Expression summary starts
open(EXPR,"$expr") or die "not able to read the file\n";
$line=<EXPR>;
chomp($line);
@expr=split("\t",$line);
undef(@array);
chomp($expr[0]);
$array[0] = $expr[0];

$line_expr_loc=<EXPR_LOC>;
chomp($line_expr_loc);
for($i=1;$i<@expr;$i++)
{
	$array[$order[$i-1]] = $expr[$i];	
}	
$line=join("\t",@array);
print WREXPR "$line\n";
while(<EXPR>)
{
	chomp($_);
	$line_expr_loc=<EXPR_LOC>;
	chomp($line_expr_loc);
	@array_expr_loc=split("\t",$line_expr_loc);
	@expr=split("\t",$_);
	undef(@array);
	$array[0] = $expr[0];
	chomp($array[0]);
	if($array[0] ne $array_expr_loc[0])
	{
		die "Gene $array[0] not in the same recode location for files expression.txt and location.txt\n"; 
	}	
	undef(@calc);
	@calc=();
	$na=0;
	for($i=1;$i<@expr;$i++)
	{
		if($expr[$i] ne $expr_na)
		{
			push(@calc,$expr[$i]);
		}
		else
		{
			$expr[$i]="NA";
			$na++;	
		}
		$array[$order[$i-1]] = $expr[$i];
	}
	#die "@calc\n";
	$ave = &average(\@calc);
	$std = &stdev(\@calc);
	$_=join("\t",@array);
	$na=$na/(@array-1);
	print WREXPR "$_\n";
	print WRGSUMM "$array[0]\t$ave\t$std\t$na\t$array_expr_loc[1]\t$array_expr_loc[2]\t$array_expr_loc[3]\n";
		
}
close(EXPR);
close(EXPR_LOC);
close(WRGSUMM);
#gene expression attach annotation
if($geneannot_attach ne "NA")
{
	$sys="head -1 $dirtemp/input_summary/Gene_Summary > $dirtemp/input_summary/tmp1";
	#print $sys."\n";
	system($sys);
	$sys="head -1 $geneannot_attach|cut -f2->$dirtemp/input_summary/tmp2";
	#print $sys."\n";
	system($sys);
	$sys="paste $dirtemp/input_summary/tmp1 $dirtemp/input_summary/tmp2 > $dirtemp/input_summary/Gene_Summary_tmp";
	#print $sys."\n";
	system($sys);
	$sys="rm $dirtemp/input_summary/tmp1 $dirtemp/input_summary/tmp2";
	#print $sys."\n";
	system($sys);
	$sys="sort -k 1,1 $geneannot_attach >$dirtemp/input_summary/tmp3";
	#print $sys."\n";
	system($sys);
	$sys="sort -k 1,1 $dirtemp/input_summary/Gene_Summary >$dirtemp/input_summary/tmp4";
	#print $sys."\n";
	system($sys);
	$sys="join  -1 1 -2 1 -t $\'\t\' $dirtemp/input_summary/tmp4 $dirtemp/input_summary/tmp3 >> $dirtemp/input_summary/Gene_Summary_tmp";
	#print $sys."\n";
	system($sys);
	$sys="rm $dirtemp/input_summary/tmp3 $dirtemp/input_summary/tmp4";
	#print $sys."\n";
	system($sys);
	$sys="mv $dirtemp/input_summary/Gene_Summary_tmp $dirtemp/input_summary/Gene_Summary";
	#print $sys."\n";
	system($sys);

}

#reading the covariates and changing the order of sample according to EXPR and plink files
if($covariates ne "")
{
	open(COVAR,"$covariates") or die "not able to read the file $covariates\n";
	open(WRCOVAR,">$dirtemp/processing/covariates.txt") or die "not able to write the file $dirtemp/processing/$covariates\n"; 
	$line=<COVAR>;
	chomp($line);
	@covar=split(/\t/,$line);
	undef(@order);
	@order=();
	for($i=1;$i<@covar;$i++)
	{	
		if(!exists($sampleorder{$covar[$i]}))
		{
			die "sample exist in covariates but not in the expression\n";
		}
		push(@order,$sampleorder{$covar[$i]});	
	}
	
	#writing the covariates file in the temp directory for matrix eqtl
	undef(@array);
	$array[0] = $covar[0];
	for($i=1;$i<@covar;$i++)
	{
		$array[$order[$i-1]] = $covar[$i];	
	}
	$line=join("\t",@array);
	print WRCOVAR "$line\n";
	
	while(<COVAR>)
	{
		chomp($_);
		@covar=split("\t",$_);
		undef(@array);
		$array[0] = $covar[0];
		chomp($array[0]);
		undef(@calc);
		@calc=();
		$na=0;
		for($i=0;$i<@covar;$i++)
		{
			push(@calc,$covar[$i]);
			$array[$order[$i-1]] = $covar[$i];
		}
		$_=join("\t",@array);
		print WRCOVAR "$_\n";
	}
}
#summary ends
print "Copying the summary files to output\n";
system("cp  $dirtemp/input_summary/Tped_Summary $outdir");
system("cp  $dirtemp/input_summary/Gene_Summary $outdir");
#print "summind uc($summ_ind)";
if(uc($summ_ind) eq "YES")
{
	exit 0;
	
}

#preparing file for the EQTL analysis
$sys = "$perl  $dir/bin/perl_convert_snp_genotype.pl $dirtemp";
print $sys."\n";
system($sys);
$sys="cp $expr_loc $dirtemp/processing/geneloc.txt";
print $sys."\n";
system($sys);

#creating script for the main EQTL job
#$sys="mkdir  $dirtemp/temp";
#print "$sys\n";
#system($sys);

open(WRSUB,">$dirtemp/processing/file_qsub")or die "unable to write the file $dirtemp/processing/file_qsub\n";
print WRSUB "$rlib $dirtemp/processing/ $dirtemp/processing/expression.txt $dirtemp/eqtl_run/PVAL.TRANS $cis_pval_in $trans_pval_in $dirtemp/eqtl_run/PVAL.CIS  $cisDist\n";
if($permutation != 0)
{
	$num_permutation=$permutation;
}


for($i=0;$i<$permutation;$i++)
{
	$var = int(time() * rand());
	#$sys="$rscript $dir/bin/rearrange_gene_exp.R $dirtemp/processing/expression.txt $dirtemp/eqtl_run/expression$var.txt";
	#print "$sys\n";
	#system($sys);
	print WRSUB "$rlib $dirtemp/processing/ $dirtemp/eqtl_run/expression$var.txt $dirtemp/eqtl_run/PVAL$var.TRANS $cis_pval_in $trans_pval_in $dirtemp/eqtl_run/PVAL$var.CIS  $cisDist\n";
}
$num_permutation=$permutation+1;
if($filter eq "YES")
{
	system("mkdir $dirtemp/eqtl_run_filter");
}
open(ARRAY_SHAPEIT,">$dirtemp/processing/main_eqtl.csh") or die "no file exists\n";

$com = '#!';
print ARRAY_SHAPEIT "$com /bin/bash\n";
$com = '#$';
if(uc($sge) ne "NO")
 {
	print ARRAY_SHAPEIT "$com -q $impute_queue\n";
	print ARRAY_SHAPEIT "$com -l h_vmem=$impute_mem\n";
	print ARRAY_SHAPEIT "$com -t 1-$num_permutation:1\n";
	print ARRAY_SHAPEIT "$com -M $email\n";
	print ARRAY_SHAPEIT "$com -m a\n";
	print ARRAY_SHAPEIT "$com -V\n";
	print ARRAY_SHAPEIT "$com -cwd\n";
	print ARRAY_SHAPEIT "$com -e $dirtemp/sgelog\n";
	print ARRAY_SHAPEIT "$com -o $dirtemp/sgelog\n";
	print ARRAY_SHAPEIT "$com -notify\n";
}
else
{
	print ARRAY_SHAPEIT 'for SGE_TASK_ID in {1..'.$num_permutation.'}'."\n";
	print ARRAY_SHAPEIT 'do'."\n";
}
$temp = 'k=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f3 -d " "`';
print ARRAY_SHAPEIT "$temp\n";
print ARRAY_SHAPEIT '				if [ !  -f '.'$k'.' ]'."\n";
print ARRAY_SHAPEIT '				then'."\n";
print ARRAY_SHAPEIT "						$rscript $dir/bin/rearrange_gene_exp.R $dirtemp/processing/expression.txt".' $k'."\n";
print ARRAY_SHAPEIT '				fi'."\n";
$temp = 'k=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1`';
print ARRAY_SHAPEIT "$temp\n";
$temp = 'trans_out=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f4 -d " "`';
print ARRAY_SHAPEIT "$temp\n";
$temp = 'cis_out=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f7 -d " "`';
print ARRAY_SHAPEIT "$temp\n";
print ARRAY_SHAPEIT '				if [[ ( !  -f $trans_out ) || ( ! -f $cis_out) ]]'."\n";
print ARRAY_SHAPEIT '				then'."\n";
print ARRAY_SHAPEIT "						$rscript ".$dir.'/bin/Matrix_eQTL_source_cis_ori.R $k'."\n";
print ARRAY_SHAPEIT '				fi'."\n";

#$temp = 'k=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f3 -d " "`';
#print ARRAY_SHAPEIT "$temp\n";
#print ARRAY_SHAPEIT 'rm $k'."\n";
$temp = 'k=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f4 -d " "`';
print ARRAY_SHAPEIT "$temp\n";
#print ARRAY_SHAPEIT 'gzip $k'."\n";
print ARRAY_SHAPEIT 'if [ -f $k ];'."\n";
print ARRAY_SHAPEIT 'then'."\n";
print ARRAY_SHAPEIT '   gzip $k'."\n";
print ARRAY_SHAPEIT 'fi'."\n";

print ARRAY_SHAPEIT '						check_sge=`gunzip -c '.'$k.gz|wc -l`'."\n";
print ARRAY_SHAPEIT '						if [ $check_sge -lt 2 ];'."\n";
print ARRAY_SHAPEIT 'then'."\n";
#print ARRAY_SHAPEIT 'if [ -f $k.gz ];'."\n";
print ARRAY_SHAPEIT '   echo "File $k.gz does not exists or has some problem"'."\n";
print ARRAY_SHAPEIT '   exit 2'."\n";
print ARRAY_SHAPEIT 'fi'."\n";
$temp = 'k1=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f7 -d " "`';
print ARRAY_SHAPEIT "$temp\n";
print ARRAY_SHAPEIT 'if [ -f $k1 ];'."\n";
print ARRAY_SHAPEIT 'then'."\n";
print ARRAY_SHAPEIT '   gzip $k1'."\n";
print ARRAY_SHAPEIT 'fi'."\n";
print ARRAY_SHAPEIT '						check_sge=`gunzip -c '.'$k.gz|wc -l`'."\n";
print ARRAY_SHAPEIT '						if [ $check_sge -lt 2 ];'."\n";
print ARRAY_SHAPEIT 'then'."\n";
print ARRAY_SHAPEIT '   echo "File $k1.gz does not exists or has some problem"'."\n";
print ARRAY_SHAPEIT '   exit 2'."\n";
print ARRAY_SHAPEIT 'fi'."\n";

print ARRAY_SHAPEIT "$perl $dir/bin/perl_100_permute_table_new.pl ".'$k  $k1'." $dirtemp/eqtl_run/file_table".'$SGE_TASK_ID'." $cis_pval_in\n";
$temp = 'k=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f4 -d " "|rev|cut -f1 -d "/"|rev`';
print ARRAY_SHAPEIT "$temp\n";
$temp = 'k1=`cat  '."$dirtemp/processing/file_qsub".' |head -$SGE_TASK_ID |tail -1|cut -f7 -d " "|rev|cut -f1 -d "/"|rev`';
print ARRAY_SHAPEIT "$temp\n";
if($filter eq "YES")
{
	$temp = "$perl $dir/bin/perl_filter_main_permutation.pl $dirtemp/input_summary/Tped_Summary $dirtemp/input_summary/Gene_Summary $dirtemp/eqtl_run/".'$k.gz '." $dirtemp/eqtl_run_filter/".'$k.gz '."$filter_geno_count $filter_geno_na $filter_geno_maf $filter_gene_mean_max $filter_gene_mean_min  $filter_gene_sd_min $filter_gene_sd_max $filter_gene_na\n";
	print ARRAY_SHAPEIT "$temp\n";
	$temp = "$perl $dir/bin/perl_filter_main_permutation.pl $dirtemp/input_summary/Tped_Summary $dirtemp/input_summary/Gene_Summary $dirtemp/eqtl_run/".'$k1.gz '." $dirtemp/eqtl_run_filter/".'$k1.gz '."$filter_geno_count $filter_geno_na $filter_geno_maf $filter_gene_mean_max $filter_gene_mean_min  $filter_gene_sd_min $filter_gene_sd_max $filter_gene_na\n";
	print ARRAY_SHAPEIT "$temp\n";
	print ARRAY_SHAPEIT "$perl $dir/bin/perl_100_permute_table_new.pl $dirtemp/eqtl_run_filter/".'$k '."$dirtemp/eqtl_run_filter/".'$k1'." $dirtemp/eqtl_run_filter/file_table".'$SGE_TASK_ID'." $cis_pval_in\n";
}
if(uc($sge) eq "NO")
{
	print ARRAY_SHAPEIT "done\n";
}
close(ARRAY_SHAPEIT);
if(uc($sge) eq "NO")
{
	system("$sh $dirtemp/processing/main_eqtl.csh");
}
 else
{
	system("$qsub $dirtemp/processing/main_eqtl.csh > $dirtemp/processing/jobid_shapeit");
	#readin job id from submit_shapeit
	open(ARRAY_SHAPEIT,"$dirtemp/processing/jobid_shapeit") or die "unable to open file $dirtemp/processing/jobid_shapeit\n";
	$shapeit = <ARRAY_SHAPEIT>;
	print "$shapeit\n";
	@shapeit =split(" ",$shapeit);
	@shapeit1 =split(/\./,$shapeit[2]);
	print "JOB ID $shapeit1[0]\n";
	$job_id_shapeit = $shapeit1[0];
	#creating script that makes the main script
	open(ARRAY_SHAPEIT,">$dirtemp/processing/ArrayJob_shapeit_wait.csh") or die "unable to create the array job wait file shape it \n";
	print ARRAY_SHAPEIT '#! /bin/bash'."\n";
	print ARRAY_SHAPEIT "$com -q $impute_queue\n";
	print ARRAY_SHAPEIT "$com -l h_vmem=$impute_mem\n";
	print ARRAY_SHAPEIT '#$ -M '."$email\n";
	print ARRAY_SHAPEIT '#$ -m a'."\n";
	print ARRAY_SHAPEIT '#$ -V'."\n";
	print ARRAY_SHAPEIT '#$ -cwd'."\n";
	$com = '#$';
	print ARRAY_SHAPEIT "$com -e $dirtemp/sgelog\n";
	print ARRAY_SHAPEIT "$com -o $dirtemp/sgelog\n";
	print ARRAY_SHAPEIT "cp $dirtemp/processing/jobid_shapeit $dirtemp/processing/waiting.txt\n";
	$sys = "$qsub -hold_jid $job_id_shapeit $dirtemp/processing/ArrayJob_shapeit_wait.csh\n";
	print $sys."\n";
	system($sys);


	$flag = 1;
	while($flag > 0)
	{
		#print "in the loop\n";
		if($flag%100000000 == 0)
		{
			$flag = 1;
			print "waiting for the array job to complete : $countit\n";
		}
		$flag++;
		#$countit=`wc -l $dirtemp/file_table`;
		#chomp($countit);
		#@countit = split(" ",$countit);
		#if($countit[0] eq "202")
		if(-e "$dirtemp/processing/waiting.txt")
		{
			$flag =0;
		}
	}
}

#combing the count files in to one
system("cat $dirtemp/eqtl_run/file_table* >$dirtemp/processing/file_table");
if($permutation > 1)
{
	system("rm  $dirtemp/eqtl_run/expression*");
}
#system("rm  $dirtemp/eqtl_run/file_table*");

#separating the cis and trans counts
open(BUFF,"$dirtemp/processing/file_table") or die " no file exists $dirtemp/processing/file_table\n";
open(WRBUFF1,">$dirtemp/processing/results_100permutation_cis_numbers") or die "not able write $dirtemp/processing/results_100permutation_cis_numbers\n";
open(WRBUFF2,">$dirtemp/processing/results_100permutation_trans_numbers") or die "not able write $dirtemp/processing/results_100permutation_trans_numbers\n";
while(<BUFF>)
{
        print WRBUFF1 $_;
        $line2 = <BUFF>;
        print WRBUFF2 $line2;
}
close(WRBUFF1);
close(WRBUFF2);
$count_cs = `head -1 $dirtemp/processing/results_100permutation_cis_numbers|wc -w`;
chomp($count_cs);
#sorting the ci and trans count files
system("sort -k$count_cs,$count_cs $dirtemp/processing/results_100permutation_cis_numbers > $dirtemp/processing/results_100permutation_cis_numbers1");
system("mv $dirtemp/processing/results_100permutation_cis_numbers1 $dirtemp/processing/results_100permutation_cis_numbers");
system("sort -k$count_cs,$count_cs $dirtemp/processing/results_100permutation_trans_numbers > $dirtemp/processing/results_100permutation_trans_numbers1");
system("mv $dirtemp/processing/results_100permutation_trans_numbers1 $dirtemp/processing/results_100permutation_trans_numbers");

$permutation=`wc -l $dirtemp/processing/results_100permutation_trans_numbers`;
chomp($permutation);
$permutation=~ s/ .+//g;
$permutation=$permutation-1;

#creating the table from cis and trans count files
printreport("$dirtemp/processing/results_100permutation_cis_numbers","$dirtemp/processing/results_100permutation_trans_numbers","$dirtemp/processing/Table_nofiltering.txt",$cis_pval_in,$permutation);
if($filter eq "YES")
{
	#creating report fo r filtered results
	system("cat $dirtemp/eqtl_run_filter/file_table* >$dirtemp/processing/file_table_filter");
	system("rm $dirtemp/eqtl_run_filter/file_table* ");
	#separating the cis and trans filtered files
	open(BUFF,"$dirtemp/processing/file_table_filter") or die " no file exists $dirtemp/file_table\n";
	open(WRBUFF1,">$dirtemp/processing/results_100permutation_cis_numbers_filter") or die "not able write $dirtemp/processing/results_100permutation_cis_numbers\n";
	open(WRBUFF2,">$dirtemp/processing/results_100permutation_trans_numbers_filter") or die "not able write $dirtemp/processing/results_100permutation_trans_numbers\n";
	while(<BUFF>)
	{
			print WRBUFF1 $_;
			$line2 = <BUFF>;
			print WRBUFF2 $line2;
	}
	close(WRBUFF1);
	close(WRBUFF2);
	#sorting the cis and trans filtered count files
	system("sort -k$count_cs,$count_cs $dirtemp/processing/results_100permutation_cis_numbers_filter > $dirtemp/processing/results_100permutation_cis_numbers1_filter");
	system("mv $dirtemp/processing/results_100permutation_cis_numbers1_filter $dirtemp/processing/results_100permutation_cis_numbers_filter");
	system("sort -k$count_cs,$count_cs $dirtemp/processing/results_100permutation_trans_numbers_filter > $dirtemp/processing/results_100permutation_trans_numbers1_filter");
	system("mv $dirtemp/processing/results_100permutation_trans_numbers1_filter $dirtemp/processing/results_100permutation_trans_numbers_filter");
	#creating the table from cis and trans filtered count files
	printreport("$dirtemp/processing/results_100permutation_cis_numbers_filter","$dirtemp/processing/results_100permutation_trans_numbers_filter","$dirtemp/processing/Table_withfiltering.txt",$cis_pval_in,$permutation);
}
print "Copying the output files to output\n";
system("cp  $dirtemp/processing/Table_nofiltering.txt $outdir");
system("cp  $dirtemp/eqtl_run/PVAL.CIS.gz $outdir");
$count=0;
$count=`gunzip -c $dirtemp/eqtl_run/PVAL.CIS.gz|wc -l`;
chomp($count);
if($count > 2)
{
	system("$rscript $dir/bin/qqplot.R $outdir PVAL.CIS.gz PVAL.CIS.png");
}
system("cp  $dirtemp/eqtl_run/PVAL.TRANS.gz $outdir");
$count=0;
$count=`gunzip -c $dirtemp/eqtl_run/PVAL.TRANS.gz|wc -l`;
chomp($count);
if($count > 2)
{
	system("$rscript $dir/bin/qqplot.R $outdir PVAL.TRANS.gz PVAL.TRANS.png");
}
if($filter eq "YES")
{
	system("cp $dirtemp/processing/Table_withfiltering.txt $outdir");
	system("cp $dirtemp/eqtl_run_filter/PVAL.CIS.gz $outdir/PVAL_FILTER.CIS.gz");
	$count=0;
	$count=`gunzip -c $dirtemp/eqtl_run_filter/PVAL.CIS.gz|wc -l`;
	chomp($count);
	if($count > 2)
	{
		system("$rscript $dir/bin/qqplot.R $outdir PVAL_FILTER.CIS.gz PVAL_FILTER.CIS.png");
	}

	system("cp $dirtemp/eqtl_run_filter/PVAL.TRANS.gz $outdir/PVAL_FILTER.TRANS.gz");
	$count=0;
	$count=`gunzip -c $dirtemp/eqtl_run_filter/PVAL.TRANS.gz|wc -l`;
	chomp($count);
	if($count > 2)
	{
		system("$rscript $dir/bin/qqplot.R $outdir PVAL_FILTER.TRANS.gz PVAL_FILTER.TRANS.png");
	}
}
###########################################SUBROUTINES##########################################################
sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}
sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub printreport
{
	my($infile_cis,$infile_trans,$outfile,$cis_pval_in,$permutation) =@_;
	open(BUFF1,$infile_cis) or die " no file exists $infile_cis\n";
	open(BUFF2,$infile_trans) or die " no file exists $infile_trans\n";
	open(WRBUFF,">$outfile") or die "no file $outfile exists\n";
	$n=0;
	for($i=$cis_pval_in;$i> 1e-74;$i = $i/2)
	{
		push(@temp,$i);
		$n++;
	}
	@temp1 = reverse(@temp);
	@temp = @temp1;
	undef(@temp1);
	for($i=0;$i<$n;$i++)
	{
		$table[0][$i]= $temp[$i];
	}
	$table[0][$n]='FILENAME';
	$n++;
	#die "$n\n";
	$i=1;
	while(<BUFF1>)
	{
		chomp($_);
		@array=split("\t",$_);
		for($j=0;$j<@array;$j++)
		{
			$table[$i][$j]=$array[$j];
		}
		$i++;
	}
	while(<BUFF2>)
	{
		chomp($_);
		@array=split("\t",$_);
		for($j=0;$j<@array;$j++)
		{
			$table[$i][$j]=$array[$j];
		}
		$i++;
	}
	for($k=0;$k<$i;$k++)
	{
		for($j=0;$j<$n;$j++)
		{
			$table1[$j][$k] = $table[$k][$j];
		}
	}
	#print header
	print WRBUFF "Pval-cutoff";
	for($z=0;$z<$permutation;$z++)
	{
		$tm_num=$z+1;
		print WRBUFF "\tnum_cis_eqtl_perm_$tm_num";
	}
	print WRBUFF "\tnum_cis_eqtl_main_run";
	for($z=0;$z<$permutation;$z++)
	{
		$tm_num=$z+1;
		print WRBUFF "\tnum_trans_eqtl_perm_$tm_num";
	}
	print WRBUFF "\tnum_trans_eqtl_perm_$tm_num";
	print WRBUFF "\tcis_FDR";
	print WRBUFF "\ttrans_FDR";
	print WRBUFF "\n";
	#print header end
	$temp=$table1[$n-1][0];
	for($j=1;$j<$i;$j++)
	{
		$temp=$temp."\t".$table1[$n-1][$j];
	}
	$temp="$temp\tFDR\tFDR";
	print WRBUFF  $temp."\n";
	$permutation_1=$permutation+1;
	$permutation_2=$permutation+2;
	$permutation_22=2*$permutation_1;
	for($k=0;$k<$n-1;$k++)
	{
		$temp=$table1[$k][0];
		$sum_cis=0;
		$sum_trans=0;
		for($j=1;$j<$permutation_1;$j++)
		{
			$sum_cis=$sum_cis+$table1[$k][$j];	
		}
		for($j=$permutation_2;$j<$permutation_22;$j++)
		{
			$sum_trans=$sum_trans+$table1[$k][$j];	
		}
		for($j=1;$j<$i;$j++)
		{
			$temp=$temp."\t".$table1[$k][$j];
		}
		#$sum_cis=$sum_cis/100;
		#$sum_trans=$sum_trans/100;
		if($table1[$k][$permutation_1] !=0)
		{
			#$sum_cis=100*($sum_cis/$table1[$k][102]);
			$sum_cis=$sum_cis/$table1[$k][$permutation_1];	
			if($permutation != 0)
			{
				$sum_cis=100*$sum_cis/$permutation;
			}
			else
			{
				$sum_cis = 0;
			}
		}	
		else
		{
			$sum_cis=0;
		}
		if($table1[$k][$permutation_22] !=0)
		{
			#$sum_trans=100*($sum_trans/$table1[$k][203]);
			$sum_trans=$sum_trans/$table1[$k][$permutation_22];
			
			if($permutation != 0)
			{
				$sum_trans=100*$sum_trans/$permutation;
			}
			else
			{
				$sum_trans = 0;
			}
		}
		else
		{
			$sum_trans=0;
		}
		print WRBUFF $temp."\t$sum_cis\t$sum_trans\n";
	}
}
