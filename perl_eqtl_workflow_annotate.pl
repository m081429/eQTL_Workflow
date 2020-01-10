
#!/usr/bin/perl
#  This script is used to get the annotations using annovar

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
	
	die "Retry with : perl  perl_eqtl_workflow_annotate.pl  -run_config <PATH TO THE RUN CONFIG FILE> -tool_config <PATH TO THE TOOL CONFIG FILE>\n";
}

require "$dir/bin/CONFIG.pl";
#reading config info parameters
getDetails($config);


my $list= $config{"snp.list"};
my $dirtemp = $config{"temp.folder"};
my $tped = $config{"tped"};
my $tfam = $config{"tfam"};
my $build = $config{"build"};
my $annot = $config{"annot"};
my $ann_db = $config{"db"};
$list =~ s/\s|\t|\r|\n//g;
$dirtemp =~ s/\s|\t|\r|\n//g;
$tped =~ s/\s|\t|\r|\n//g;
$tfam =~ s/\s|\t|\r|\n//g;
$build =~ s/\s|\t|\r|\n//g;
$annot =~ s/\s|\t|\r|\n//g;
$ann_db =~ s/\s|\t|\r|\n//g;

#reading tool config info parameters
getDetails($toolconfig);
my $plink= $config{"path_plink"};
my $perl=$config{"path_perl"};
my $annovar=$config{"path_annovar"};

$perl =~ s/\s|\t|\r|\n//g;
$plink =~ s/\s|\t|\r|\n//g;
$annovar =~ s/\s|\t|\r|\n//g;
print "***********RUN INFO CONFIG ARGUMENTS***********\n";
print "SNP_LIST: $list\n";
print "TEMP_FOLDER: $dirtemp\n";
print "TPED : $tped\n";
print "TFAM : $tfam\n";  
print "BUILD :hg$build\n";
print "ANNOT : $annot\n";
print "DB : $ann_db\n";

print "***********TOOLINFO CONFIG ARGUMENTS***********\n";
print "PERL : $perl\n";
print "PLINK : $plink\n";
print "ANNOVAR : $annovar\n";

if( $list !~ m/\w/ | $dirtemp !~ m/\w/ |$tped !~ m/\w/ | $tfam !~ m/\w/ |$build !~ m/\w/ |$annot !~ m/\w/  |$ann_db !~ m/\w/)
{
	print "**********************PARAMETERS TO BE IN THE CONFIG FILE***********************************\n";
	print "(1) SNP_LIST:List of RSIDs to be annotated\n";
	print "(2) TEMP_FOLDER: PATH to the TEMP FOLDER\n";
	print "(3) TPED : PATH to the TPED  FILE\n";
	print "(4) TFAM : PATH to the TFAM  FILE\n";
	print "(5) BUILD :Build Version 18/19\n";
	print "(6) ANNOT :Multiple annotations can be separated by ',' and below are the different annotations available\n";
	print "GENE_ANNOT_REFGENE-> 'geneanno__gene'\n";
	#print "GENE_ANNOT_UCSC-> 'geneanno__knowngene'\n";
	#print "GENE_ANNOT_ENSEMBL-> 'geneanno__ensgene'\n";
	print "REGION_BASED_ANNOT_MOST_CONSERVED_ELEMENT-> 'regionanno__mce44way'\n";
	print "REGION_BASED_ANNOT_TRANSCRIPT_FACTOR_BINDING_SITE-> 'regionanno__tfbs'\n";
	print "REGION_BASED_ANNOT_CYTO_BAND-> 'regionanno__band'\n";
	print "REGION_BASED_ANNOT_SEG_DUP-> 'regionanno__segdup'\n";
	print "REGION_BASED_ANNOT_STRUCT_VAR_DGV-> 'regionanno__dgvMerged'\n";
	print "REGION_BASED_ANNOT_PUBLISHED_GWAS-> 'regionanno__gwascatalog'\n";
	print "FILTER_BASED_ANNOT_1000genome-> 'filter__1000g2012apr'\n";
	print "FILTER_BASED_ANNOT_dbSNP annotations-> 'filter__snp130__webfrom'\n";
	print "FILTER_BASED_ANNOT_SIFT_POLYPHEN annotations-> 'filter__avsift__webfrom'\n"; 
	print "FILTER_BASED_ANNOT_PolyPhen2, MutationTaster, LRT, PhyloP annotations-> 'filter__ljb_pp2__webfrom'\n"; 
	print "FILTER_BASED_ANNOT_exome sequencing project annotations-> 'filter__esp6500_all__webfrom'\n";
	print "FILTER_BASED_ANNOT_GERP++ annotations-> 'filter__gerp++gt2__webfrom'\n";
	die "some of the input parameters are Empty.Retry with new parameters\nperl perl_eqtl_workflow_annotate.pl -run_config run.config -tool_config tool.config\n";
}
if(!(-e $perl))
{
	die "specified perl $perl does not exist\n";
}
if(!(-e $plink))
{
	die "specified plink $plink does not exist\n";
}
if(!(-e "$annovar/annotate_variation.pl"))
{
	die "specified ANNOVAR direcoty if wrong does not exist\n";
}
unless(-d $dirtemp)
{
    system("mkdir -p $dirtemp");
}
open(WRLOG,">$dirtemp/FINAL.log") or die "no file found $dirtemp/FINAL.log\n";
#preparing DB for all the annotation
print WRLOG "preparing DB for annotation module\n";
unless(-d $ann_db)
{
    system("mkdir -p $ann_db");
}
print WRLOG "downloading reference to figure out reference allele\n";
unless(-d "$ann_db/UCSC_REF_HG$build")
{
    system("mkdir -p $ann_db/UCSC_REF_HG$build");
	chdir("$ann_db/UCSC_REF_HG$build");
	
	@array=('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','M');
	for($i=0;$i<@array;$i++)
	{
		 print WRLOG "downloading with chr$array[$i]\n";
		if( !(-e "$ann_db/UCSC_REF_HG$build/chr$array[$i].fa.gz"))
		{	
			$sys = "wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg$build/chromosomes/chr$array[$i].fa.gz";
			system($sys);	
		}
	}

}
print WRLOG "Done preparing DB for annotation module\n";
my @ls=();
my @annotation=();
my @annotation_fil=();
print WRLOG "Downloading annotation DB for annotation module\n";
print WRLOG "downloading reference for annotations\n";
unless(-d "$ann_db/ANN_HG$build")
{
    system("mkdir -p $ann_db/ANN_HG$build");
}	
chdir("$ann_db/ANN_HG$build");
@annotations =split(',',$annot);

for($i=0;$i<@annotations;$i++)
{
	print WRLOG "Downloading files for $annotations[$i] request\n";
	$ann_dir = "$ann_db/ANN_HG$build/ann_$annotations[$i]";
	unless(-d $ann_dir)
	{
		system("mkdir $ann_dir");
	}	
	chdir($ann_dir) or die "can't change directory $ann_dir\n";
	@annot=split("__",$annotations[$i]);
	if($annot[2] eq "webfrom")
	{
		$annot[2] ="-webfrom annovar";
	}
	else
	{
		$annot[2] ="";
	}
	$ls=`ls $ann_dir`;
	chomp($ls);
	undef(@ls);
	@ls=split("\n",$ls);
	if(@ls < 2 || !(-e "$ann_dir/sucess"))
	{	
		#print "downloading required files\n";
		#$date=`date`;
		#chomp($date);
		#$date =~ s/ /_/g;
		print "executing remove command .Just to do clean up in case any file are present before download\n";
		system("rm $ann_dir/*");
		$sys = "$annovar/annotate_variation.pl -downdb $annot[1] . $annot[2] -buildver hg$build \n";
		print WRLOG "executing $sys\n";
		system($sys);
	}
	$ls=`ls $ann_dir`;
	chomp($ls);
	undef(@ls);
	@ls=split("\n",$ls);
	if(@ls < 2)
	{
		print WRLOG "\n\n\n*****download failed for $annotations[$i] directory $ann_dir****\nUpdate the version and retry.";
		print WRLOG "******removing request for this annotation********\n\n\n\n\n\n\n";
		system("rm $ann_dir/*");
		system("rmdir $ann_dir/*");
	}
	else
	{
		print WRLOG "\n\n\n*****download succes for $annotations[$i] directory $ann_dir****\nUpdate the version and retry.";
		print WRLOG "******including request for this annotation********\n\n\n\n\n\n\n";
		system("touch $ann_dir/sucess");
		push(@annotation_fil,$annotations[$i]);
	}
}
if(@annotation_fil == 0)
{
	die "no annotations avail;able after filter.Please try with new annotations request\n";
}
#die;
print WRLOG "Done downloading annotation DB for annotation module\n";
print WRLOG "Done preparing DB for annotation module\n";	
print WRLOG "Done preparing DB for all the annotation\n";

print WRLOG "creating the temp directory if not exists\n";

print WRLOG "copying the SNP list to the temp directory\n";
system("cp $list $dirtemp/snplist.txt");

print WRLOG "copying the Plink files to the temp directory\n";
system("cp $tped $dirtemp/input.tped");
system("cp $tfam $dirtemp/input.tfam");

print WRLOG "extracting the list of snps from the dataset to convertto bim file\n";
system("$plink --tfile $dirtemp/input --extract $dirtemp/snplist.txt --make-bed --out $dirtemp/extracted");  
system("rm $dirtemp/extracted.bed");
system("rm $dirtemp/extracted.fam");


print WRLOG "extracting the positions from the UCSC reference and switcing the alleles\n";
open(BIM,"$dirtemp/extracted.bim") or die "no file found $dirtemp/extracted.bim\n";		
open(WRBIM,">$dirtemp/input_prep.bim") or die "no file found $dirtemp/input_prep.bim\n";
$prev_chr=0;
$num_fasta=50;
while(<BIM>)
{
	chomp($_);
	@bim=split("\t",$_);
	if($prev_chr ne $bim[0])
	{
		$chr=$bim[0];
		$chr =~ s/23/X/g;
		$chr =~ s/24/Y/g;
		$chr =~ s/26/M/g;
		#print $chr."\n";
		#$sys = "wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg$build/chromosomes/chr$chr.fa.gz";
		#print $sys."\n";
		#system($sys);
		open(CHR,"gunzip -c $ann_db/UCSC_REF_HG$build/chr$chr.fa.gz |") or die " no file found $ann_db/UCSC_REF_HG$build/chr$chr.fa.gz\n";
		$line = <CHR>;
		$length =0;
		$prevlength =0;
	}
	while($bim[3] > $length)
	{
		if($line=<CHR>)
		{
			chomp($line);
			$prevlength = $length;
			$length=$length+length($line);
		}
		else
		{
			die "$_ line exceeds current chromsome position $length\n";
		}
	}
	@current=split("",$line);
	$current=$bim[3]-$prevlength-1;
	$current[$current]=uc($current[$current]);
	$bim[4]=uc($bim[4]);
	$bim[5]=uc($bim[5]);
	if($current[$current] eq $bim[4] || $current[$current] eq $bim[5])
	{
		#print "sucess $_ $current[$current]\n";
		if($current[$current] eq $bim[4])
		{
			print WRBIM "$bim[0] $bim[3] $bim[3] $bim[4] $bim[5] comments: $bim[1]\n";
		}
		else
		{
			print WRBIM "$bim[0] $bim[3] $bim[3] $bim[5] $bim[4] comments: $bim[1]\n";
		}
	}
	else
	{
		print WRBIM "$bim[0] $bim[3] $bim[3] $bim[5] $bim[4] comments: $bim[1]\n";
		print WRLOG "warning $_ as not matching reference position $current[$current]\n";
	}
	#die "$prevlength $length\n";
	$prev_chr=$bim[0];	
}
close(BIM);
close(WRBIM);
chdir("$dirtemp") or die "can't change directory $dirtemp\n";
#system("rm $dirtemp/UCSC_refdir/*");
#system("rmdir $dirtemp/UCSC_refdir");
system("rm $dirtemp/input.tfam");
system("rm $dirtemp/input.tped");

#@annotations =split(',',$annot);
print WRLOG "preparing FINAL.txt file\n";
open(FINAL,"$dirtemp/input_prep.bim") or die "no file found $dirtemp/input_prep.bim\n";
open(WRFINAL,">$dirtemp/FINAL.txt") or die "no file found 	$dirtemp/FINAL.txt\n";
print WRFINAL "chr start_pos stop_pos ALL1 ALL2 comments: marker\n";
while(<FINAL>)
{
	print WRFINAL $_;
}
close(WRFINAL);
close(FINAL);
for($i=0;$i<@annotation_fil;$i++)
{
	print WRLOG "processing request for $annotations[$i]\n";
	$ann_dir = "$dirtemp/ann_$annotation_fil[$i]";	
	system("mkdir $ann_dir");
	chdir($ann_dir) or die "can't change directory $ann_dir\n";
	@annot=split("__",$annotation_fil[$i]);
	if($annot[2] eq "webfrom")
	{
		$annot[2] ="-webfrom annovar";
	}
	else
	{
		$annot[2] ="";
	}
	#print "downloading required files\n";
	#$sys = "$perl $dir/bin/ANNOVAR_SCRIPTS/annotate_variation.pl -downdb $annot[1] . $annot[2] -buildver hg$build\n";
	#print "$sys\n";
	#system($sys);
	print WRLOG "Annotating $annotation_fil[$i]...\n";
	system("cp $dirtemp/input_prep.bim $ann_dir");
	$sys ="$perl $annovar/annotate_variation.pl -buildver hg$build -$annot[0] -dbtype $annot[1] $ann_dir/input_prep.bim $ann_db/ANN_HG$build/ann_$annotation_fil[$i] ";
	print WRLOG "$sys\n";
	#die;
	system($sys);
	#die;
	
	system("rm $ann_dir/input_prep.bim");	
	$ls=`ls $ann_dir`;
	@ls=split("\n",$ls);
	if(-e "$ann_dir/input_prep.bim.invalid_input")
	{
		print WRLOG "Invalid Input file found as annotation failed.So skipping this annotation file in the FINAL.txt\n";
		next;
	}
	for($f=0;$f<@ls;$f++)
	{
		if($ls[$f] !~ m/log$/ && $ls[$f] !~ m/FINAL.txt/ && $ls[$f] !~ m/_filtered$/)
		{
			#die "$ls[$f]\n";
			open(FINAL,"$dirtemp/FINAL.txt") or die "no file found $dirtemp/FINAL.txt\n";
			open(BUF1,"$ann_dir/$ls[$f]") or die "no file found $ann_dir/$ls[$f]\n";
			open(WRBUF1,">$ann_dir/TEMP.txt") or die "no no file found $ann_dir/TEMP.txt\n";
			while(<BUF1>)
			{
				chomp($_);
				@buf=split("\t",$_);
				$buf=pop(@buf);
				print WRBUF1 $buf."\t".$_."\n";
			}
			close(BUF1);
			close(WRBUF1);
			$sys="sort -k1,1n -k2,2n $ann_dir/TEMP.txt > $ann_dir/$ls[$f]";
			system($sys);
			system("rm $ann_dir/TEMP.txt");
			open(TEMP,"$ann_dir/$ls[$f]") or die "no file found $ann_dir/$ls[$f]\n";
			open(WRFINAL,">$dirtemp/FINAL1.txt") or die "no file found $dirtemp/FINAL1.txt \n";
			$final=<FINAL>;
			chomp($final);
			$ls[$f]=~s/input_prep\.bim\.//g;
			print WRFINAL "$final\t$ls[$f]\n";
			$temp=<TEMP>;
			@temp=split("\t",$temp);
			shift(@temp);
			if($temp[0] =~ m/^line/)
			{
				shift(@temp);
			}
			$info=pop(@temp);
			$anno_temp=join('_',@temp);
			$anno_temp =~ s/\t/ /g;
			@ini=split(" ",$info);
			$chr=shift(@ini);
			$pos=shift(@ini);
			while($final=<FINAL>)
			{
				chomp($final);
				@final=split(" ",$final);
				$final_chr=shift(@final);
				$final_pos=shift(@final);
				if("$final_chr.$final_pos" eq "$chr.$pos" && $chr ne "NA" && $pos ne "NA")
				{
					print WRFINAL $final."\t";
					$flag=0;
					while("$final_chr.$final_pos" eq "$chr.$pos" && $chr ne "NA" && $pos ne "NA")
					{
						if($flag==0)
						{
							print WRFINAL "$anno_temp";
							$flag++;
						}
						else
						{
							print WRFINAL "___$anno_temp";
						}
						$temp=<TEMP>;
						if($temp =~ m/\w/)
						{
							@temp=split("\t",$temp);
							@temp=split("\t",$temp);
							shift(@temp);
							if($temp[0] =~ m/^line/)
							{
								shift(@temp);
							}
							$info=pop(@temp);
							$anno_temp=join("_",@temp);
							$anno_temp =~ s/\t/ /g;
							@ini=split(" ",$info);
							$chr=shift(@ini);
							$pos=shift(@ini);
						}
						else
						{
							$chr="NA";
							$pos="NA";
						}
					}
					print WRFINAL "\n";
				}
				else
				{
					print WRFINAL $final."\tNA\n";	
				}
			}
			close(WRFINAL);
			close(FINAL);
			close(TEMP);
			system("mv $dirtemp/FINAL1.txt $dirtemp/FINAL.txt");
		}
	}
	print "\n\nsee $ann_dir for results for $annotation_fil[$i] annotations\n\n\n\n\n\n\n\n\n\n\n";
	
}
#die "@annotations\n";
close(WRLOG);