#!/usr/bin/perl
# This script is used to calculate the LD blocks

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
	
	die "Retry with : perl  perl_eqtl_cal_LD.pl  -run_config <PATH TO THE RUN CONFIG FILE> -tool_config <PATH TO THE TOOL CONFIG FILE>\n";
}

require "$dir/bin/CONFIG.pl";
#reading config info parameters
getDetails($config);


my $list= $config{"snp.list"};
my $dirtemp = $config{"temp.folder"};
my $ld_block=$config{"ld.window.size"};
my $ld_window_r2=$config{"ld.window.r2"};
my $rounded = $config{"run.dir"};
my $outdir = $config{"output.dir"};
my $outfile = $config{"outputfile"};

#reading tool config info parameters
getDetails($toolconfig);
my $plink=$config{"path_plink"};

$list =~ s/\s|\t|\r|\n//g;
$dirtemp =~ s/\s|\t|\r|\n//g;
$ld_block =~ s/\s|\t|\r|\n//g;
$ld_window_r2 =~ s/\s|\t|\r|\n//g;
$outdir =~ s/\s|\t|\r|\n//g;
$plink =~ s/\s|\t|\r|\n//g;
$rounded =~ s/\s|\t|\r|\n//g;
$outfile =~ s/\s|\t|\r|\n//g;

print "***********INPUT ARGUMENTS***********\n";
print "snp.list: $list\n";
print "temp.folder: $dirtemp\n";
print "run.dir : $rounded\n";
print "ld.window.size : $ld_block\n";
print "ld.window.r2 : $ld_window_r2\n";
print "output.dir : $outdir\n";
print "outputfile : $outfile\n";
print "path_plink : $plink\n";

if(!(-e $plink))
{
	die "Specified path Plink not exists\n";
}

if(!(-e $list))
{
	die "Input file $list not exists\n";
}

if((!(-e "$dirtemp/$rounded/processing/processed_inputdata.tped")) || (!(-e "$dirtemp/$rounded/processing/processed_inputdata.tfam")))
{
	die "Plink TPED file $dirtemp/$rounded/processing/processed_inputdata.tped or Plink TFAM file $dirtemp/$rounded/processing/processed_inputdata.tfam not found\n";
}

if( $plink !~ m/\w/ |$list !~ m/\w/ | $dirtemp !~ m/\w/ |$outfile !~ m/\w/ | $outdir !~ m/\w/ |$ld_block !~ m/\w/  | $ld_window_r2 !~ m/\w/ | $rounded !~ m/\w/)
{
	
	die "Some of the input parameters are empty\n";
}



#creating the temp directory if not exists
unless(-d "$outdir/$rounded")
{
    system("mkdir -p $outdir/$rounded");
}

#extracting the list of snps from the dataset to extract bim file
$sys = "$plink --tfile $dirtemp/$rounded/processing/processed_inputdata  --extract $list --blocks --ld-window $ld_block --ld-window-r2 $ld_window_r2 --out $outdir/$rounded/$outfile";
print "executing : $sys\n";
system($sys);
print "Output files : $outdir/$rounded/$outfile.blocks,$outdir/$rounded/$outfile.blocks.det are generated if no errors occured\n"; 