#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use Cwd 'abs_path';

my ($test_hit, $ref_hit);
my $log2 = 0.6;
my $pvalue = 0.00001;
my $min_window=3;
my $bigger = 1.5;
my $Rexe = '/usr/bin/env R';
my $cw='TRUE';
my $hmm = 'TRUE';
my $cbs = 'FALSE';

my $refine = 0;

my $usage = <<"USAGE";

	usage: $0 [options]

		### for stage 1 ###

		--test = test.hits.file
		   (sorted, only one chromosome)
		--ref = ref.hits.file
		   (sorted, only one chromosome)
		
		--log2-threshold = number
		  (default=$log2)
		--p-value = number
		  (default=$pvalue)
		--bigger-window = number
		  (default=$bigger)
		--minimum-windows-required = number 
		  (default=$min_window)

		--cw    (default; simple consequtive window annotation)
		--no-cw
		--hmm   (default; HMM segmentation)
		--no-hmm
		--cbs     (CBS segmentation)
		--no-cbs  (default)

		### for stage 2 ###
		
		--refine = stage.1.log
		  (all other options will be overwritten)

		### others ###
		--myR = path to your R program
		  (default: "$Rexe")
		--help

USAGE

GetOptions(
	"test=s" => \$test_hit,
	"ref=s" => \$ref_hit,
	"log2-threshold=f" => \$log2,
	"p-value=f" => \$pvalue,
	"bigger-window=f" => \$bigger,
	"no-cw" => sub{$cw='FALSE'},
	"cbs" => sub{$cbs='TRUE'},
	"no-hmm" => sub{$hmm='FALSE'},
	"myR=s" => \$Rexe,
	"minimum-windows-required=i" => \$min_window,

	"refine=s" => \$refine,

	"help|?" => sub{print $usage; exit}
);

die $usage unless $test_hit and $ref_hit;

my $seghmm_path = abs_path($0);
$seghmm_path =~ s/\/[^\\\/]*$//;

my $return = system(qq(echo 'source("$seghmm_path/lib.segHMM.r"); check.library(hmm = $hmm, cbs = $cbs, refine = "$refine");q()' | $Rexe --vanilla --slave 2>/dev/null));
die "missing required R package\n" if $return;

if($refine)
{
	open(SL, $refine) or die "can not open --refine (stage 1 log file: $refine)\n";
	my $log = join('', <SL>);

	my $outfile = $1 if $log =~ m/rough HMM output:\t(.*)$/m;
	my $rough_file = $outfile;
	$outfile =~ s/(\.rough|)$/.refined/;
	my $unreliablefile = $outfile;
	$unreliablefile =~ s/(\.refined)$/.unreliable/;

	my $test_hit = $1 if $log =~ m/test hit file:\t(.*)$/m;
	my $ref_hit = $1 if $log =~ m/ref hit file:\t(.*)$/m;

	die "can not find file $test_hit (from $refine)\n" unless -f $test_hit;
	die "can not find file $ref_hit (from $refine)\n" unless -f $ref_hit;

	stage2($test_hit, $ref_hit, $rough_file, $outfile, $unreliablefile, $seghmm_path);

	open(SL, ">>$refine") or die;
	print SL "raw reads around CNV candidates:\t$outfile.cand\n";
	print SL "unreliable HMM output:\t$unreliablefile\n";
	print SL "refined HMM output:\t$outfile\n";
	close(SL);
	exit;
}

die $usage unless($test_hit && $ref_hit);
die "can not find $test_hit\n" unless -f $test_hit;
die "can not find $ref_hit\n" unless -f $ref_hit;

my ($chrom, $total_test, $total_ref) = check_file();
my ($window_test, $window_ref) = window_size();
my $cnvout = outfile_name();

# reading hit files
# I hope there are no funny characters in your chromosome names, such as *[]+?\
my $pat = qr/^$chrom\t(\d+)$/;
my $which = $total_test / ($total_test + $total_ref);

open(TEST,$test_hit) or die;
open(REF, $ref_hit) or die;
open(OUT, ">$cnvout.count") or die;
print OUT "chromosome\tstart\tend\ttest\tref\n";
my ($t, $r, $test, $ref) = (0,0,0,0);
my @win;
# might miss a few reads at the very end
for(1..($total_test + $total_ref))
{
	# randomize the order of test and ref files,
	# to prevent clustering of reads at a single spot:
	# eg: 10000 test and ref reads at the same position x
	if($r || (!($t) && rand() < $which))
	{
		my $temp = <TEST>;
		if($temp =~ m/$pat/)
		{
			$t = $1;
		}
	}else
	{
		my $temp = <REF>;
		if($temp =~ m/$pat/)
		{
			$r = $1;
		}
	}

	# if both $t and $r are initialized, 
	if($t && $r)
	{
		if($t < $r)
		{
			$test++;
			push @win, $t;
			$t = 0;
		}elsif($t > $r)
		{
			$ref++;
			push @win, $r;
			$r = 0;
		}else
		{
			if(rand() < $which)
			{
				$test++;
				push @win, $t;
				$t = 0;
			}else
			{
				$ref++;
				push @win, $r;
				$r = 0;
			}
		}
		if($test >= $window_test || $ref >= $window_ref)
		{
			my $start = min(@win);
			my $stop = max(@win);
			print OUT "$chrom\t$start\t$stop\t$test\t$ref\n";
#			print  "$chrom\t$start\t$stop\t$test\t$ref\n";
			$test = 0;
			$ref = 0;
			@win = ();
		}
	}
}
close OUT;

open(PATH, ">$cnvout.log") or die;
print PATH <<"PATH";
test hit file:	$test_hit
ref hit file:	$ref_hit
hit count file:	$cnvout.count
raw CNV output:	$cnvout.cnv.raw
consequetive window output:	$cnvout.cw
CBS output:	$cnvout.cbs
rough HMM output:	$cnvout.hmm.rough
PATH
close(PATH);

my $outcmd  = qq|cnv.print(cnv, "hmm", file="$cnvout.hmm.rough")\n| if $hmm eq 'TRUE';
$outcmd .= qq|cnv.print(cnv, "cw", file="$cnvout.cw")\n| if $cw eq 'TRUE';
$outcmd .= qq|cnv.print(cnv, "cbs", file="$cnvout.cbs")\n| if $cbs eq 'TRUE';

print "check out the log for a list of output:\n\t$cnvout.log\n";
open(R, ">$cnvout.r") or die;
print R <<"R";
	source("$seghmm_path/lib.segHMM.r");
	cnv<-cnv.cal("$cnvout.count", log2=$log2, min=$min_window, cw=$cw, hmm=$hmm, cbs=$cbs);
	write.table(cnv,"$cnvout.cnv.raw",row.names=FALSE,sep="\\t");
	$outcmd
R
system("$Rexe --vanilla --slave < $cnvout.r");

sub outfile_name
{
	my $temp = $test_hit;
	$temp =~ s/.+\///;
	my $out = $temp;
	$temp = $ref_hit;
	$temp =~ s/.+\///;
	$out .= "-vs-$temp";
	$out .= ".log2-$log2.pvalue-$pvalue";
	$out .= ".minw-$min_window";
	$out .= ".cw" if $cw eq 'TRUE';
	$out .= ".cbs" if $cbs eq 'TRUE';
	$out .= ".hmm" if $hmm eq 'TRUE';

	return($out);
}

# determining genome size
# and briefly checking file format
# only checking head and tail
# sorry, no full file check due to speed reason
sub check_file
{
	my $total_test = `wc -l $test_hit`;
	my $total_ref = `wc -l $ref_hit`;

	my @temp;
	my @temp_chr;
	for my $f($test_hit, $ref_hit)
	{
		my $temp = `head -n 100 $f`;
		if($temp =~ m/^(\S+)\t(\d+)$/m)
		{
			push @temp, $2;
			push @temp_chr, $1;
		}else
		{
			die "please check your input file ($f): no valid entry in first 100 lines?\n";
		}
		$temp = `tail -n 100 $f`;
		$temp =~ s/\s+$//;
		if($temp =~ m/(\S+)\t(\d+)$/)
		{
			push @temp, $2;
			push @temp_chr, $1;
		}else
		{
			die "please check your input file ($f): no valid entry in the last 100 lines?\n";
		}
	}
	@temp_chr = sort @temp_chr;
	die "only one chromosome is allowed, please check your input files\n" if $temp_chr[0] ne $temp_chr[3];
	die "best-hit file must be sorted\n" unless $temp[0] < $temp[1] && $temp[2] < $temp[3];

	return(
		$temp_chr[0],                    # chromosome name
		$total_test,                     # total reads in test genome
		$total_ref                      # total reads in ref genome
	)
}
	
sub window_size
{
	my $bt = `echo 'options(digits=16); qnorm(1-0.5*$pvalue)' | $Rexe --vanilla --slave`;
	die "\n Error:\ncan not find program $Rexe" unless $bt;
	my $st = `echo 'options(digits=16); qnorm(0.5*$pvalue)' | $Rexe --vanilla --slave`;
	$bt = $1 if $bt =~ m/^\[1\]\s+([\d.e\+\-]+)/m;
	$st = $1 if $st =~ m/^\[1\]\s+([\d.e\+\-]+)/m;
	
	$log2 = abs($log2);
	my $brp = 2**$log2;
	my $srp = 1/(2**$log2);

	my $bw = ($total_test * $brp**2 + $total_ref) * $bt**2 / ((1-$brp)**2);
	my $bwt = $bw / $total_ref;
	my $bwr = $bw / $total_test;
	my $sw = ($total_test * $srp**2 + $total_ref) * $st**2 / ((1-$srp)**2);
	my $swt = $sw / $total_ref;
	my $swr = $sw / $total_test;
	my $window_test = max($bwt, $swt);
	my $window_ref  = max($bwr, $swr);
	printf "... minimum reads for detecting log2>= $log2 should be %.0f(test) and %.0f(ref)\n", $bwt, $bwr;
	printf "... minimum reads for detecting log2<=-$log2 should be %.0f(test) and %.0f(ref)\n", $swt, $swr;
	$window_test = sprintf("%.0f", $window_test * $bigger);
	$window_ref  = sprintf("%.0f", $window_ref  * $bigger);
	print "... each window should contain, after X $bigger:\n  \tminimum $window_test test reads or minimum $window_ref ref reads\n";
	return($window_test, $window_ref);
}

sub stage2
{
	my ($testf, $reff, $infile, $outfile, $unreliablefile, $seghmm_path) = @_;
	my $temp = $outfile;

	my $total_test = `wc -l $testf`;
	$total_test = $1 if $total_test =~ m/^(\d+)/;
	my $total_ref = `wc -l $reff`;
	$total_ref  = $1 if $total_ref  =~ m/^(\d+)/;
	die "please check your $testf and $reff\n" unless $total_test>0 && $total_ref>0;

	open(INF, $infile) or die;
	my %seg;
	my %temp;
	my $min = 1e10;
	while(<INF>)
	{
#hmm	chromosome	start	end	hmm.size	hmm.log2
#hmm_1	chr1	10180656	11913413	1717915	-1
		if(m/^(?:\w*?)(\d+)\t(?:chr)?(\S+)\t(\d+)\t(\d+)\t(\d+)\t(\S+)$/)
		{
			my ($id, $chr, $start, $end, $size, $log2) = ($1, $2, $3, $4, $5, $6);
			$seg{$start} = {end=>$end, log2=>$log2, id=>$id};
			$temp{$chr} = 1;
			$min = $size if $min > $size;
		}
	}
	my @temp = keys %temp;
	die "one chromosome a time please\n" if $#temp > 0;
	die "please check $infile\n" unless @temp;
	my $chr  = shift @temp;

	mkdir "$temp.cand" or die "can not create dir $temp.cand: $!\n";
	open(RT, ">$temp.cand/$temp.total") or die;
	print RT "chromosome\ttest\tref\n";
	print RT "$chr\t$total_test\t$total_ref\n";
	close RT;

# maximum extension from a potential boundary
	$min *= 2;
#print "$min\n";
	my @starts = sort {$a <=> $b} keys %seg;
	unshift @starts, 0;
	push @starts, 3e10;
	$seg{0} = {end=>0, log2=>0, id=>0};
	$seg{3e10} = {end=>3e10, log2=>0, id=>0};

	open(R, ">$temp.cand/$temp.r") or die;
	print R "source('$seghmm_path/lib.segHMM.r');\n";
	my $header = 'TRUE';
	my @cand;
	my $cheader = "id\tside\tchromosome\tfrom\tto\tmid\tlog2.a\tlog2.b\n";
	my $rheader = "chromosome\tlocation\n";
	my @region;
	for my $i(1..($#starts-1))
	{
		my $left = $starts[$i];
		my $right = $seg{$starts[$i]}->{end};
		my $log2 = $seg{$starts[$i]}->{log2};
		my $id = $seg{$starts[$i]}->{id};

		my $left_to = min(($left + $right)/2, $left + $min);
		my $right_from = max(($left + $right)/2, $right - $min);

		my $pre_start = $starts[$i-1];
		my $pre_end = $seg{$pre_start}->{end};
		my $left_from = $left - $min;
		my $left_log2 = $seg{$pre_start}->{log2};
		# if gap between two segs
		if($left - $pre_end > $min/10 || $log2 == $left_log2)
		{
			$left_from = max($left_from, ($left + $pre_end)/2);
			$left_log2 = 0;
		}else
		{
			$left_from = max($left_from, ($left + $pre_start)/2);
		}

		my $next_start = $starts[$i+1];
		my $next_end = $seg{$next_start}->{end};
		my $right_to = $right + $min;
		my $right_log2 = $seg{$next_start}->{log2};
		# if gap between two segs
		if($next_start - $right > $min/10 || $log2 == $right_log2)
		{
			$right_to = min($right_to, ($right + $next_start)/2);
			$right_log2 = 0;
		}else
		{
			$right_to = min($right_to, ($right + $next_end)/2);
		}

		$left_from = int($left_from);
		$left_to = int($left_to);
		$right_from = int($right_from);
		$right_to = int($right_to);

		push @region, [$id, $left_from, $left_to];
		push @region, [$id, $right_from, $right_to];

		open(CC, ">$temp.cand/$id.cand") or die;
		print CC $cheader;
		print CC "$id\tleft\t$chr\t$left_from\t$left_to\t$left\t$left_log2\t$log2\n";
		print CC "$id\tright\t$chr\t$right_from\t$right_to\t$right\t$log2\t$right_log2\n";
		close(CC);

		open(RT, ">$temp.cand/$id.test") or die;
		print RT $rheader;
		close RT;
		open(RR, ">$temp.cand/$id.ref") or die;
		print RR $rheader;
		close RR;

		print R "hmm.bound('$temp.cand/$id.cand', '$temp.cand/$id.test', '$temp.cand/$id.ref', '$temp.cand/$temp.total', '$outfile', '$unreliablefile', header = $header)\n";
		$header = 'FALSE';
	}
	close(R);

# extracting reads

	my %rfile = ('test' => $testf, 'ref' => $reff); 
	my $pat = qr/^$chr\t(\d+)$/;
	for my $g(qw(test ref))
	{

		my @starts = @region;
		my ($id, $start, $end) = @{shift @starts};
		open(T, $rfile{$g}) or die;
		R: while(<T>)
		{
			if(m/$pat/)
			{
				my $pos = $1;
				while($pos > $end)
				{
					my $temp = shift @starts;
					if($temp)
					{
						($id, $start, $end) = @$temp;
					}else
					{
						last R;
					}
				}

				if($pos >= $start)
				{
					open(R, ">>$temp.cand/$id.$g") or die;
					print R $_;
				}

			}
		}
	}

	system("R --vanilla --slave < $temp.cand/$temp.r");
}
