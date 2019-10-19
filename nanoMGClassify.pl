#!/usr/bin/env perl

use Cwd qw(realpath); #https://perldoc.perl.org/Cwd.html get pathname of current working directory; resolve full path with filesystem check
use Fcntl 'O_RDONLY'; #https://perldoc.perl.org/Fcntl.html load the C Fcntl.h defines; needed for Tie::File
use File::Basename; #https://perldoc.perl.org/File/Basename.html Parse file paths into directory, filename and suffix.
use File::Spec; #https://perldoc.perl.org/File/Spec.html portably perform operations on file names and paths; more extensive than File::Basename
use Getopt::Long; #https://perldoc.perl.org/Getopt/Long.html Extended processing of command line options
use IO::Select; #https://perldoc.perl.org/IO/Select.html OO interface to the select system call; allows the user to see what IO handles, are ready for reading, writing or have an exception pending
use IPC::Open3; #https://perldoc.perl.org/IPC/Open3.html open a process for reading, writing, and error handling using open3()
use POSIX qw(strftime); #https://perldoc.perl.org/POSIX.html Perl interface to IEEE Std 1003.1
use strict; #https://perldoc.perl.org/strict.html Perl pragma to restrict unsafe constructs
use Symbol; #https://perldoc.perl.org/Symbol.html manipulate Perl symbols and their names; needed for working with open3
use Tie::File; #https://perldoc.perl.org/Tie/File.html access the lines of a disk file via a Perl array
use Time::HiRes qw(time); #https://perldoc.perl.org/Time/HiRes.html High resolution alarm, sleep, gettimeofday, interval timers
use Time::Seconds; #https://perldoc.perl.org/Time/Seconds.html a simple API to convert seconds to other date values
use v5.10;
use warnings; #https://perldoc.perl.org/warnings.html Perl pragma to control optional warnings


use diagnostics; #https://perldoc.perl.org/diagnostics.html produce verbose warning diagnostics; at least during development

#########################
######### Usage #########
#########################

sub Info
	{
	say "
################################################################################
## nanoMGClassify.pl
## Accurate metagenomic classification of long-read Nanopore DNA sequences using
## ensemble methods.
##
## Input: FILL INPUT DETAILS
##
## Output: FILL OUTPUT DETAILS
##
## Dependencies: MENTION DEPENDENCIES POSSIBLY WITH VERSION
##
## Author: piyuranjan\@gmail.com
## Source: https://github.com/piyuranjan/NanoMGClassify/blob/master/nanoMGClassify.pl
################################################################################
	";
	}
sub Usage
	{
	say "
Usage:\n$0 [options] <arguments> <*.fastq[.gz]|fastqFile1..n> >[outputFile]

At least one of the following is required
 *.fastq[.gz]   Feel free to use a pattern for .fastq or .fastq.gz files.
 fastqFile1..n  Specify files as positional arguments.

Arguments: required
 -b|blIdx    <string> Path to Blast database with index prefix
               (minus trailing .nXX).
 -c|cfIdx    <string> Path to Centrifuge database with index prefix
               (minus trailing .X.cf).

Options:
 --bestHit   Force reporting of best hit from both Centrifuge and Blast.
               Using this option disables LCA by design. Default: Off
 -f|force    Force overwrite outFile if it exists. Depends on -o|outFile.
 --noTaxa    Do not use taxonomy information for faster classification.
               Default: Off
 -o|outFile  <string> Send output to a file. STDOUT otherwise.
 -t|threads  <integer> Number of threads/processors to run. Default: 1
 --track     <integer> Track progress after every N sequences. Default: 50
 -h|help     Show more help and exit 0.
 -v|verbose  Use multiple times to increase verbosity.
 -q|quiet    Silent execution except for significant/fatal errors (ERR).
               No warnings (WARN), statistics (STAT) and information (INFO)
               when -q engaged. Make sure you know what you're doing.
 --debug     Set max verbosity for debugging.
	";
	}
sub Examples
	{
	say "
Examples:

### Newbie: HOW TO USE FOR A NEW PERSON.
$0 *.fastq.gz

### Expert: AN EXPERT USE CASE.
$0 -v -v -o summary.txt Path/To/Fastq/*.fastq.gz *.fastq

### Wizard: BEYOND EXPERT OR FULL SCALE USE, FOR EXAMPLE: Pair it with GNU Parallel.
parallel -j 6 $0 -v {} ::: Path/to/Fastq/*.fastq.gz | sort | uniq >summary.txt
	";
	}
sub StatusInfo
	{
	say "
Status message categories:
 - DBG:  Debug statements only when debugging mode is on.
 - ERR:  Fatal errors!
 - INFO: Information about the program only when verbose mode is on.
 - STAT: Statistics about time or other parameters.
 - STEP: Mark for important steps in the execution.
 - WARN: Non-fatal warnings!
 

Exit codes:
 0  Successful
 1  Unsuccessful or need arguments
 2  Problem with a dependency
 3  Problem with source FastQ file
 4  Problem with running Centrifuge
 5  Problem with running Blastn
 x  Maybe: Problem creating local copy of the fastq files
	";
	}


#########################
######### Main ##########
#########################

my $time0=time();

## Option reading and configuration ##

### Pre-process script specific arguments and defaults
my($bestHit,$blIdx,$cfIdx,$noTaxa,$threads,$trackAtSeq,$workDir);
$threads=1; #default thread/proc count is 1
$trackAtSeq=50; #track progress at every 50 fastq records
$workDir=realpath()."/NanoMGClassify"; #default work directory is PWD/NanoMGClassify (-yyyymmddhhmmss appended later) with full path resolved

### Pre-process generic arguments and defaults
my($debug,$quiet,$help,$verbose); #generic arguments
$verbose=1; #set default verbose to 1

### Read-in arguments
unless(GetOptions(
				#script specific arguments
				'bestHit' => \$bestHit,
				'b|blIdx=s' => \$blIdx,
				'c|cfIdx=s' => \$cfIdx,
				'noTaxa' => \$noTaxa,
				't|threads=i' => \$threads,
				'track=i' => \$trackAtSeq,
				'w|workDir=s' => \$workDir,
				#generic arguments
				'debug' => \$debug,
				'q|quiet' => \$quiet,
				'h|help' => \$help,
				'v|verbose+' => \$verbose,
				))
	{Usage;exit 1;} #quit with error code
if($help) #quit with full help
	{Info;Usage;Examples;StatusInfo;exit 0;}

### Set necessary arguments
unless(defined $ARGV[0])
	{warn TimeStamp($verbose)."ERR: Need fastq files or file pattern.\n";Usage;exit 1;}
unless(defined $blIdx)
	{warn TimeStamp($verbose)."ERR: Need Blast index path.\n";Usage;exit 1;}
unless(defined $cfIdx)
	{warn TimeStamp($verbose)."ERR: Need Centrifuge index path.\n";Usage;exit 1;}

### Post-process option behaviors, arguments and defaults
chomp($workDir);
$workDir=~s/\/$//; #remove trailing slash from user provided path
$workDir=$workDir."-".strftime("%Y%m%d%H%M%S",localtime(time)); #default work directory appended with local timestamp
$verbose=0 if($quiet);
if($debug) #override verbosity to highest value for debugging code
	{$verbose=100;warn TimeStamp($verbose)."INFO: Debug mode enabled. Setting verbosity to max.\n";}
if(defined $bestHit)
	{warn TimeStamp($verbose)."WARN: --bestHit is enabled! This is not recommended. Please make sure you know what it does by using -h|help.\n" if($verbose);}
if(defined $noTaxa)
	{warn TimeStamp($verbose)."WARN: --noTaxa is enabled! This is not recommended. Please make sure you know what it does by using -h|help.\n" if($verbose);}

### Time stats after option reading
my $timeReadOptions=time();
warn TimeStamp($verbose)."DBG: STAT: Time to read options: ".TimeDiffPretty($time0,$timeReadOptions)."\n" if($debug);


## Check package dependencies ##

### Not needed right now ###
### Check for GNU Parallel
# my $parallelVersion="20161222"; #min version needed for this application
# unless(CheckParallel($verbose,$parallelVersion))
	# {warn TimeStamp($verbose)."ERR: Need GNU Parallel. Exiting\n";exit 2;}

### Check for Centrifuge
my $centrifugeVersion="1.0.4"; #min version needed for this application
unless(CheckCentrifuge($verbose,$centrifugeVersion))
	{warn TimeStamp($verbose)."ERR: Need Centrifuge. Exiting\n";exit 2;}

### Check for Blast
my $blastVersion="2.8.1"; #min version needed for this application
if(CheckBlast($verbose,$blastVersion)<1)
	{warn TimeStamp($verbose)."ERR: Need NCBI Blast higher than $blastVersion. Exiting\n";exit 2;}

### Time stats after checking dependencies
my $timeCheckDependency=time();
warn TimeStamp($verbose)."DBG: STAT: Time to check dependencies: ".TimeDiffPretty($timeReadOptions,$timeCheckDependency)."\n" if($debug);


## Scan and prepare sequencing files for classification ##

### Establish work directory
warn TimeStamp($verbose)."INFO: Establishing work directory.\n" if($verbose>2);
my $workDirExists=IsValidDirectory($workDir,$verbose);
unless($workDirExists)
	{warn TimeStamp($verbose)."INFO: Creating directory $workDir\n" if($verbose>1);mkdir $workDir or die $!;}
else
	{
	warn TimeStamp($verbose)."ERR: Work directory already exists, $workDir\n";
	warn TimeStamp($verbose)."DBG: Directory code: \'$workDirExists\'\n" if($debug);
	exit 0;
	}
$workDir=realpath($workDir); #resolve full path for work dir if prefix provided by user
warn TimeStamp($verbose)."INFO: Work directory established: $workDir\n" if($verbose>2);

### Scan all fastq[.gz]
warn TimeStamp($verbose)."INFO: Scanning all arguments for fastq[.gz] files.\n" if($verbose>2);
my @fqgzFiles;
foreach my $arg(@ARGV) #read all user arguments and scan all matching files
	{
	warn TimeStamp($verbose)."DBG: User argument returned by shell is \'$arg\'\n" if($debug); #reports arguments passed to the script by the user and if unix shell expanded anything; will help pattern debugging
	$arg=File::Spec->catfile($arg,'*') if(-d $arg); #if arg is a directory, add /* so it matches all files when globbed
	my @files=glob("$arg"); #scan everything that matches the pattern
	### Make sure ONLY \.f(ast)?q(\.gz)? extension files are pushed (no folder or other extensions)
	foreach(@files)
		{if(/\.f(ast)?q(\.gz)?$/){chomp;push(@fqgzFiles,$_);}}
	}
warn TimeStamp($verbose)."STAT: Found ".scalar(@fqgzFiles)." file(s) to classify.\n" if($verbose);
if(scalar(@fqgzFiles)<1){warn TimeStamp($verbose)."ERR: No files to classify. Exiting.\n";exit 3;} #exit if no files found
if($verbose>1) #list all files found
	{warn TimeStamp($verbose)."INFO: List of files:\n";warn "$_\n" foreach @fqgzFiles;}

### Extract any .gz if fastq are compressed
### If they are already uncompressed, copy them over.
### Right now, this is a single threaded operation and can be transformed on multiple threads later.
warn TimeStamp($verbose)."STEP: Starting to build local copies of uncompressed fastq.\n" if($verbose);
my $fqDir=$workDir."/00.FastqFiles"; #this direcotry will hold all extracted fastq
mkdir $fqDir or die $!;
warn TimeStamp($verbose)."INFO: Succesfully created directory for local fastq: $fqDir\n" if($verbose>2);
my @fqFiles;
foreach my $file(@fqgzFiles)
	{
	my(undef,undef,$newFile)=File::Spec->splitpath($file); #extract .fq.gz fileName from path
	$newFile=~s/\.gz$//;
	$newFile=File::Spec->catfile($fqDir,$newFile); #prepare .fq fileName with full path
	my $cmd;
	if($file=~/\.gz$/) #extract if extension is .gz
		{$cmd="gzip -cd $file >$newFile";}
	else #copy otherwise
		{$cmd="cp $file $newFile";}
	warn TimeStamp($verbose)."INFO: Executing command:\n`$cmd`\n" if($verbose>2);
	my $cmdOut=`$cmd`; #not inserting 2>&1 because one command is using 1>file
	warn TimeStamp($verbose)."DBG: Command return status: \'$?\'\n" if($debug);
	if($?!=0){warn TimeStamp($verbose)."ERR: Failed to create local fastq.\n";exit 1;}
	push(@fqFiles,$newFile);
	}
warn TimeStamp($verbose)."INFO: Succesfully created local fastq for ".scalar(@fqFiles)." file(s).\n" if($verbose>2);
if($verbose>3) #list all local fastq created
	{warn TimeStamp($verbose)."INFO: List of all files created:\n"; warn "$_\n" foreach(@fqFiles);}

### Time stats after creating local fastq
my $timeLocalFastq=time();
warn TimeStamp($verbose)."DBG: STAT: Time to create local fastq: ".TimeDiffPretty($timeCheckDependency,$timeLocalFastq)."\n" if($debug);


## Execute centrifuge and blast on each fastq ##

######## add package name variables and modify statements accordingly.

### Prepare work directories for centrifuge and blast
warn TimeStamp($verbose)."STEP: Starting to run Centrifuge and BLAST on each fastq.\n" if($verbose);
# my $cfDir=$workDir."/01.CentrifugeClassification"; #this directory will hold all centrifuge outputs
my $l1Dir=$workDir."/01.CentrifugeClassification"; #this directory will hold all centrifuge outputs
mkdir $l1Dir or die $!;
warn TimeStamp($verbose)."INFO: Succesfully created directory for Centrifuge output: $l1Dir \n" if($verbose>2);
# my $blDir=$workDir."/02.BlastClassification"; #this directory will hold all blast outputs
my $l2Dir=$workDir."/02.BlastClassification"; #this directory will hold all blast outputs
mkdir $l2Dir or die $!;
warn TimeStamp($verbose)."INFO: Succesfully created directory for Blast output: $l2Dir \n" if($verbose>2);

### Cycle on each fastq for classification
my $fqFileCount=0;
foreach my $fqFile(@fqFiles)
	{
	my $timeFastq0=time(); #start time working for current fastq
	$fqFileCount++;
	warn TimeStamp($verbose)."INFO: Classifying sequences in fastq $fqFileCount: ".basename($fqFile)."\n" if($verbose>1);
	
	### Run Level-1 classification with the fqFile
	my $l1Idx=$cfIdx; #change this statement before this level to make l1Idx standard terminology
	my %seqTaxaL1;
	my $l1FqCount=ExecCentrifuge($fqFileCount,$fqFile,\%seqTaxaL1,$l1Dir,$threads,$l1Idx,$bestHit,$verbose); #use Centrifuge for L1 classification
	
	### Time stats after forking and parsing L1 output
	my $timeExecL1=time();
	warn TimeStamp($verbose)."INFO: STAT: Time to process $l1FqCount sequences at level-1 with Centrifuge: ".TimeDiffPretty($timeFastq0,$timeExecL1)."\n" if($verbose>2);
	
	### Run Level-2 classification with the fqFile
	my $l2Idx=$blIdx; #change this statement before this level to make l2Idx standard terminology
	my $l2FqCount=ExecBlast($fqFileCount,$fqFile,\%seqTaxaL1,$l2Dir,$threads,$l2Idx,$noTaxa,$bestHit,$trackAtSeq,$verbose); #use Blast for L2 classification
	
	### Time stats after forking and parsing L2 output
	my $timeExecL2=time();
	warn TimeStamp($verbose)."INFO: STAT: Time to process $l2FqCount sequences at level-2 with Blast: ".TimeDiffPretty($timeExecL1,$timeExecL2)."\n" if($verbose>2);
	
	### Prepare Krona report from Blast results
	
	
	### Time stats after all steps for current fastq
	my $timeFastqAll=time();
	# warn TimeStamp($verbose)."STAT: Finished sequence classification in ".TimeDiffPretty($timeFastq0,$timeFastqAll)." from fastq $fqFileCount.\n" if($verbose);
	warn TimeStamp($verbose)."STAT: Time to finish sequence classification from fastq $fqFileCount: ".TimeDiffPretty($timeFastq0,$timeFastqAll)."\n" if($verbose);
	}


## Print stats of operation ##

### Time stats after everything
my $timeEnd=time();
warn TimeStamp($verbose)."STAT: Time to finish classification for ".scalar(@fqFiles)." fastq file(s): ".TimeDiffPretty($time0,$timeEnd)."\n" if($verbose);


#########################
###### Subroutines ######
#########################

sub ExampleFuntion #BRIEF DESCRIPTION OF THE FUNCTION
	{
	### Order of Input: DESCRIBE THE INPUTS AND THEIR ORDER. For example: fileName-req (required), inBuffer-opt (optional)
	### Order of Output: DESCRIBE THE OUTPUTS AND THEIR ORDER.
	}

sub CheckBlast #Checks if the NCBI Blast package is available on PATH, returns good/bad status
	{
	### Order of Input: verboseLevel-opt; versionOfBlastToCheck-opt
	### Order of Output: blastStatus (0:notFound, 1:good, -1:lower)
	my ($verbose,$needVersion,@rest)=@_;
	$verbose//=1;
	$needVersion//="2.8.1";
	warn TimeStamp($verbose)."INFO: Executing function CheckBlast\n" if($verbose>2);
	
	## Check blast and extract version
	my $blastCmd="blastn -version";
	warn TimeStamp($verbose)."INFO: Running NCBI Blast+ command:\n`$blastCmd`\n" if($verbose>3);
	my @blastOut=`$blastCmd`;
	if($verbose>3){warn $_ foreach @blastOut;}
	# if($verbose>3)
		# {warn TimeStamp($verbose)."INFO: Finished Blast command:\n";warn $_ foreach @blastOut;}
	### If blast is not found in PATH, suggest easiest installation and return with 0
	if(!@blastOut && $verbose)
		{
		warn TimeStamp($verbose)."WARN: NCBI Blast+ not found in PATH\n";
		warn TimeStamp($verbose)."INFO: Please install NCBI Blast+ via `conda install blast` or `sudo apt install ncbi-blast+` \n";
		warn TimeStamp($verbose)."INFO: Or please check out https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download \n";
		return (0);
		}
	my $version=(split(/\h+/,$blastOut[0]))[-1]; chomp($version);
	warn TimeStamp($verbose)."INFO: Found NCBI Blast version: $version\n" if($verbose>2);
	
	## Compare blast version
	if($version ge $needVersion)
		{
		warn TimeStamp($verbose)."INFO: NCBI Blast $version found is same/higher than $needVersion needed. Good to go!\n" if($verbose>2);
		return (1);
		}
	else
		{
		warn TimeStamp($verbose)."WARN: NCBI Blast $version found is lower than $needVersion needed.\n" if($verbose);
		warn TimeStamp($verbose)."WARN: Program may not work as intended.\n" if($verbose);
		return (-1);
		}
	}

sub CheckCentrifuge #Checks if the Centrifuge package is available on PATH, returns good/bad status
	{
	### Order of Input: verboseLevel-opt; versionOfCentrifugeToCheck-opt
	### Order of Output: centrifugeStatus (0:notFound, 1:good, -1:lower)
	my ($verbose,$needVersion,@rest)=@_;
	$verbose//=1;
	$needVersion//="1.0.4";
	warn TimeStamp($verbose)."INFO: Executing function CheckCentrifuge\n" if($verbose>2);
	
	## Check centrifuge and extract version
	my $centrifugeCmd="centrifuge --version";
	warn TimeStamp($verbose)."INFO: Running Centrifuge command:\n`$centrifugeCmd`\n" if($verbose>3);
	my @centrifugeOut=`$centrifugeCmd`;
	if($verbose>3){warn $_ foreach @centrifugeOut;}
	# if($verbose>3)
		# {warn TimeStamp($verbose)."INFO: Finished Centrifuge command:\n";warn $_ foreach @centrifugeOut;}
	### If centrifuge is not found in PATH, suggest easiest installation and return with 0
	if(!@centrifugeOut && $verbose)
		{
		warn TimeStamp($verbose)."WARN: Centrifuge package not found in PATH\n";
		warn TimeStamp($verbose)."INFO: Please install Centrifuge via `conda install centrifuge` \n";
		warn TimeStamp($verbose)."INFO: Or please check out https://ccb.jhu.edu/software/centrifuge/ \n";
		return (0);
		}
	my $version=(split(/\h+/,$centrifugeOut[0]))[-1]; chomp($version);
	warn TimeStamp($verbose)."INFO: Found Centrifuge version: $version\n" if($verbose>2);
	
	## Compare centrifuge version
	if($version ge $needVersion)
		{
		warn TimeStamp($verbose)."INFO: Centrifuge $version found is same/higher than $needVersion needed. Good to go!\n" if($verbose>2);
		return (1);
		}
	else
		{
		warn TimeStamp($verbose)."WARN: Centrifuge $version found is lower than $needVersion needed.\n" if($verbose);
		warn TimeStamp($verbose)."WARN: Program may not work as intended.\n" if($verbose);
		return (-1);
		}
	}

sub CheckParallel #Checks if GNU parallel is available on PATH, returns good/bad status
	{
	### Order of Input: verboseLevel-opt; versionOfParallelToCheck-opt
	### Order of Output: parallelStatus (0:notFound, 1:good, -1:lower)
	my ($verbose,$needVersion,@rest)=@_;
	$verbose//=1;
	$needVersion//="20161222";
	warn TimeStamp($verbose)."INFO: Executing function CheckParallel\n" if($verbose>2);
	
	## Check parallel and extract version
	my $parallelCmd="parallel -V";
	warn TimeStamp($verbose)."INFO: Running GNU Parallel command:\n`$parallelCmd`\n" if($verbose>3);
	my @parallelOut=`$parallelCmd`;
	if($verbose>3){warn $_ foreach @parallelOut;}
	# if($verbose>3)
		# {warn TimeStamp($verbose)."INFO: Finished Parallel command:\n";warn $_ foreach @parallelOut;}
	### If parallel is not found in PATH, suggest easiest installation and return with 0
	if(!@parallelOut && $verbose)
		{
		warn TimeStamp($verbose)."WARN: GNU Parallel not found in PATH\n";
		warn TimeStamp($verbose)."INFO: Please install GNU Parallel via `conda install parallel` or `sudo apt install parallel`\n";
		warn TimeStamp($verbose)."INFO: Or please check out https://www.gnu.org/software/parallel/ \n";
		return (0);
		}
	my $version=(split(/\h+/,$parallelOut[0]))[-1]; chomp($version);
	warn TimeStamp($verbose)."INFO: Found GNU Parallel version: $version\n" if($verbose>2);
	
	## Compare parallel version
	if($version ge $needVersion)
		{
		warn TimeStamp($verbose)."INFO: GNU Parallel $version found is same/higher than $needVersion needed. Good to go!\n" if($verbose>2);
		return (1);
		}
	else
		{
		warn TimeStamp($verbose)."WARN: GNU Parallel $version found is lower than $needVersion needed.\n" if($verbose);
		warn TimeStamp($verbose)."WARN: Program may not work as intended.\n" if($verbose);
		return (-1);
		}
	}

sub ExecBlast #Executes blastn (bl) with parameters and returns parsed sequence taxa
	{
	### Order of Input: fqFileCount-req, fqFile-req, \%seqIDTaxID-req (a ref to seqID:taxa->taxID1,taxID2,... and seqID->taxID->refID1,... mappings), blWorkDir-req, threads-opt, blIndexPath-req, noTaxaUseAlgoFlag-opt, bestHitAlgoFlag-opt, trackFqProgressAfterNSeq-opt, verboseLevel-opt
	### Order of Output: fqRecordsProcessedCount-opt
	my $timeSub0=time();
	my $packageName="Blast";
	my ($fqFileCount,$fqFile,$seqTaxaRef,$cmdWorkDir,$threads,$idxPath,$noTaxa,$bestHit,$trackAtSeq,$verbose,@rest)=@_;
	die "ERR: Incomplete parameters for function to work." unless($fqFileCount && $fqFile && $seqTaxaRef && $cmdWorkDir && $idxPath);
	# die "ERR: Incomplete parameters for function to work." if(!$fqFile||!$cmdWorkDir);
	$threads//=1;
	$trackAtSeq//=50; #default progress tracking is at every 50 fastq
	$verbose//=1;
	my $debug; $debug=1 if($verbose>=100);
	
	### Prepare blast output file names
	my $outName;
	(undef,undef,$outName)=File::Spec->splitpath($fqFile); #extract .fq fileName from path
	$outName=~s/\.f(ast)?q.*$//; #delete extension
	my $outBase=File::Spec->catfile($cmdWorkDir,$outName); #fileBase with full path
	my $cmdOutFile=$outBase."_blastMatches.txt"; #blast top hits fileName with full path
	
	### Prepare fasta file buffer from fastq and feed to blast with taxonomy
	warn TimeStamp($verbose)."INFO: Starting to buffer fastq to fasta and feeding $packageName.\n" if($verbose>2);
	tie my @fqFileReader, 'Tie::File', $fqFile, autochomp => 0, memory => 0, mode => O_RDONLY or die $!; #open fastq as an array
	open(my $CMDEXP,">$cmdOutFile") or die $!; #open file to write the output as is
	my @blCols=qw(qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid); #this will be the order of columns in blast output
	print $CMDEXP "#".join("\t#",@blCols)."\n"; #print column headers to blast files
	my $fqMaxLines=$#fqFileReader; #store max lines to avoid recaching
	my $fqCountTotal=($fqMaxLines+1)/4;
	my $fqCount=0; my $fqFailed=0;
	my $fqUsedTaxa=0;
	my $fqLine=0;
	my $taxa="taxa"; #key to address taxonomy IDs for seqID
	while($fqLine<=$fqMaxLines)
		{
		my $timeSeq0=time();
		
		### Prepare fasta buffer
		my ($seqID)=$fqFileReader[$fqLine]=~/^@(.*?)\s/; #isolate seqID
		unless($seqID) #check for fastq format problems
			{warn TimeStamp($verbose)."ERR: Format corruption. Files are probaly not fastq format.\n";exit 3;}
		my $faBuffer=">$seqID\n".$fqFileReader[$fqLine+1]; #read sequence and store as fasta buffer
		$fqLine+=4; #advance line count to next fastq record
		# warn TimeStamp($verbose)."DBG: Fasta buffer:\n---\n$faBuffer\n---\n" if($debug); #uncomment to see how sequences are being parsed
		
		### Construct package command based on taxonomy ID from L1 classification
		my $cmd="blastn -db $idxPath -evalue 0.01 -num_threads 1 -outfmt \"6 @blCols\"";
		if($$seqTaxaRef{$seqID}{$taxa} && !defined($noTaxa)) #append taxonomy data for blastn search if the taxa ID is non-zero and is not disabled
			{$cmd.=" -taxids $$seqTaxaRef{$seqID}{$taxa}";$fqUsedTaxa++;}
		$cmd.=" -max_target_seqs 1" if(defined $bestHit); #engage reporting best hit if enabled.
		
		### Execute package with forking and feed the fasta buffer
		my($CMDWT,$CMDRD,$CMDERR);
		$CMDERR=gensym(); #create a symbol for error filehandle (FH); open3 doesn't do this automatically
		warn TimeStamp($verbose)."DBG: Executing command: `$cmd`\n" if($debug); #uncomment to see command for each fasta record buffered
		my $cmdPid;
		eval {$cmdPid=open3($CMDWT,$CMDRD,$CMDERR,$cmd);}; #execute open3 in eval to capture errors
		die "ERR: Couldn't fork $packageName\n$@" if($@); #throw out the error
		print $CMDWT $faBuffer; #feed in the sequence
		close($CMDWT); #finish feeding the sequence
		
		### Start reading and storing the package output
		my ($cmdOutBuffer,$cmdErrBuffer)=('','');
		my $cmdSelect=new IO::Select; #create a select object to notify read on FHs
		$cmdSelect->add($CMDRD,$CMDERR); #add FHs
		warn TimeStamp($verbose)."DBG: Reading from $packageName:\n" if($debug);
		while(my @ready=$cmdSelect->can_read) #find all FHs that are ready for reading
			{
			foreach my $FH(@ready) #loop through each FH that is ready
				{
				my $buffer;
				my $size=sysread($FH,$buffer,4096); #read upto 4096 bytes in this iteration
				unless(defined $size) #quit if error reading
					{die "Error reading from $packageName child pid $cmdPid \n$!";}
				elsif($size==0) #finish reading and remove if the FH is empty
					{$cmdSelect->remove($FH);}
				else #otherwise read data
					{
					warn TimeStamp($verbose)."DBG: Read $size bytes from $FH \n" if($debug);
					if($FH==$CMDRD) #read output
						{
						$cmdOutBuffer.=$buffer;
						print $CMDEXP $buffer; #export to file as is
						}
					elsif($FH==$CMDERR) #read error
						{$cmdErrBuffer.=$buffer;}
					else {die "ERR: Extra filehandle detected $FH ";}
					}
				}
			}
		waitpid($cmdPid,0);
		
		### Time stats after package execution
		my $timeCmdExec=time();
		# warn TimeStamp($verbose)."DBG: STAT: Time to execute $packageName: ".TimeDiffPretty($timeSeq0,$timeCmdExec)."\n" if($debug); #uncomment to see package execution time for each seq
		
		### Error handling and output parsing for package command
		warn TimeStamp($verbose)."DBG: Command return status: \'$?\'\n" if($debug); #uncomment to see package command exit status for each execution
		if($?!=0) #if exit code non-zero, print error
			{
			warn TimeStamp($verbose)."WARN: Failed to run $packageName properly on $seqID from fastq $fqFileCount!\n" if($verbose);
			warn TimeStamp($verbose)."WARN: INFO: $packageName command used was: `$cmd`\n" if($verbose>1);
			warn TimeStamp($verbose)."WARN: $packageName said:\n---\n$cmdErrBuffer---\n" if($verbose);
			$fqFailed++; #count failed/unsuccessful executions.
			$fqUsedTaxa-- if($$seqTaxaRef{$seqID} && !defined($noTaxa)); #undo count of fastq using taxa if they did
			# exit 5; #uncomment if termination is desired on Blast fail
			}
		else #otherwise, parse output from command
			{}
		
		### Time stats after parsing package output
		my $timeCmdParse=time();
		# warn TimeStamp($verbose)."DBG: STAT: Time to parse $packageName output: ".TimeDiffPretty($timeCmdExec,$timeCmdParse)."\n" if($debug); #uncomment to see package parsing time for each seq
		
		$fqCount++; #count fastq record executions.
		
		### Status of fastq progress
		unless($fqCount % $trackAtSeq)
			{warn TimeStamp($verbose)."STAT: Sequences classified from fastq $fqFileCount: $fqCount of $fqCountTotal \n" if($verbose);}
		}
	close($CMDEXP);
	untie @fqFileReader; #close fastq filehandle
	
	### Stats for package performance after all executions for the current fastq
	if($verbose && $fqFailed) #failed stats
		{
		warn TimeStamp($verbose)."WARN: $fqFailed of $fqCountTotal sequences FAILED in $packageName for fastq $fqFileCount! : ".basename($fqFile)."\n";
		warn TimeStamp($verbose)."WARN: Please check the log above for messages from $packageName.\n";
		}
	warn TimeStamp($verbose)."STAT: $fqUsedTaxa of $fqCountTotal sequences in fastq $fqFileCount used taxonomy.\n" if($verbose); #taxonomy stats
	
	return($fqCount); #return number of fq seq records processed
	}

sub ExecCentrifuge #Executes centrifuge (cf) with parameters and returns parsed sequence taxa
	{
	### Order of Input: fqFileCount-req, fqFile-req, \%seqIDTaxID-req (a hash ref to store seqID:taxa->taxID1,taxID2,... and seqID->taxID->refID1,... mappings), cfWorkDir-req, threads-opt, cfIndexPath-req, bestHitAlgoFlag-opt, verboseLevel-opt
	### Order of Output: uniqueFqRecordsProcessedCount-opt
	### Order of Output: \%seqMatches-req (a reference to seqID:taxa->taxID1,taxID2,... and seqID->taxID->refID1,... mappings)
	my $timeSub0=time();
	my $packageName="Centrifuge";
	my ($fqFileCount,$fqFile,$seqTaxaRef,$cmdWorkDir,$threads,$idxPath,$bestHit,$verbose,@rest)=@_;
	die "ERR: Incomplete parameters for function to work." unless($fqFileCount && $fqFile && $seqTaxaRef && $cmdWorkDir && $idxPath);
	# die "ERR: Incomplete parameters for function to work." if(!$fqFile||!$cmdWorkDir);
	$threads//=1;
	$verbose//=1;
	my $debug; $debug=1 if($verbose>=100);
	# $bestHit//=0; #no default condition needed as this flag is inherited; becomes defined if set
	
	### Prepare centrifuge output file names
	my $outName;
	(undef,undef,$outName)=File::Spec->splitpath($fqFile); #extract .fq fileName from path
	$outName=~s/\.f(ast)?q.*$//; #delete extension
	my $outBase=File::Spec->catfile($cmdWorkDir,$outName); #fileBase with full path
	my $cfReport=$outBase."_centrifugeReport.txt"; #centrifuge report fileName with full path
	my $cmdOutFile=$outBase."_centrifugeClassify.txt"; #classification output fileName with full path
	
	### Construct centrifuge command
	my $cmd="centrifuge -p $threads -x $idxPath -U $fqFile --report-file $cfReport";
	$cmd.=" -k 1" if($bestHit); #engage reporting best hit if enabled.
	
	### Execute package with forking on the fastq
	my($CMDRD,$CMDERR);
	$CMDERR=gensym(); #create a symbol for error filehandle (FH); open3 doesn't do this automatically
	warn TimeStamp($verbose)."INFO: Executing command:\n`$cmd`\n" if($verbose>2);
	my $cmdPid;
	eval {$cmdPid=open3(undef,$CMDRD,$CMDERR,$cmd);}; #execute open3 in eval to capture errors
	die "ERR: Couldn't fork $packageName!\n$@" if($@); #throw out the error
	
	### Start reading and storing the package output
	open(my $CMDEXP,">$cmdOutFile") or die $!; #open file to write the output as is
	my ($cmdOutBuffer,$cmdErrBuffer)=('','');
	my $cmdSelect=new IO::Select; #create a select object to notify read on FHs
	$cmdSelect->add($CMDRD,$CMDERR); #add FHs
	warn TimeStamp($verbose)."DBG: Reading from $packageName:\n" if($debug);
	while(my @ready=$cmdSelect->can_read) #find all FHs that are ready for reading
		{
		foreach my $FH(@ready) #loop through each FH that is ready
			{
			my $buffer;
			my $size=sysread($FH,$buffer,4096); #read upto 4096 bytes in this iteration
			unless(defined $size) #quit if error reading
				{die "Error reading from $packageName child pid $cmdPid \n$!";}
			elsif($size==0) #finish reading and remove if the FH is empty
				{$cmdSelect->remove($FH);}
			else #otherwise read data
				{
				warn TimeStamp($verbose)."DBG: Read $size bytes from $FH \n" if($debug);
				if($FH==$CMDRD) #read output
					{
					$cmdOutBuffer.=$buffer;
					print $CMDEXP $buffer; #export to file as is
					}
				elsif($FH==$CMDERR) #read error
					{$cmdErrBuffer.=$buffer;}
				else {die "ERR: Extra filehandle detected $FH ";}
				}
			}
		}
	waitpid($cmdPid,0);
	close($CMDEXP);
	
	### Time stats after package execution
	my $timeCmdExec=time();
	warn TimeStamp($verbose)."DBG: STAT: Time to execute $packageName: ".TimeDiffPretty($timeSub0,$timeCmdExec)."\n" if($debug);
	
	### Error handling and output parsing for package command
	# my %seqMatches; #stores seqID:taxa->taxID1,taxID2,... and seqID->taxID->refID1,...
	warn TimeStamp($verbose)."DBG: Command return status: \'$?\'\n" if($debug); #uncomment to see package command exit status for each execution
	if($?!=0) #if exit code non-zero, print error
		{
		warn TimeStamp($verbose)."ERR: Failed to run $packageName properly on fastq $fqFileCount! Exiting.\n";
		warn TimeStamp($verbose)."WARN: INFO: $packageName command used was: `$cmd`\n" if($verbose);
		warn TimeStamp($verbose)."WARN: $packageName said:\n---\n$cmdErrBuffer---\n" if($verbose);
		exit 4;
		}
	else #otherwise, parse output from command
		{
		### Convert bulk output buffer in an array for parsability and destroy scalar
		my @cmdOut;
		@cmdOut=split(/\n/,$cmdOutBuffer); #breaking the buffer this way strips line endings as well
		undef $cmdOutBuffer; #destroy scalar to save memory
		
		### Parse classification output to store assigned taxonomy ID(s)
		my $header=shift(@cmdOut); #remove headers from output
		# warn TimeStamp($verbose)."DBG: header from $packageName: \'$header\'\n" if($debug); #uncomment to see package output header on each execution
		warn TimeStamp($verbose)."DBG: Number of classification records from $packageName: ".scalar(@cmdOut)."\n" if($debug);
		foreach(0..$#cmdOut)
			{
			my($seqID,$refID,$taxID,undef,undef,undef,undef,undef)=split(/\t/,shift(@cmdOut));
			### Store taxonomy as just taxIDs or taxIDs->refIDs
			if($taxID==9606) #store exact reference refIDs for specific taxonomies
				{
				### The if statement above can be changed to include more taxIDs
				### 9606=HomoSapiens
				unless(defined($$seqTaxaRef{$seqID}{$taxID})) #assign the refID if not established
					{$$seqTaxaRef{$seqID}{$taxID}=$refID;}
				else #otherwise append
					{$$seqTaxaRef{$seqID}{$taxID}.=",$refID";}
				}
			else #otherwise store just the taxIDs
				{
				unless(defined($$seqTaxaRef{$seqID})) #assign the taxID if not established
					{$$seqTaxaRef{$seqID}{"taxa"}=$taxID}
				else #otherwise append
					{$$seqTaxaRef{$seqID}{"taxa"}.=",$taxID"}
				}
			}
		}
	my $uniqFqCount=scalar(%$seqTaxaRef);
	warn TimeStamp($verbose)."DBG: Number of unique classification records from $packageName: $uniqFqCount \n" if($debug);
	
	### Time stats after parsing package output
	my $timeCmdParse=time();
	warn TimeStamp($verbose)."DBG: STAT: Time to parse $packageName output: ".TimeDiffPretty($timeCmdExec,$timeCmdParse)."\n" if($debug);
	
	return ($uniqFqCount); #return number of fq seq records processed
	}

### No need for the following function at this time. May be incompletely written.
# sub IsReadableFile #Checks if a file exists and is readable. Terminate 
	# {
	#### Order of Input: filePath-req; verboseLevel-opt
	#### Order of Output: fileStatus (0:doesNotExist, 1:exists, 2:existsAndReadable)
	# my($filePath,$verbose,@rest)=@_;
	# $verbose//=1;
	# die TimeStamp($verbose)."ERR: No filePath specified." unless($filePath);
	
	# warn TimeStamp($verbose)."INFO: Running function IsReadableFile\n" if($verbose>3);
	
	# my $fileStatus=0;
	# $fileStatus++ if(-e $filePath);
	# $fileStatus++ if(-r $filePath);
	# warn TimeStamp($verbose)."WARN: File verification failed \n" if($verbose>3);
	# return($fileStatus);
	# }

sub IsValidDirectory #Checks if the directory exists and is empty
	{
	### Order of Input: dirPath-req; verboseLevel-opt
	### Order of Output: directoryStatus (0:doesNotExist, 1:existsAndEmpty, 2:existsWithContent)
	my($dirPath,$verbose,@rest)=@_;
	$verbose//=1;
	warn TimeStamp($verbose)."INFO: Executing function IsValidDirectory\n" if($verbose>2);
	die TimeStamp($verbose)."ERR: No dirPath specified\n$!" unless($dirPath);
		
	unless(-d $dirPath)
		{warn TimeStamp($verbose)."INFO: Path doesn't exist: $dirPath\n" if($verbose>1);return (0);}
	
	opendir(my $DH,$dirPath) or die TimeStamp($verbose)."ERR: Can't open $dirPath \n$!";
	my $dirContent=scalar(grep {$_ ne "." && $_ ne ".."} readdir($DH)); #count anything except . and ..
	unless($dirContent)
		{warn TimeStamp($verbose)."INFO: Path exists but is empty, $dirPath\n" if($verbose>1);return (1);}
	else
		{warn TimeStamp($verbose)."WARN: Path exists with contents, $dirPath\n" if($verbose);return (2);}
	}

sub PrintHeader #Prints \t spaced header column names to STDOUT (or a fileHandle) with a '#' apended before names
	{
	### Order of Input: \@columnNames-req (a reference to array); $fileHandle-opt
	### Order of Output: none
	my @colNames=@{$_[0]};
	die "ERR: No column names provided.\n$!" unless(@colNames);
	my $FH=$_[1];
	
	my $header="#".join("\t#",@colNames)."\n";
	
	if(defined $FH) #print to filehandle if provided
		{print $FH $header;}
	else #default to STDOUT
		{print $header;}
	}

sub TimeDiffPretty #Provides difference in timepoints in a pretty format
	{
	### Order of Input: startTime-req, endTime-opt
	### Order of Output: timeDiffPretty
	my ($t0,$t1,@rest)=@_;
	die "ERR: No startTime provided.\n$!" unless($t0);
	$t1//=time();
	
	my $diffTime=Time::Seconds->new(int(($t1-$t0)+0.5));
	my $spentTime=$diffTime->pretty;
	return $spentTime;
	}

sub TimeStamp #Provides the current time in an organized format
	{
	### Order of Input: verboseLevel-opt
	### Order of Output: "[formattedCurrentTimeAsString] "
	my ($verbose,@rest)=@_;
	$verbose//=1;
	
	## Find current time, format according to verbosity level
	my $currentTime=time();
	my $formattedTime=strftime("%Y%m%d %H:%M:%S", localtime $currentTime);
	if($verbose>=3) #add microseconds if verbose 3
		{$formattedTime.=sprintf(".%06d", ($currentTime-int($currentTime))*1000000);}
	elsif($verbose>=2) #add miliseconds if verbose 2
		{$formattedTime.=sprintf(".%03d", ($currentTime-int($currentTime))*1000);}
	$formattedTime="[$formattedTime] ";
	return($formattedTime);
	}

