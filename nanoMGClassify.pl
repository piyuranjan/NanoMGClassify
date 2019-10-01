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
## A script for accurate metagenomic classification of long-read Nanopore DNA
## sequences.
##
## Input: FILL INPUT DETAILS
##
## Output: FILL OUTPUT DETAILS
##
## Dependencies: MENTION DEPENDENCIES POSSIBLY WITH VERSION
##
## Author: piyuranjan\@gmail.com
## Source: https://github.com/piyuranjan/PerlScripts/blob/master/nanoMGClassify.pl
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
 -f|force    Force overwrite outFile if it exists. Depends on -o|outFile.
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
 3  Problem with a source FastQ file
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
my($blIdx,$cfIdx,$threads,$trackAtSeq,$workDir);
$threads=1; #default thread/proc count is 1
$trackAtSeq=50; #track progress at every 50 fastq records
$workDir=realpath()."/NanoMGClassify"; #default work directory is PWD/NanoMGClassify (-yyyymmddhhmmss appended later) with full path resolved

### Pre-process generic arguments and defaults
my($debug,$quiet,$help,$verbose); #generic arguments
$verbose=1; #set default verbose to 1

### Read-in arguments
unless(GetOptions(
				#script specific arguments
				'b|blIdx=s' => \$blIdx,
				'c|cfIdx=s' => \$cfIdx,
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

### Time stats after option reading
my $timeReadOptions=time();
warn TimeStamp($verbose)."DBG: Time to read options: ".TimeDiffPretty($time0,$timeReadOptions)."\n" if($debug);


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
warn TimeStamp($verbose)."DBG: Time to check dependencies: ".TimeDiffPretty($timeReadOptions,$timeCheckDependency)."\n" if($debug);


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
warn TimeStamp($verbose)."STAT: Found ".scalar(@fqgzFiles)." file(s) to classify\n" if($verbose);
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
	# warn TimeStamp($verbose)."INFO: Success\n" if($verbose>3);
	}
warn TimeStamp($verbose)."INFO: Succesfully created local fastq for ".scalar(@fqFiles)." file(s).\n" if($verbose>2);
if($verbose>3) #list all local fastq created
	{warn TimeStamp($verbose)."INFO: List of all files created:\n"; warn "$_\n" foreach(@fqFiles);}

### Time stats after creating local fastq
my $timeLocalFastq=time();
warn TimeStamp($verbose)."DBG: Time to create local fastq: ".TimeDiffPretty($timeCheckDependency,$timeLocalFastq)."\n" if($debug);


## Execute centrifuge and blast on each fastq ##

### Prepare work directories for centrifuge and blast
warn TimeStamp($verbose)."STEP: Starting to run Centrifuge and BLAST on each fastq.\n" if($verbose);
my $cfDir=$workDir."/01.CentrifugeClassification"; #this directory will hold all centrifuge outputs
mkdir $cfDir or die $!;
warn TimeStamp($verbose)."INFO: Succesfully created directory for Centrifuge output: $cfDir\n" if($verbose>2);
my $blDir=$workDir."/02.BlastClassification"; #this directory will hold all blast outputs
mkdir $blDir or die $!;
warn TimeStamp($verbose)."INFO: Succesfully created directory for Blast output: $blDir\n" if($verbose>2);

### Cycle on each fastq for centrifuge and blast
my $fqFileCount=0;
foreach my $fqFile(@fqFiles)
	{
	my $fqT0=time();
	$fqFileCount++;
	# warn TimeStamp($verbose)."INFO: Classifying sequences in fastq $fqFileCount: $fqFile\n" if($verbose>2);
	warn TimeStamp($verbose)."INFO: Classifying sequences in fastq $fqFileCount: ".basename($fqFile)."\n" if($verbose>1);
	
	### Prepare centrifuge and blast output file names
	my $fileName=$fqFile;
	(undef,undef,$fileName)=File::Spec->splitpath($fqFile); #extract .fq fileName from path
	$fileName=~s/\.f(ast)?q.*$//; #delete extension
	my $cfFileBase=File::Spec->catfile($cfDir,$fileName); #centrifuge fileBase with full path
	my $cfReport=$cfFileBase."_centrifugeReport.txt"; #centrifuge report fileName
	my $cfClassify=$cfFileBase."_centrifugeClassify.txt"; #classification output fileName with full path
	my $blMatches=File::Spec->catfile($blDir,$fileName."_blastMatches.txt"); #blast top hits fileName with full path
	
	### Execute centrifuge on the fastq and store output
	my $cfCmd="centrifuge -p $threads -x $cfIdx -U $fqFile --report-file $cfReport";
	warn TimeStamp($verbose)."INFO: Executing command:\n`$cfCmd`\n" if($verbose>2);
	my @cfOut=`$cfCmd`;
	warn TimeStamp($verbose)."DBG: Command return status: \'$?\'\n" if($debug);
	if($?!=0){warn TimeStamp($verbose)."ERR: Failed to run Centrifuge.\n";exit 4;}
	warn TimeStamp($verbose)."INFO: Success\n" if($verbose>3);
		
	### Parse classification output to store assigned taxonomy ID(s)
	open(my $CFOUT,">$cfClassify") or die $!; #open file to write the classification output as is
	print $CFOUT shift(@cfOut); #remove header line and print to file
	warn TimeStamp($verbose)."DBG: Number of classification records from Centrifuge: ".scalar(@cfOut)."\n" if($debug);
	my %seqTaxa; #stores seqID->taxID1,taxID2,...
	foreach (0..$#cfOut)
		{
		my $buffer=shift(@cfOut);
		print $CFOUT $buffer; #print the classification output as it is on file
		my($seqID,undef,$taxID,undef,undef,undef,undef,undef)=split(/\t/,$buffer);
		unless(defined($seqTaxa{$seqID})) #if taxID for a seqID is not established, assign the taxID
			{$seqTaxa{$seqID}=$taxID;}
		else #otherwise append the taxID with already existing taxID
			{$seqTaxa{$seqID}.=",$taxID";}
		}
	close($CFOUT);
	warn TimeStamp($verbose)."DBG: Number of unique classification records from Centrifuge: ".scalar(%seqTaxa)."\n" if($debug);
	
	### Prepare fasta file buffer from fastq and feed to blast with taxonomy
	warn TimeStamp($verbose)."INFO: Starting to buffer fastq to fasta and feeding Blastn.\n" if($verbose>2);
	tie my @fqReader, 'Tie::File', $fqFile, autochomp => 0, memory => 0, mode => O_RDONLY or die $!; #open fastq as an array
	open(my $BLOUT,">$blMatches") or die $!; #open file to write blast top hits output as is
	my $fqMax=$#fqReader; #store max lines to avoid recaching
	my $fqCount=0;
	my $fqLine=0;
	while($fqLine<=$fqMax)
		{
		### Prepare fasta buffer
		my ($seqID)=$fqReader[$fqLine]=~/^@(.*?)\s/; #isolate seqID
		unless($seqID) #check for fastq format problems
			{warn TimeStamp($verbose)."ERR: Format corruption. Files are probaly not fastq format.\n";exit 3;}
		my $faBuffer=">$seqID\n".$fqReader[$fqLine+1]; #read sequence and store as fasta buffer
		$fqLine+=4; #advance line count to next fastq record
		# warn TimeStamp($verbose)."DBG: Fasta buffer:\n$faBuffer\n" if($debug); #uncomment to see how sequences are being parsed
		
		### Construct blastn command
		my $blCmd="blastn -db $blIdx -evalue 0.01 -num_threads 1 -outfmt \"6 qaccver saccver staxid pident length mismatch gapopen qstart qend sstart send evalue bitscore\"";
		$blCmd.=" -taxids $seqTaxa{$seqID}" if($seqTaxa{$seqID}); #appends taxonomy data for blastn search if the taxa ID is non-zero
		
		### Execute blastn command and parse output
		my($blWriter,$blReader,$blErr);
		warn TimeStamp($verbose)."DBG: Executing command:\n`$blCmd`\n" if($debug); #uncomment to see blast command for each sequence
		my $blPid=open3($blWriter,$blReader,$blErr,$blCmd) or die "Couldn't fork blastn: $!";
		print $blWriter $faBuffer; #feed in the sequence to blastn
		close($blWriter); #finish feeding the sequence
		###################################### needs work!
		###error handler for blast.
		while(<$blReader>) #export blast output
			{print $BLOUT $_;}
		waitpid($blPid,0); #wait in perl parent until blastn child closes; reduces possibility of a zombie process
		warn TimeStamp($verbose)."DBG: Command return status: \'$?\'\n" if($debug);
		if($?!=0){warn TimeStamp($verbose)."ERR: Failed to run Blast.\n";exit 5;}
		# warn TimeStamp($verbose)."DBG: Success\n" if($debug);
		
		$fqCount++; #count successful fastq record executions.
		
		### Status of fastq progress
		unless($fqCount%$trackAtSeq)
			{warn TimeStamp($verbose)."INFO: Sequences classified from fastq $fqFileCount: $fqCount \n" if($verbose>1);}
		}
	close($BLOUT);
	untie @fqReader; #close fastq filehandle
	
	warn TimeStamp($verbose)."STAT: Finished sequence classification in ".TimeDiffPretty($fqT0)." from fastq $fqFileCount.\n" if($verbose);
	}



## Print stats of operation ##
# my $diffTime=Time::Seconds->new(int((time()-$time0)+0.5));
# my $spentTime=$diffTime->pretty;
warn TimeStamp($verbose)."STAT: Finished everything in ".TimeDiffPretty($time0)."\n" if($verbose);


#########################
###### Subroutines ######
#########################

sub ExampleFuntion #BRIEF DESCRIPTION OF THE FUNCTION
	{
	### Order of Input: DESCRIBE THE INPUTS AND THEIR ORDER.
	### Order of Output: DESCRIBE THE OUTPUTS AND THEIR ORDER.
	}

sub CheckBlast #Checks if the NCBI Blast package is available on PATH, returns good/bad status
	{
	### Order of Input: verboseLevel-optional; versionOfBlastToCheck-optional
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
	### Order of Input: verboseLevel-optional; versionOfCentrifugeToCheck-optional
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
	### Order of Input: verboseLevel-optional; versionOfParallelToCheck-optional
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

### No need for the following function at this time. May be incompletely written.
# sub IsReadableFile #Checks if a file exists and is readable. Terminate 
	# {
	#### Order of Input: filePath-required; verboseLevel-optional
	#### Order of Output: fileStatus (0:doesNotExist, 1:exists, 2:existsAndReadable)
	# my($filePath,$verbose,@rest)=@_;
	# $verbose//=1;
	# die TimeStamp($verbose)."ERR: No filePath specified\n" unless($filePath);
	
	# warn TimeStamp($verbose)."INFO: Running function IsReadableFile\n" if($verbose>3);
	
	# my $fileStatus=0;
	# $fileStatus++ if(-e $filePath);
	# $fileStatus++ if(-r $filePath);
	# warn TimeStamp($verbose)."WARN: File verification failed \n" if($verbose>3);
	# return($fileStatus);
	# }

sub IsValidDirectory #Checks if the directory exists and is empty
	{
	### Order of Input: dirPath-required; verboseLevel-optional
	### Order of Output: directoryStatus (0:doesNotExist, 1:existsAndEmpty, 2:existsWithContent)
	my($dirPath,$verbose,@rest)=@_;
	$verbose//=1;
	warn TimeStamp($verbose)."INFO: Executing function IsValidDirectory\n" if($verbose>2);
	die TimeStamp($verbose)."ERR: No dirPath specified\n$!" unless($dirPath);
		
	unless(-d $dirPath)
		{warn TimeStamp($verbose)."INFO: Path doesn't exist: $dirPath\n" if($verbose>1);return (0);}
	
	opendir(my $DH,$dirPath) or die TimeStamp($verbose)."ERR: Can't open $dirPath\n".$!;
	my $dirContent=scalar(grep {$_ ne "." && $_ ne ".."} readdir($DH)); #count anything except . and ..
	unless($dirContent)
		{warn TimeStamp($verbose)."INFO: Path exists but is empty, $dirPath\n" if($verbose>1);return (1);}
	else
		{warn TimeStamp($verbose)."WARN: Path exists with contents, $dirPath\n" if($verbose);return (2);}
	}

sub PrintHeader #Prints \t spaced header column names to STDOUT (or a fileHandle) with a '#' apended before names
	{
	### Order of Input: \@columnNames-required (a reference to array); $fileHandle-optional
	### Order of Output: none
	my @colNames=@{$_[0]};
	die "ERR: No column names provided.\n$!" unless(@colNames);
	my $FH=$_[1];
	
	if(defined $FH) #print to filehandle if provided
		{print $FH "#$_\t" foreach @colNames;say "";}
	else #default to STDOUT
		{print "#$_\t" foreach @colNames;say "";}
	}

sub TimeDiffPretty #Provides difference in timepoints in a pretty format
	{
	### Order of Input: startTime-required, endTime-optional
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
	### Order of Input: verboseLevel-optional
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

