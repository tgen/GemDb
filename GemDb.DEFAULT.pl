#!/usr/bin/perl -w
#/*************************************************************************
# *
# * TGen CONFIDENTIAL
# * __________________
# *
# *  [2010] - [2015] Translational Genomics Research Institute (TGen)
# *  All Rights Reserved.
# *
# * NOTICE:  All information contained herein is, and remains
# * the property of Translational Genomics Research Institute (TGen).
# * The intellectual and technical concepts contained herein are proprietary
# * to  TGen and may be covered by U.S. and Foreign Patents,
# * patents in process, and are protected by trade secret or copyright law.
# * Dissemination of this information, application or reproduction
# * is strictly forbidden unless prior written permission is obtained
# * from TGen.
# *
# * Major Contributor(s):
#    David Craig
#    Release 09/15/15
#
#  Dependencies
#     MongoDB,Time::localtime,File::stat, Scalar::Util, Storable, FindBin,
#/
#######   Determines which Configuration Library is needed ###########################

##### DEFAULTS ######
our $runType = "DEFAULT";
our $VERSION = '3.91';
##### DEFAULTS ######


######################################################################################
###                     These should not need to be changed                       ####
###                      Please edit conf file                                    ####
######################################################################################
$| = 1;
our $GemStatus = 0;
use FindBin '$Bin';
use Time::localtime;
require "$Bin/lib/Maintain.$VERSION.pm";
require "$Bin/lib/InsertVar.$VERSION.pm";
require "$Bin/lib/Annotate.$VERSION.pm";
require "$Bin/lib/Heuristics.$VERSION.pm";
our $date          = localtime->year * 365 * 3600 * 24 + localtime->yday * 3600 * 24 + localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
our $start         = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
our $stdout        = *STDOUT;
our $stderr        = *STDERR;
our $diff          = $start-$start;
our $RunParameters = { 'runDir' => $Bin, 'version' => $VERSION, 'binDir' => $Bin };
chdir( $RunParameters->{'runDir'} );


######################################################################################    
###                      MAIN BLOCK                                               ####
######################################################################################

{  
##### HELP PRINTOUT######
if (join(@ARGV)=~/--h/) {
    print <<OUT
    _____________GemDb3 $VERSION, Copyright 2015 TGen, David W. Craig_____________________
    Requires "GemDb3.PROD.conf.pm" which contains run parameters (RunParameters)
    Example Commands:
         ./GemDb3.pl     
         ./GemDb3.pl threads=16   
         ./GemDb3.pl runType=SU2C    
         ./GemDb3.pl clean=1 
        
    Options:   
          clean=1       | Removes all partial or incomplete runs, and retries
          nolog=1       | No log files (print to STDOUT/STDERR      
          runType=SU2C  | Parameter file is GemDb3.SU2C.conf.txt, default runType is PROD
          CreateOnPath= | Full Path for Deletion of all contained VCFs, then adding VCFs
          threads=16    | Start with 16 threads
      
    Examples:
          ./GemDb3.pl                                  |  Rebuilds from scratch
          ./GemDb3.pl clean=1                          |  Cleans from crash
          ./GemDb3.pl merge=1                          |  Consolidates records
          ./GemDb3.pl AddOnly=su2c                     |  Adds specific files with 'su2c' in name  
          ./GemDb3.pl ForcePath=su2c runType=SU2C      |  Removes any entrys and reinserts with new runType    
          ./GemDb3.pl CreateOnPath=/reports/C017/0001/20150609 /T3/A1STX/analysisResults/ 

    Files For Communication:
          PROD.RUNNING        |  At least one sample has been added.  If not present, database is purged.
          PROD.pid093.ADDING  |  Markers being added by process PID 123
          PROD.INACTIVE       |  No records are being added.    
          PROD.TOCLEAN        |  Flags that records have been added, and clean should start when nothing left to add   
          PROD.pid123.CLEANING|  Cleaning process is underway    
          PROD.INACTIVE       |  No records are being added.    
    -------END OF HELP-----


___________ STARTING GemDb3 $VERSION, Copyright 2015 TGen, David W. Craig ______
            #######-- Configuration File is GemDb3.$runType.conf.pm -#########
_______________________________________________________________________________   

OUT
      ;
}
###############################################################    
## Clean OLD files
###############################################################    
    unless (-e "$Bin/log") {system("mkdir $Bin/log")}
    unless (-e "$Bin/runDir") {system("mkdir $Bin/runDir")}
    system("find ./log/ -name \"GemDb3.$runType.*log\" -mtime +1 -exec cat {} >> ./log/$runType.old.log \\;");
    system("find ./log/ -name \"GemDb3.$runType.*log\" -mtime +1 -exec rm {} \\;");
    print <<OUT
++---------------------------------------------++
     Log printing to $RunParameters->{'runDir'}/log/GemDb3.$runType.$date.pid$$.log. 
     To print log to STDOUT use -nolog=1.

     +For help use --h.  
       To suspend and push to background use control-z and 'bg;disown;'
       To clean crashed runs, use ./GemDb.pl clean=1 or rm runDir/DEFAULT.pid*
++---------------------------------------------++
OUT
;
    open my $LOG, ">", "$RunParameters->{'runDir'}/log/GemDb3.$runType.$date.pid$$.log" or die "Can't open $RunParameters->{'runDir'}/GemDb3.$date.log";
    $RunParameters->{'LOG'} = \$LOG;

###############################################################    
##### LOAD COMMANDLINE ARGS ######
###############################################################
  $RunParameters->{'threads'}=1;
  ARGS: foreach my $arg (@ARGV) {
        chomp($arg);
        if ( $arg =~ /(.*?)=(.*)/ ) {
            if ( length($1) > 0 && length($2) > 0 ) {
                chomp($2);
                $RunParameters->{$1} = $2;
                $RunParameters->{'StartupRunParameters'}->{$1} = $2;
            }
        }
    }
    $RunParameters->{'StartupRunParameters'}->{'version'}=$VERSION;
    if ( exists( $RunParameters->{runType} ) ) { $runType = $RunParameters->{runType}; }
    print "\n\n\n++------------------------------------------------++\n";
    print "++---GemDb3 $VERSION started $date and pid:$$---++\n";
    print "++------------------------------------------------++\n\n\n";    
###############################################################
##### Threading
###############################################################
    if ($RunParameters->{'threads'} > 1) {
        print "\t+Starting Threading\n";
        for ($j=1;$j<=$RunParameters->{'threads'};++$j) {
            sleep(10);
            my $username = getpwuid( $< );
            $procCount=`ps -ef | grep GemDb | grep perl | grep $username | wc -l`; chomp($procCount); 
            --$procCount;
            if ($procCount >= $RunParameters->{'threads'}) { 
                 print "\t$procCount Threads already started, $RunParameters->{'threads'} requested.\n";
                 goto END;
            }
            print "\t+Forking thread $j, $procCount threads running for user $username\n";
            system("nohup $Bin/GemDb3.$runType.pl >& /dev/null &");
        }
        goto END;
    }
    unless ( join( ":", @ARGV ) =~ /NoLog/i ) { *STDERR = $LOG; *STDOUT = $LOG; }
###############################################################
##### LOAD RunParamaters######
###############################################################
    print "\n\n #######---Configuration File is GemDb3.$runType.conf.pm--#########\n";
    system("rm -f $RunParameters->{'runDir'}/$runType.INACTIVE");
    if ( -e "$Bin/GemDb3.$runType.conf.pm" ) {
        require "$Bin/GemDb3.$runType.conf.pm";
        $started                = "$runType.RUNNING";
        $RunParameters          = Conf->loadDefaults($RunParameters);
        $RunParameters->{'LOG'} = \$LOG;
        foreach $arg (@ARGV) {
            chomp($arg);
            if ( $arg =~ /(.*?)=(.*)/ ) { $RunParameters->{$1} = $2 }
        }
        Maintain->mongoConnect($RunParameters);
    }
    else { die "!!!!Can't find GemDb3.$runType.conf.pm!!!\n"; }
    if (exists($RunParameters->{'maxProcesses'})) {
       $procCount=`ps -ef | grep GemDb | grep perl | wc -l`; chomp($procCount);
       if ($procCount>=$RunParameters->{'maxProcesses'}) {
              print "\tMax proccesses of $RunParameters->{'maxProcesses'} is met, no more processes\n";
              goto END;
       }
    }

    Maintain->printSystemDefaults($RunParameters);

###############################################################
##### STARTUP ######
###############################################################
    unless ( -e "$RunParameters->{'runDir'}/$started" ) {
        print "\n+++-------------- Replacing a DEFAULT.RUNNING file, did you mean to erase database.  Its still there-------------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$started");
        $RunParameters->{'largeRun'} = int(1);
    #    Maintain->dropDb($RunParameters);
        $GemStatus = 1;
    }
###############################################################
##### CLEAN UP AFTER CRASH ######
###############################################################
    if ( exists( $RunParameters->{'clean'} ) ) {
        print "\n+++-------------Cleaning and Restarting Due to Crash-------+++ \n";
        $RunParameters->{'scrubEnabled'}=1;
        $RunParameters->{'eraseOld'}=1;
        $RunParameters->{'sync'}=1;
       sleep(60);
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
        system("rm -f $RunParameters->{'runDir'}/$runType.*.ADDING");        
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        Maintain->closeAllConnections($RunParameters);
        $GemStatus = 1;        
        Maintain->purgePartials($RunParameters);             
        #$RunParameters->{'resetCount'}=1;
        Maintain->mongoConnect($RunParameters);
        Maintain->buildBufferConnections($RunParameters);
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        if ( InsertVar->Insert($RunParameters) > 0 ) {
            Annotate->annotate( $RunParameters);
            $GemStatus = 1;
            system("echo > $RunParameters->{'runDir'}/$runType.TOCLEAN");
        }
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        Maintain->calculateDbFreq($RunParameters);        
        Heuristics->indexCollections($RunParameters);
        Heuristics->runInheritance($RunParameters,'germline');
        Heuristics->joinVCF($RunParameters);
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
       print "Completed purging partials\n";
    }
    if ( exists( $RunParameters->{'merge'} ) ) {
        print "\n+++------------- Merging-------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        $RunParameters->{'scrubEnabled'}=1;
        $RunParameters->{'resetCount'}=1;
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
        system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
        system("rm -f $RunParameters->{'runDir'}/$runType.*.ADDING");
        Maintain->closeAllConnections($RunParameters);
        $GemStatus = 1;
        Maintain->purgePartials($RunParameters);
        #$RunParameters->{'resetCompact'}=1;
        #$RunParameters->{'resetCount'}=1;
        Maintain->cleanAndMerge($RunParameters);
        Maintain->calculateDbFreq($RunParameters);
        Heuristics->indexCollections($RunParameters);
        Heuristics->joinVCF($RunParameters);
    }
    if ( glob("$RunParameters->{'runDir'}/$runType.pid*") ) { 
        $RunParameters->{'pidRunning'} = 1; 
    } else {
         $RunParameters->{'pidRunning'} = 0;
    }
    
###############################################################
#### INSERTING ######
###############################################################
    if ( -e "$RunParameters->{'runDir'}/$started" && !( glob("$RunParameters->{'runDir'}/$runType.pid*.CLEANING") ) ) {
        print "\n+++--------------GemDb3 Adding Files--------------+++ \n";
        system("echo > $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        Maintain->buildBufferConnections($RunParameters);
        if ( InsertVar->Insert($RunParameters) > 0 ) {
            Annotate->annotate( $RunParameters);
            $GemStatus = 1;
            system("echo > $RunParameters->{'runDir'}/$runType.TOCLEAN");
        }
        system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.ADDING");
        Maintain->closeBufferConnections($RunParameters);
    }

###############################################################
##### CLEANING MAINTAINANCE & FIRST RUN ######
###############################################################

    unless ( glob("$RunParameters->{'runDir'}/$runType.pid*.ADDING") || glob("$RunParameters->{'runDir'}/$runType.pid*.CLEANING") ) {
        if ( -e "$RunParameters->{'runDir'}/$runType.TOCLEAN" ) {
            print "\n+++--------------GemDb3 Count, Indexing--------------+++ \n";
            system("echo > $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
            Maintain->calculateDbFreq($RunParameters);
            Heuristics->indexCollections($RunParameters);
            Heuristics->runInheritance($RunParameters,'germline');
            Heuristics->joinVCF($RunParameters);
            $status=1;
            $GemStatus = 1;
            system("rm -f $RunParameters->{'runDir'}/$runType.pid$$.CLEANING");
            if ($status==1) {
                system("rm -f $RunParameters->{'runDir'}/$runType.TOCLEAN");
            }
        }
        else {
            print "\t+Nothing to Do, calculating DbFreq\n";
            Maintain->calculateDbFreq($RunParameters);
            print "\t+Nothing to Do\n";
        }
    }
    system("echo > $RunParameters->{'runDir'}/$runType.INACTIVE");


###############################################################
##### SPAWN IF NOT FINISHED ######
###############################################################
    $diff = sprintf( "%3.2f", ( localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start ) / 60 );
    if ( $diff < 0 ) { $diff = $diff + 1440 };
    print "\n+++-----GemDb3 Finished--($diff minutes)--------+++ \n";
    close $LOG;
}

#####  EXIT STRATEGY  #####
sleep(2);
END: if ( $GemStatus == 0 ) {
    system("rm  ./log/GemDb3.$runType.$date.pid$$.log");
} elsif ( $diff > 5 ) {
    system("cat ./log/$runType.old.log ./log/GemDb3.$runType.$date.pid$$.log > ./log/$runType.old.log");
    system("rm  ./log/GemDb3.$runType.$date.pid$$.log");
    print $stdout "\t+Spawning new process\n";
    system("$Bin/GemDb3.$runType.pl&");
    print "--Ending This Process pid$$ Gracefully--\n";
} else {
    system("cat ./log/$runType.old.log ./log/GemDb3.$runType.$date.pid$$.log > ./log/$runType.old.log");
    system("rm  ./log/GemDb3.$runType.$date.pid$$.log");
}
