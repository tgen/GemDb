#!/usr/bin/perl
#/*************************************************************************
# *
# * TGen CONFIDENTIAL
# * __________________
# *
# *  [2010] - [2014] Translational Genomics Research Institute (TGen)
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
# * Author (s):
#    David Craig
#    Release 9/1/15
#
#/
######################################################################################
package InsertVar;
use strict;
use Time::localtime;
use File::stat;
use Scalar::Util qw(looks_like_number);
use MongoDB;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);

$| = 1;


###########################
#### GLOBAL VARIABLES #####
###########################
our $lastGlobalTime   = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
our $isExon           = {};
our $isGene           = {};
our $excludedVariants = {};
##########################
#### MAIN ################
##########################
sub Insert {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters->{'collectionDb'} = 'noDb';
    Maintain->reloadStartupRunParameters($RunParameters);
    my $tCounter = 0;
    my $start    = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    $RunParameters->{'GlobalStart'} = $start;
    print "\n\t----------Adding Files--------- \n";
    findFiles->($RunParameters);
    if ( scalar( @{ $RunParameters->{'filepaths'} } ) == 0 ) {
        return 0;
    }
    if ( exists( $RunParameters->{'loadLocations'} ) ) {
        if ( $RunParameters->{'loadLocations'} eq "1" ) {
            ( $isExon, $isGene ) = Maintain->loadLocations($RunParameters);
        }
    } else { print '\t\t+Not loading Ensemble\n'; }
    $RunParameters->{'isGene'} = $isGene;
    print "\t\t+Checking for nonReportableVariants...";
    if ( exists( $RunParameters->{'nonReportableVariants'} ) ) {
        excludeNonReportableVariants($RunParameters);
    }
    else { print "-No excluded variants\n"; }
    if (exists ( $RunParameters->{'largeRun'} )) {
        if ( $RunParameters->{'largeRun'}==0) {
            for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
                my $collectionName=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};   
                $RunParameters->{$collectionName}->ensure_index( { 'biomarker' => 1 } );
                $RunParameters->{$collectionName}->ensure_index( { 'info.toDelete'  => 1 } );
                $RunParameters->{$collectionName}->ensure_index( { 'info.group'     => 1 } );
                $RunParameters->{$collectionName}->ensure_index( { 'variants.filepath' => 1 } );            
            }    
           undef $RunParameters->{'largeRun'};
           delete $RunParameters->{'largeRun'};            
        }    
    } 
    if (exists ( $RunParameters->{'multipledVariant'} )) {
        if ( $RunParameters->{'multipledVariant'}==0) {    
           undef $RunParameters->{'multipledVariant'};
           delete $RunParameters->{'multipledVariant'};   
        }
    }        
    my $temp=scalar(@{ $RunParameters->{'filepaths'} });
    print "\t+ Total files to process: $temp\n";
    if ($temp>35) {
        print "\t\t-Files to add is $temp, switching to largeRun=1\n";
        $RunParameters->{'largeRun'}=1;
    }
  ALLVCFS_LOOP: while ( my $filepath = shift( @{ $RunParameters->{'filepaths'} } ) ) {
        my $temp=scalar(@{ $RunParameters->{'filepaths'} });
        print "\t+ Remaining files to process: $temp \n";  
        my $filedate   = shift( @{ $RunParameters->{'filedates'} } );
        my $curVcfPath = shift( @{ $RunParameters->{'vcfPathsHash'} } );
        $tCounter += processVcfFile( $RunParameters, $filepath, $filedate, $curVcfPath );
        if ( $tCounter > $RunParameters->{'maxInsertsPerRun'} ) { last ALLVCFS_LOOP }
    }
    if ( $tCounter > 0 ) { print "\t+Total Inserted: $tCounter\n"; }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "------ Completed Adding ($diff seconds)----  \n";
    return 1;
}


############## VCF Processing ##########################
sub processVcfFile {
    my ( $RunParameters, $filepath, $filedate, $curVcfPath ) = @_;
    my $counter    = 0;
    my $inserts    = 0;
    my $parsedLine = "";
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    if ( $RunParameters->{'filesCollection'}->count( { 'filepath' => $filepath } ) > 0 ) {
        return 0;
    }
    else {
        $RunParameters->{'filesCollection'}->insert( { 'filepath' => $filepath, 'status' => int(0), 'pid'=>0 }, { safe => 1 } );
        sleep(1);
        my $holder = $RunParameters->{'markersDb'}->run_command(
            {
                findAndModify => 'files',
                query         => { 'filepath' => $filepath, 'pid'=>0 },
                update        => { '$set' => { 'pid' => $$ } },
                'new'         => 1
            }
        );
        my $holderRef = $holder->{'value'};
        unless (exists($holderRef->{'pid'})) {
            print "\t\t+Dup submission:  $filepath\n";
            $RunParameters->{'filesCollection'}->remove({'pid'=>$$,'filepath' => $filepath});
            return 0;        
        }
        unless ( $holderRef->{'pid'} eq $$ ) {
            print "\t\t+Dup submission:  $filepath\n";
            $RunParameters->{'filesCollection'}->remove({'pid'=>$$,'filepath' => $filepath});
            return 0;
        }
        my $cursor=$RunParameters->{'filesCollection'}->find( { 'filepath' => $filepath}, { multiple=>1,safe => 1 } );
        while (my $doc3=$cursor->next){
            if ($doc3->{'pid'} ne $$) {
                 $RunParameters->{'filesCollection'}->remove({'pid'=>$$,'filepath' => $filepath});
                 return 0;                
            }
        }
    }
    sleep(3);
    if ( $RunParameters->{'filesCollection'}->count( { 'filepath' => $filepath } ) > 1 ) {
        $RunParameters->{'filesCollection'}->remove({'pid'=>$$,'filepath' => $filepath});
        return 0;
    }
    my $fileInfoRef = {};
    $fileInfoRef = validateFile( $RunParameters, $filepath, $filedate, $curVcfPath );
    $fileInfoRef->{'info'}->{'version'} = $RunParameters->{'version'}; 
    $RunParameters->{'filesCollection'}->update( { 'filepath' => $filepath }, $fileInfoRef, { safe => 1 } );
    if ( $fileInfoRef->{'valid'} == 0 || $fileInfoRef->{'filetype'} eq "notKnown"  || length($RunParameters->{'uidName'} )<1 || !(defined($RunParameters->{'uidName'})) || !(exists($RunParameters->{'uidName'}))) {
        $RunParameters->{'filesCollection'}->update( { 'filepath' => $filepath }, { '$set' => { 'status' => int(1) } }, { safe => 1 } );
        return 0;
    }
    ( $fileInfoRef, $parsedLine ) = openVcf( $RunParameters, $fileInfoRef, $filepath, $curVcfPath );
    push( @{ $RunParameters->{'addedFilepaths'} }, $filepath );
    push( @{ $RunParameters->{'addedUids'} },      $fileInfoRef->{ $RunParameters->{'uidName'} } );
    my $batchInsertArrayRef = [];
    
  VCF_LOOP: while ($parsedLine) {
        chomp($parsedLine);
        my @temp = split( /\t/, $parsedLine );
        if ( scalar(@temp) < 5 ) { $parsedLine = (<FILE>); next VCF_LOOP; }
        $temp[4] =~ s/,<NON_REF>//g;        
        my @alts=split( /,/, $temp[4] );
        foreach my $alt ( @alts  ) {
            $temp[4] = $alt;
            my (@geneList)   = ();
            my $variantRef = {
                'filepath'=>$fileInfoRef->{'filepath'},
                'md5sum'=>$fileInfoRef->{'md5sum'},
                'assay'=>$fileInfoRef->{'assay'},
                'filename'=>$fileInfoRef->{'filename'},                
                'filetype'=>$fileInfoRef->{'filetype'},                                
                'patient'=>$fileInfoRef->{'patient'},
                'sample'=>$fileInfoRef->{'sample'},
                'study'=>$fileInfoRef->{'study'},
                'varDatabase'=>$fileInfoRef->{'varDatabase'},
                'fileDir'=>$fileInfoRef->{'fileDir'},
                'filepath'=>$fileInfoRef->{'filepath'},
                'vcfPath'=>$fileInfoRef->{'vcfPath'},
                'studyPatient'=>$fileInfoRef->{'studyPatient'},                
                'studyPatientTissue'=>$fileInfoRef->{'studyPatientTissue'},
                'style'=>$fileInfoRef->{'style'},
                'tissue'=>$fileInfoRef->{'tissue'},
                'collectionDb'=>$fileInfoRef->{'collectionDb'},
                "$RunParameters->{'uidName'}"=>$fileInfoRef->{$RunParameters->{'uidName'}}
            };
            $variantRef->{'altCount'}=scalar(@alts); 
            $variantRef->{'altCount'}=int($variantRef->{'altCount'});
            ++$counter;
            if ( $counter % 100000 == 0 || $inserts + 1 % 1000 == 0 ) { print "\t\t+Insert attempts at $counter, Inserts at $inserts for $temp[0]:$temp[1]\n"; }
            my $parsed = join( "\t", @temp );
            if ( parseVcfLine( $variantRef, $RunParameters, $fileInfoRef, $parsed ) == 1 ) {
                if ( exists( $variantRef->{'genes'}->[1] ) ) {
                    @geneList = @{ $variantRef->{'genes'} };
                    @geneList = join( ":", @geneList );
                }
                elsif ( scalar( @{ $variantRef->{'genes'} } ) == 1 ) {
                    @geneList = $variantRef->{'genes'}->[0];
                }
                else {
                    @geneList = 'none';
                }
                foreach my $gene (@geneList) {
                    $variantRef->{'gene'} = $gene;
                    if ($variantRef->{'varDatabase'} eq "PASS" ) { 
                        InsertVar->FilterVariantPipeline( $variantRef, $RunParameters );
                    }
                    validateVar($variantRef);
                    unless ( $variantRef->{'filter'} eq "PASS" || $variantRef->{'filter'} eq "QUALIFIED" || $variantRef->{'filter'} eq "LOWQC" ) {
                        $variantRef->{'varDatabase'} = 'FAIL';
                    }                      
                    if ( $variantRef->{'varDatabase'} eq "PASS" ) {
                        $inserts = $inserts + insertBiomarker( $variantRef, $gene, $RunParameters, $batchInsertArrayRef );
                    }
                }
            }
        }
        $parsedLine = (<FILE>);
    }
    my $stat = {};
    my $rParameter = {};    
    $stat->{'Start_ps'} = `ps -o cmd,%cpu,%mem`;
    $stat->{'Start_df'} = `df -hk`;
    $stat->{'Start_date'} = `date`;        
    $rParameter->{'filetypesFromFilename'}     = $RunParameters->{'filetypesFromFilename'};
    $rParameter->{'filterOutGenes'}            = $RunParameters->{'filterOutGenes'};
    $rParameter->{'reported_aberration_names'} = $RunParameters->{'reported_aberration_names'};
    $rParameter->{'FilterVariantPipeline'}            = $RunParameters->{'FilterVariantPipeline'};
    $rParameter->{'FilterAnnotatePipeline'}            = $RunParameters->{'FilterAnnotatePipeline'};
    $rParameter->{'snpeffValidAnnotations'}            = $RunParameters->{'snpeffValidAnnotations'};
    $rParameter->{'FilterGeneInfoPipeline'}            = $RunParameters->{'FilterGeneInfoPipeline'};    
    $rParameter->{'AnnotateByCoordFields'}            = $RunParameters->{'AnnotateByCoordFields'};    
    $rParameter->{'AnnotateByGeneFields'}            = $RunParameters->{'AnnotateByGeneFields'};    
    $rParameter->{'MongodbHost'}               = $RunParameters->{'MongodbHost'};
    $rParameter->{'AnnotateHost'}            = $RunParameters->{'AnnotateHost'};
    $rParameter->{'StartupRunParameters'}      = $RunParameters->{'StartupRunParameters'};
    $rParameter->{'nonReportableVariants'} = $RunParameters->{'nonReportableVariants'};    
    if ( exists( $RunParameters->{'annotateForceType'} ) ) {
        $rParameter->{'annotateForceType'} = $RunParameters->{'annotateForceType'};
    }    
    $fileInfoRef->{'RunParameters'}        = $rParameter;        
    $fileInfoRef->{'compute'}              = { 'start' => $stat };    
    if ( scalar( @{$batchInsertArrayRef} ) > 0 ) {
        $inserts = $inserts + InsertVar->batchInsertFlush( $RunParameters, $batchInsertArrayRef );
    }
    close(FILE);
    my $end  = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $diff = $end - $start;
    $fileInfoRef->{'vcfCount'}     = int($counter);
    $fileInfoRef->{'insertsCount'} = int($inserts);
    print "\t\t-$RunParameters->{'uidName'}:$fileInfoRef->{$RunParameters->{'uidName'}} InsertAttempts: $inserts in $diff s in $fileInfoRef->{'collectionDb'} (filetype: $fileInfoRef->{'filetype'})\n";
    $RunParameters->{'filesCollection'} ->update( { 'filepath' => $filepath }, { '$set' => { 'vcfCount' => int($counter), 'insertsCount' => int($inserts), 'status' => int(4) } }, { safe => 1 } );
    return $inserts;
}

sub insertBiomarker {
    my $variantRef          = shift;
    my $varGene             = shift;
    my $RunParameters       = shift;
    my $batchInsertArrayRef = shift;
    my $collectionIndex=$RunParameters->{$variantRef->{'collectionDb'} . '-index'};    
    my $connMarkers = $RunParameters->{ $RunParameters->{'collectionDb'} };
    undef( $RunParameters->{''} );
    delete( $RunParameters->{''} );
    undef( $variantRef->{''} );
    delete( $variantRef->{''} );
    $variantRef->{ $RunParameters->{'uidName'} } =~ s/\./_/g;
    my $biomarkerRef = assignToBiomarker( $variantRef, $varGene, $RunParameters );
    my ( $het, $hom, $oth, $all, $ref ) = countGenotypes($variantRef);
    undef( $variantRef->{'genes'} );
    delete( $variantRef->{'genes'} );
    $biomarkerRef->{'variants'}->[0]   = $variantRef;    
    Maintain->assignGroup( $biomarkerRef, $RunParameters );    
    if ( $variantRef->{'varDatabase'} eq "PASS" ) {
        if ( !(exists( $RunParameters->{'largeRun'} ) ) && exists($RunParameters->{'multipledVariant'})) {
            my $cur = $connMarkers->find( Tie::IxHash->new( 'info.group' => $biomarkerRef->{'info'}->{'group'}, 'info.varFilled' => 0  ), { limit => 1 } );    #here
            if ( my $doc = $cur->next ) {
                my $id             = $doc->{'_id'};
                my $hets           = $doc->{'count'}->{'hets'};
                my $homs           = $doc->{'count'}->{'homs'};
                my $lastArrayCount = scalar( @{ $doc->{'variants'} } );
                $connMarkers->update( { '_id' => $id }, { '$push' => { 'variants' => $variantRef } }, { safe => 1 } );
                $connMarkers->update(
                    { '_id' => $id },
                    {
                        '$push' => { $RunParameters->{'UidName'} => $variantRef->{ $RunParameters->{'uidName'} }  
                                   },
                        '$inc' => { 'count.hets' => $het, 'count.all' => $all, 'count.refs' => $ref, 'count.homs' => $hom, 'info.varArray' => 1, 'biomarkerCount' => 1 },
                        '$set' => {'info.checkCount'=>0}
                    },
                    { safe => 1 }
                );
                if ( $doc->{'info'}->{'varArray'} + 1 >= $doc->{'info'}->{'varMax'} ) {
                    $connMarkers->update( { '_id' => $id }, { '$set' => { 'info.varFilled' => int(1) } }, { 'safe' => 1 } );
                }
                return 1;
            }
        }
        $biomarkerRef->{'info'}->{'varArray'}        = int(1);
        $biomarkerRef->{'info'}->{'varMax'}          = $RunParameters->{'varGlobalMax'};
        $biomarkerRef->{'count'}->{'hets'} = $het;
        $biomarkerRef->{'count'}->{'homs'} = $hom;
        $biomarkerRef->{'count'}->{'oths'} = $oth;
        $biomarkerRef->{'count'}->{'all'}  = $all;
        $biomarkerRef->{'count'}->{'refs'} = $ref;

        if ( exists( $variantRef->{'studyPatient'} ) ) {
            $biomarkerRef->{'studyPatient'} = $variantRef->{'studyPatient'};
        }
        $biomarkerRef->{'info'}->{'toDelete'}  = int(0);
        $biomarkerRef->{'info'}->{'checkCount'}  = int(1);  
        $biomarkerRef->{'info'}->{'compact'}  = int(0);                      
        $biomarkerRef->{'info'}->{'version'}   = $RunParameters->{'version'};
        $biomarkerRef->{'info'}->{'varFilled'} = int(1);
        InsertVar->forceVariantsType( $biomarkerRef, $RunParameters );        
        if ( $biomarkerRef->{'info'}->{'varArray'} < $biomarkerRef->{'info'}->{'varMax'} && exists( $RunParameters->{'multipledVariant'} ) ) {
            $biomarkerRef->{'info'}->{'varFilled'} = int(0);
        }
        if ( exists( $variantRef->{ $RunParameters->{'uidName'} } ) ) {
            if (length($variantRef->{ $RunParameters->{'uidName'}})==0) {return 0;}
            if (length($RunParameters->{'uidName'})==0) {return 0;}
            undef( $variantRef->{''} );
            delete( $variantRef->{''} );
            $biomarkerRef->{ $RunParameters->{'UidName'} }->[0]=$variantRef->{ $RunParameters->{'uidName'} };
            $biomarkerRef->{ $RunParameters->{'uidName'} } = $variantRef->{ $RunParameters->{'uidName'} };
            if (  $biomarkerRef->{'toAnnotateByCoord'} == 1 ) {
                $RunParameters->{'vcfCollection'}->{'nakedInsert_' . $biomarkerRef->{'collectionDb'}}->insert( $biomarkerRef, { safe => 1 } );
            }
            elsif (  $biomarkerRef->{'toAnnotateByGene'} == 1 ) {
                $RunParameters->{'vcfCollection'}->{'geneInsert_' . $biomarkerRef->{'collectionDb'}}->insert( $biomarkerRef, { safe => 1 } );
            }
            else { die "!!! Shouldn't be here.  This is here $biomarkerRef->{'biomarker'}"; }
            return 1;
        }
    }
    return 0;
}

sub validateFile {
##################################
    # validateFile
##################################
    my ( $RunParameters, $filepath, $filedate, $curVcfPath ) = @_;
    my %fileInfo = ();
    $fileInfo{'status'}   = int(0);
    $fileInfo{'pid'}   = $$; 
    $fileInfo{'filepath'} = $filepath;
    $fileInfo{'filedate'} = $filedate;
    $fileInfo{'vcfPath'}  = $curVcfPath;
    my $md5sum=`md5sum $fileInfo{'filepath'}`;
    chomp($md5sum);
    if ($md5sum=~/^MD5.*= /) {
      $md5sum=~s/MD5.*= //g;
    }else{
       $md5sum =~s/ .*//g;
    }
    $fileInfo{'md5sum'} = $md5sum;
    print "\n\t\t+Checking File $filepath in $curVcfPath\n";
    #$RunParameters->{'filesCollection'}-ensure_index({'md5sum'=>1});    
    my $cursor = $RunParameters->{'filesCollection'}->find({'md5sum'=>$md5sum});
    if (my $doc3=$cursor->next) {
      $fileInfo{'duplicateFile'}=int(1);
      print "\t\t+!WARN: Adding duplicate file:  $fileInfo{'filepath'} to $doc3->{'filepath'} and $fileInfo{'filedate'} and $fileInfo{'vcfPath'}\n";
    }
    $fileInfo{'collectionDb'} = $RunParameters->{'Collections'}->[0]->{'CollectionName'}; 
    $fileInfo{'valid'}        = int(0);
    $fileInfo{'filetype'}     = 'notKnown';
    if ( $fileInfo{'filepath'} =~ /(.*\/)(.*?)$/ ) {
        $fileInfo{'filename'} = $2;
        $fileInfo{'valid'}    = int(0);
        $fileInfo{'fileDir'}  = $1;
        my $filetypeFromFilename = "unknown";
        if ( exists( $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename} ) ) {
            if ( exists( $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'regex'} ) ) {
                my $regex = $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'regex'};
                if ( $fileInfo{'filename'} =~ /$regex/i ) {
                    $fileInfo{'filetype'} = $filetypeFromFilename;
                    foreach my $parameter ( keys %{ $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename} } ) {
                        unless ( $parameter eq "regex" ) {
                            $fileInfo{$parameter} = $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{$parameter};
                        }
                    }
                }
            }
        }
      LOOP: foreach my $filetypeFromFilename ( keys %{ $RunParameters->{'filetypesFromFilename'} } ) {
            if ( $filetypeFromFilename eq "unknown" ) { next LOOP }
            if ( exists( $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'regex'} ) ) {
                my $regex = $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'regex'};
                if ( $fileInfo{'filename'} =~ /$regex/i ) {
                    $fileInfo{'filetype'} = $filetypeFromFilename;
                    if (exists($RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'joinGeneInMarkers'})) {
                        if ($RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{'joinGeneInMarkers'}==1) {
                            $fileInfo{'joinGeneInMarkers'}=1;
                        }
                    }                    
                    foreach my $parameter ( keys %{ $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename} } ) {
                        unless ( $parameter eq "regex" ) {
                            $fileInfo{$parameter} = $RunParameters->{'filetypesFromFilename'}->{$filetypeFromFilename}->{$parameter};
                        }
                    }
                    last LOOP;
                }
            }
        }
    }
    $RunParameters->{'collectionDb'} = $fileInfo{'collectionDb'};

    ## Validate uidName
    if ( exists( $RunParameters->{'CreateOnUid'} ) ) { 
        if ( exists( $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} } ) ) {
            $fileInfo{'patient'}                     = $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} };
            $fileInfo{'study'}                       = $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} };
            $fileInfo{ $RunParameters->{'uidName'} } = $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} };
        }  else { die "\t!! Please specify CreateOnUid' from commandline\n";  }
    } else {
        Heuristics->discoverFileFields($RunParameters, $filepath, $filedate, $curVcfPath,\%fileInfo);
    }
    if (exists($RunParameters->{'style'})) {
        if ($RunParameters->{'style'} eq 'reports') {
            if ( $fileInfo{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.+?)_(.+?)_(.+?)_(.+?)_(ps.*)/ ) {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}= $1;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $2;}
                unless(exists($fileInfo{'recordId'})) {$fileInfo{'recordId'}                        = $3;}
                unless(exists($fileInfo{'sampleSubgroupVisitSampleNumber'})) {$fileInfo{'sampleSubgroupVisitSampleNumber'} = $4;}
                if ( exists( $RunParameters->{'reportableGenes'}->{$5} ) ) {
                    unless(exists($fileInfo{'assay'})) {$fileInfo{'assay'} = $5;}
                }
                unless(exists($fileInfo{'modifier'})) {$fileInfo{'modifier'} = $6;}
            }
            elsif ( $fileInfo{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.+?)_(.+?)_(.+?)_(.*)/ ) {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}                           = $1;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $2;}
                unless(exists($fileInfo{'recordId'})) {$fileInfo{'recordId'}                        = $3;}
                unless(exists($fileInfo{'sampleSubgroupVisitSampleNumber'})) {$fileInfo{'sampleSubgroupVisitSampleNumber'} = $4;}
                if ( exists( $RunParameters->{'reportableGenes'}->{$5} ) ) {
                    unless(exists($fileInfo{'assay'})) {$fileInfo{'assay'} = $5;}
                }
            }elsif ( $fileInfo{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.+?)_(.+?)_(.*)/ ) {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}                           = $1;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $2;}
                unless(exists($fileInfo{'recordId'})) {$fileInfo{'recordId'}                        = $3;}
                unless(exists($fileInfo{'misc'})) {$fileInfo{'misc'}                            = $4;    }
            }elsif ( $fileInfo{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.+?)_(.*)/ ) {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}                           = $1;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $2;}
                unless(exists($fileInfo{'tissue'})) {$fileInfo{'tissue'}                          = $3;  }  
            }elsif ( $fileInfo{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.*)/ ) {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}                           = $1;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $2;}            
            }else {
                unless(exists($fileInfo{'study'})) {$fileInfo{'study'}                           = $fileInfo{ $RunParameters->{'uidName'} } ;}
                unless(exists($fileInfo{'patient'})) {$fileInfo{'patient'}                         = $fileInfo{ $RunParameters->{'uidName'} } ;        }
            }
        }
    } else {
        unless (exists($fileInfo{'study'})) {
            $fileInfo{'study'} = "noStudy";
        }
        unless (exists($fileInfo{'patient'})) {
            $fileInfo{'patient'} = "noPatient";
        }
    }
    unless ( exists( $fileInfo{'collectionDb'} ) ) {
        $fileInfo{'valid'}        = int(0);
        $fileInfo{'collectionDb'} = $RunParameters->{'Collections'}->[0]->{'CollectionName'}; 
    }
    if (exists( $RunParameters->{'assay'}) && !(exists($fileInfo{'assay'}))) {
        $fileInfo{'assay'} = $RunParameters->{'assay'};     
    }
    unless (exists($fileInfo{'tissue'})) {
        if (exists($fileInfo{'sampleSubgroupVisitSampleNumber'})) {
            $fileInfo{'tissue'}=$fileInfo{'sampleSubgroupVisitSampleNumber'};
        } else {
            $fileInfo{'tissue'}='unknown';
        }
    }
    $fileInfo{'studyPatientTissue'}=join("_",$fileInfo{'study'},$fileInfo{'patient'},$fileInfo{'tissue'});    
    $fileInfo{ $RunParameters->{'uidName'} } =~ s/\./_/g;
    if ( $fileInfo{'valid'} == 1 ) {
        print "\t+ValidFile=$fileInfo{'filename'}\nvcfPath=$fileInfo{'vcfPath'}\n\t\t+$RunParameters->{'uidName'}=$fileInfo{$RunParameters->{'uidName'}}\n\t\t+filetype:$fileInfo{'filetype'} +collectionDb=$fileInfo{'collectionDb'} assay=$fileInfo{'assay'} +valid: $fileInfo{'valid'}\n\t\t+filepath=$fileInfo{'filepath'}\n";
    }
    return \%fileInfo;    
}

sub assignToBiomarker {
    my $variantRef    = shift;
    my $varGene       = shift;
    my $RunParameters = shift;
    my $biomarkerRef = {};
    $biomarkerRef->{'coord'}  = $variantRef->{'coord'};
    if (exists($variantRef->{'chrPos'})) {
       $biomarkerRef->{'chrPos'}=$variantRef->{'chrPos'};
    }    
    $biomarkerRef->{'biomarker'} = $variantRef->{'coord'};
    $variantRef->{'biomarker'} = $biomarkerRef->{'biomarker'};    
    $biomarkerRef->{'biomarker' . $RunParameters->{'UidName'}} = $variantRef->{'coord'} . ":" . $variantRef->{ $RunParameters->{'uidName'} };
    my ( $chr, $pos, $ref, $alt ) = split( ":", $variantRef->{'coord'} );
    $biomarkerRef->{'toAnnotateByCoord'} = int(1);
    $biomarkerRef->{'toAnnotateByChrPos'} = int(1);  
    unless ( $variantRef->{'filter'} eq "PASS" || $variantRef->{'filter'} eq "QUALIFIED" || $variantRef->{'filter'} eq "LOWQC" ) {
        $variantRef->{'varDatabase'} = 'FAIL';
    }           
    if ( ( $ref eq "A" || $ref eq "T" || $ref eq "C" || $ref eq "G" ) && ( $alt eq "A" || $alt eq "T" || $alt eq "C" || $alt eq "G" ) ) {
        ### SNV
        $biomarkerRef->{'type'}              = 'SNV';
        $biomarkerRef->{'molecule'}          = 'DNA';
    }  elsif ( $alt !~ /\>|\</ ) {
        ### Point Mutation
        $biomarkerRef->{'molecule'}          = 'DNA';
        if ( length($ref) > length($alt) ) {
            #Deletion
            $biomarkerRef->{'type'} = 'smallDeletion';
        }
        elsif ( length($ref) < length($alt) ) {
            #Insertion
            $biomarkerRef->{'type'} = 'smallInsertion';
        }
        elsif ( length($ref) == length($alt) ) {
            #Equivical Block substitution
            $biomarkerRef->{'type'} = 'blockSubstitution';
        }
    } else {
        ### Gene-wise
        $biomarkerRef->{'gene'}              = $varGene;
        $biomarkerRef->{'molecule'}          = 'gene';
        $biomarkerRef->{'toAnnotateByGene'}  = int(1);
        $biomarkerRef->{'toAnnotateByCoord'} = int(0);
        $biomarkerRef->{'toAnnotateByChrPos'} = int(0);                

        if ( $alt eq "<SV>" || $alt eq "<TRA>" ) {
            $biomarkerRef->{'type'}     = 'StructuralVariant';
            $biomarkerRef->{'molecule'} = 'gene';
        }
        elsif ( $alt eq "<EXP>" ) {
            $biomarkerRef->{'molecule'}  = 'gene';
            $biomarkerRef->{'biomarker'} = $varGene;
            if ( looks_like_number( $variantRef->{'LOG2FC'} ) ) {
                if ( $variantRef->{'LOG2FC'} < -0.1 ) {
                    $biomarkerRef->{'type'} = 'UnderExpressed';
                }
                elsif ( $variantRef->{'LOG2FC'} > 0.1 ) {
                    $biomarkerRef->{'type'} = 'OverExpressed';
                }
                else {
                    $biomarkerRef->{'type'}      = 'NoExpressionChange';
                    $variantRef->{'varDatabase'} = 'FAIL';
                }
            }
            else {
                $variantRef->{'varDatabase'} = 'FAIL';
            }
        }
        elsif ( $alt eq "<DEL>" ) {
            $biomarkerRef->{'type'} = 'FocalDeletion';
        }
        elsif ( $alt eq "<INS>" ) {
            $biomarkerRef->{'type'} = 'NovelInsertion';
        }      
        elsif ( $alt eq "<INV>" ) {
            $biomarkerRef->{'type'} = 'Inversion';
        }   
        elsif ( $alt eq "<CNV>" ) {
            $biomarkerRef->{'type'} = 'CnvChange';
        }          
        elsif ( $alt eq "<DUP:TANDEM>" ) {
            $biomarkerRef->{'type'} = 'TandemDuplication';
        }   
        elsif ( $alt eq "<DEL:ME>" ) {
            $biomarkerRef->{'type'} = 'MobileElementionDeletion';
        }          
        elsif ( $alt eq "<INS:ME>" ) {
            $biomarkerRef->{'type'} = 'MobileElementionInsertion';
        }                                    
        elsif ( $alt eq "<DUP>" || $alt eq "<GAIN>" ) {
            $biomarkerRef->{'type'} = 'FocalGain';
        }
        elsif ( $alt eq "<FUSION>" ) {
            $biomarkerRef->{'molecule'}  = 'fusion';
            $biomarkerRef->{'gene'}      = $variantRef->{'FUSED_GENE'};
            $biomarkerRef->{'biomarker'} = $variantRef->{'FUSED_GENE'};
            $biomarkerRef->{'type'}      = 'FusedGenes';
        }
        elsif ( $alt eq "<TRANSCRIPT_VARIANT>" ) {
            $biomarkerRef->{'molecule'}  = 'gene';        
            $biomarkerRef->{'gene'}      = $variantRef->{'TRANSCRIPT'};
            $biomarkerRef->{'biomarker'} = $variantRef->{'TRANSCRIPT'};
            $biomarkerRef->{'type'}      = 'TranscriptVariant';
        }
        elsif ( $alt eq "<INSERT_NUC>" && $ref eq "<DELETE_NUC>" ) {
            $biomarkerRef->{'type'} = 'LargeBlockSubstitution';
        }
        elsif ( $alt eq "<INSERT_NUC>" ) {
            $biomarkerRef->{'type'} = 'LargeInsertion';
        }
        elsif ( $ref eq "<DELETE_NUC>" ) {
            $biomarkerRef->{'type'} = 'LargeDeletion';
        } 
        elsif ( $ref eq "<GENE>" ) {
            $biomarkerRef->{'type'} = 'GeneEvent';
            $biomarkerRef->{'molecule'}  = 'gene';
            $biomarkerRef->{'biomarker'} = $varGene;            
        } elsif ( $alt eq "<NOCALL>" && exists($variantRef->{'LOG2FC'})  && exists($variantRef->{'END'})) {
            if ($variantRef->{'LOG2FC'} > 0) {
                $biomarkerRef->{'type'} = 'FocalGain';
            } 
            if ( $variantRef->{'LOG2FC'}<0 ) {
                $biomarkerRef->{'type'} = 'FocalDeletion';
                $variantRef->{'filter'} = 'LOWQC';                
            } 
        }elsif ( $ref=~/\<.*\>/ || $alt=~/\<.*\>/) {
            $biomarkerRef->{'type'} = 'variant';
        } else {
            $variantRef->{'varDatabase'} = 'FAIL';
        }
    }
    $biomarkerRef->{'collectionDb'}=$variantRef->{'collectionDb'};
    $biomarkerRef->{'vcfPath'}=$variantRef->{'vcfPath'};    
    return $biomarkerRef;
}

sub parseVcfLine {
    ##################################
    # parsedLine
    ##################################
    my ( $variantRef, $RunParameters, $fileInfoRef, $line ) = @_;
    my ( $i, $j ) = 0;
    $variantRef->{'degree'}      = int(0);
    $variantRef->{'varDatabase'} = "PASS";
    chomp($line);
    $line =~ s/\"//g;
    my @vcfFields = split( /\t+/, $line );
    if ( $#vcfFields < 7 ) { return -1 }
    $variantRef->{'chr'} = $vcfFields[0];
    if ( $vcfFields[1] eq "POS" )      { return -1 }
    if ( $vcfFields[1] =~ /[A-Za-z]/ ) { return -1 }
    $variantRef->{'hg19pos'} = int( $vcfFields[1] );
    $variantRef->{'varName'} = $vcfFields[2];
    if ( length( $vcfFields[4] ) > 100 ) {
        my $end = $variantRef->{'hg19pos'} + length( $vcfFields[4] );
        $variantRef->{'END'} = "$end";
        $vcfFields[4] = "<INSERT_NUC>";
    }
    if ( length( $vcfFields[3] ) > 100 ) {
        my $end = $variantRef->{'hg19pos'} + length( $vcfFields[3] );
        $variantRef->{'END'} = "$end";
        $vcfFields[3] = "<DELETE_NUC>";
    }
    $variantRef->{'ref'} = $vcfFields[3];
    $variantRef->{'alt'} = $vcfFields[4];
    $variantRef->{'alt'} =~ s/,<NON_REF>//g;
    if ($variantRef->{'alt'} eq "<NON_REF>") {
       return-1
    }
    if ( length( $variantRef->{'alt'} ) == 0 ) { return -1 }
    $variantRef->{'qual'}   = $vcfFields[5];
    $variantRef->{'filter'} = $vcfFields[6];
    $variantRef->{'filterOrig'} = $vcfFields[6];
    $variantRef->{'infoField'}   = $vcfFields[7];
    $variantRef->{'genes'}  = [];
    if ( $variantRef->{'alt'} eq '.' ) { return -1; }
    $variantRef->{'chr'} =~ s/chr//;
    $variantRef->{'chr'} = "chr" . $variantRef->{'chr'};
    if ( $variantRef->{'chr'} !~ /chr/ ) {
        $variantRef->{'chr'} = join( "", "chr", $variantRef->{'chr'} );
    }
    if ($variantRef->{'chr'} eq "chr23") {
       $variantRef->{'chr'} ='chrX';
    } elsif ($variantRef->{'chr'} eq "chr24") {
       $variantRef->{'chr'} ='chrY';
    } elsif ($variantRef->{'chr'} eq "chr25") {
       $variantRef->{'chr'} ='chrMt';
    }
    if ( exists( $RunParameters->{'geneSpecific'} ) ) {
        if ( $RunParameters->{'geneSpecific'} == 1 ) {
            unless ( $variantRef->{'alt'} =~ /\<|\>/ ) {
                my $trimPos = int( $variantRef->{'hg19pos'} / 10 ) * 10;
                my $lookup  = $variantRef->{'chr'} . ":" . $trimPos;
                unless ( exists( $isExon->{$lookup} ) ) {
                    $variantRef->{'varDatabase'} = "FAIL";
                    return -1;
                }
            }
        }
    }
    unless ( exists( $RunParameters->{'InitVariantStatus'} ) ) {
        $RunParameters->{'InitVariantStatus'} = "PASS";
    }    
    if ( $variantRef->{'filter'} eq "QUALIFIED" ) { $variantRef->{'filter'} = $RunParameters->{'InitVariantStatus'}; $variantRef->{'note'} = "Qualified" }
    $variantRef->{'coord'} = join( ":", $variantRef->{'chr'}, $variantRef->{'hg19pos'}, $variantRef->{'ref'}, $variantRef->{'alt'} );
    $variantRef->{'chrPos'} = join( ":", $variantRef->{'chr'}, $variantRef->{'hg19pos'} );    
    if ( exists( $RunParameters->{'nonReportableVariants'} ) ) {
        if ( exists( $excludedVariants->{ $variantRef->{'coord'} } ) ) { return -1 }
    }
    if ( exists( $vcfFields[9] ) ) {
        my @gtVcfFields = @vcfFields[ 9 .. $#vcfFields ];
        $vcfFields[8] =~ s/\./P/g;
        my @gtSubFieldNames = split( /:/, lc( $vcfFields[8] ) );
        for ( $i = 0 ; $i <= $#gtVcfFields ; ++$i ) {
            my @gtSubFields = split( /:/, $gtVcfFields[$i] );
            if ( $#gtSubFieldNames != $#gtSubFields ) { return -1; }
            my $thash = { 'sampleName' => $fileInfoRef->{'sampleNames'}->[$i] };
            for ( $j = 0 ; $j <= $#gtSubFields ; ++$j ) {
                if ( exists( $gtSubFields[$j] ) ) {
                    if ( $gtSubFieldNames[$j] eq "dp" || $gtSubFieldNames[$j] eq "DP") {
                        if ( $gtSubFields[$j] eq "." ) {
                            $thash->{ $gtSubFieldNames[$j] } = int(0);
                        } elsif ( $gtSubFields[$j] =~ /\D/ ) {
                            $thash->{ $gtSubFieldNames[$j] } = int(0);
                        } else {
                            $thash->{$gtSubFieldNames[$j] } = int( $gtSubFields[$j] );
                        }
                    } else {
                        $thash->{$gtSubFieldNames[$j] } = $gtSubFields[$j];
                    }
                    if ( $gtSubFieldNames[$j] eq "gq" ) {
                        if ( $gtSubFields[$j] eq "." ) {
                            $thash->{ $gtSubFieldNames[$j] } = int(0);
                        }elsif ( $gtSubFields[$j] =~ /\D/ ) {
                            $thash->{ $gtSubFieldNames[$j] } = int(0);
                        }else {
                            $thash->{$gtSubFieldNames[$j] } = int( $gtSubFields[$j] );
                        }
                    } else {
                        $thash->{$gtSubFieldNames[$j] } = $gtSubFields[$j];
                    }
                    if ( $gtSubFieldNames[$j] eq "GT" ||  $gtSubFieldNames[$j] eq "gt") {
                       if ($gtSubFields[$j]=~/(\.|0|1|2|3|4|5|6|7|8|9)+/) {
                          my $c=0;
                          while (my $gen=shift @_) {
                              ++$c;
                              if ($gen eq ".") {
                                $gen=-1;
                              }
                              if (length($gen)==1) {
                                 $thash->{ 'g' . $c}=int($gen);
                              }
                          }
                          if ($gtSubFields[$j]=~/\|/) {
                              $thash->{ 'phased'}=int(0);                              
                          } else {
                              $thash->{ 'phased'}=int(1);                                                        
                          }
                       }
                    }                    
                } else {
                    if ( $gtSubFieldNames[$j] eq "dp" || $gtSubFieldNames[$j] eq "DP") {
                        if ($thash->{ $gtSubFieldNames[$j] } eq ".") {
                             $thash->{ $gtSubFieldNames[$j] }=int(0);
                        }
                        $thash->{ $gtSubFieldNames[$j] } = int( $gtSubFields[$j] );
                    } elsif ( $gtSubFieldNames[$j] eq "gq" ) {
                        if ($thash->{ $gtSubFieldNames[$j] } eq ".") {
                             $thash->{ $gtSubFieldNames[$j] }=int(0);
                        }                    
                        $thash->{ $gtSubFieldNames[$j] } = int( $gtSubFields[$j] );
                    }
                }
            }
            my $pos = $i + 1;
            $variantRef->{'genotype'}->{'s' . $pos} = $thash;
        }
    }
    if ( exists( $variantRef->{'swap'} ) ) { $variantRef = swapGt($variantRef) }
    if ( parseInfo( $RunParameters, $variantRef ) == 1 ) {
        if ( $variantRef->{'varDatabase'} eq "FAIL" ) {
            return -1;
        }
        else { return 1; }
    }
    else { return -1; }
}    #    parseVcfLine (variantRef)=variantRef

sub parseInfo {
    my ( $RunParameters, $variantRef ) = @_;
    my @words = split( /\;+/, $variantRef->{'infoField'} );
  WORD: foreach my $word (@words) {
        my ( $key, $val ) = "";
        if ( $word =~ /=/ ) {
            if ( scalar( $word =~ /=/ ) == 1 ) {
                ( $key, $val ) = split( /=/, $word );
            }
            else { $key = $word; $val = "1"; }
        }
        $key =~ s/\./\_/g;
        if ( $key eq "TYPE" ) { $val =~ s/somatic_//g; }
        if ( $key eq "TYPE" && $val =~ /LOH/ ) { $variantRef->{'varDatabase'} = "FAIL-LOH"; return -1; }
        if ( $key eq "END" ) {
            $variantRef->{$key} = $val;
            if ( $variantRef->{'filetype'} eq "cnvFile" ) {
                my $val2 = abs( int( $val - $variantRef->{'hg19pos'} ) );
                my $key2 = "cnvLength";
                $variantRef->{$key2} = $val2;
                if ( $val2 > 350000000 ) { $variantRef->{'varDatabase'} = "FAIL"; return -1; }
            }
        }
        if (   ( $key eq "gene" || $key eq "GENE" )
            && ( $variantRef->{'alt'} eq "<FUSION>" || $variantRef->{'alt'} eq "<TRANSCRIPT_VARIANT>" || $variantRef->{'alt'} eq "<EXP>" ) )
        {
            if ( $val =~ /\,/ ) {
                push( @{ $variantRef->{'genes'} }, split( /\,/, $val ) );
            }
            else {
                push( @{ $variantRef->{'genes'} }, $val );
            }
        }
        if ( exists( $RunParameters->{'IncludeInfoFields'} ) ) {
            if ( exists( $RunParameters->{'IncludeInfoFields'}->{'ALL'} ) ) {
                if ( $RunParameters->{'IncludeInfoFields'}->{'ALL'} == 1) {
                    $variantRef->{$key} = $val;
                }
            }
            if ( exists( $RunParameters->{'IncludeInfoFields'}->{$key} ) ) {
                if ( $RunParameters->{'IncludeInfoFields'}->{$key} > 0 ) {
                    $variantRef->{$key} = $val;
                } elsif ($RunParameters->{'IncludeInfoFields'}->{$key} ==0 )  {
                    undef $variantRef->{$key};
                    delete $variantRef->{$key};
                }
            }
        }
    }
    return 1;
}    #parseInfo(variantRef)=variantRef;

sub validateVar {
    my ($variantRef) = shift;
    if ( exists( $variantRef->{'qual'} ) ) {
        if ( $variantRef->{'qual'} eq "." || $variantRef->{'qual'} =~ /[[:alpha:]]/ ) {
            $variantRef->{'qual'} = 255;
        }
        $variantRef->{'qual'} = $variantRef->{'qual'};
        MongoDB::force_double( $variantRef->{'qual'} );
    }
    if ( exists( $variantRef->{'AC'} ) ) {
        unless ( looks_like_number( $variantRef->{'AC'} ) ) {
            $variantRef->{'AC'} = 255;
        }
        $variantRef->{'AC'} = int( $variantRef->{'AC'} );
    }
    if ( exists( $variantRef->{'DP'} ) ) {
        unless ( looks_like_number( $variantRef->{'DP'} ) ) {
            $variantRef->{'DP'} = 255;
        }
        $variantRef->{'DP'} = int( $variantRef->{'DP'} );
    }
    if ( exists( $variantRef->{'AR1'} ) && exists( $variantRef->{'AR2'} ) ) {
        unless ( looks_like_number( $variantRef->{'AR1'} ) || looks_like_number( $variantRef->{'AR2'} ) ) {
            $variantRef->{'AR1'} = 1;
            $variantRef->{'AR2'} = 1;
        }
        else {
            $variantRef->{'ARDIFF'} = $variantRef->{'AR2'} - $variantRef->{'AR1'};
            MongoDB::force_double( $variantRef->{'ARDIFF'} );
        }
        MongoDB::force_double( $variantRef->{'AR1'} );
        MongoDB::force_double( $variantRef->{'AR2'} );
    }
    else { undef( $variantRef->{'AR1'} ); delete( $variantRef->{'AR1'} ); undef( $variantRef->{'AR2'} ); delete( $variantRef->{'AR2'} ); }
    if ( exists( $variantRef->{'MQ'} ) ) {
        unless ( $variantRef->{'MQ'} =~ /[[:alpha:]]/ ) {
            $variantRef->{'MQ'} = 255;
        }
        MongoDB::force_double( $variantRef->{'MQ'} );
    }
    undef( $variantRef->{''} );
    #undef($variantRef->{'snpEffVariantRef'});
    undef( $variantRef->{'sampleNames'} );    
    undef( $variantRef->{'info'} );
    undef( $variantRef->{'infoField'} );    
    undef( $variantRef->{'degree'} );
    undef( $variantRef->{'info'} );
    #delete($variantRef->{'snpEffVariantRef'});
    delete( $variantRef->{'sampleNames'} );
    delete( $variantRef->{'EFF'} );    
    delete( $variantRef->{'info'} );
    delete( $variantRef->{'infoField'} );  
    delete( $variantRef->{''} );
    delete( $variantRef->{'degree'} );
    return $variantRef;
}

sub batchInsertFlush {
    my $function            = shift;
    my $RunParameters       = shift;
    my $batchInsertArrayRef = shift;
    my $num                 = scalar( @{$batchInsertArrayRef} );
    $RunParameters->{'vcfCollection'}->{'nakedInsert'}->batch_insert( $batchInsertArrayRef, { safe => 1 } );
    @{$batchInsertArrayRef} = ();
    return $num;
}

############## File Maintenance ##########################
sub findFiles {
    my $RunParameters = shift;
    if ( exists( $RunParameters->{'DeleteOnPath'} ) ) {
        InsertVar->DeleteOnPath($RunParameters);
    }
    if ( exists( $RunParameters->{'CreateOnPath'} ) ) {
        InsertVar->CreateOnPath($RunParameters);
    }
    my ( $vars, $filepath ) = "";
    my (%filepathHash) = ();
    my (%vcfPathHash)  = ();
    my @list = ();
    @{ $RunParameters->{'filepaths'} }    = ();
    @{ $RunParameters->{'filedates'} }    = ();
    @{ $RunParameters->{'vcfPathsHash'} } = ();
    unless ( exists( $RunParameters->{'vcfPath'} ) || exists( $RunParameters->{'vcfPaths'} ) ) { die "\t!Unable to find a vcfPath\n"; }
    if ( exists( $RunParameters->{'vcfPath'} ) ) {
        push( @{ $RunParameters->{'vcfPaths'} }, $RunParameters->{'vcfPath'} );
    }
    if ( exists( $RunParameters->{'vcfPaths'} ) ) {
        for ( my $i = 0 ; $i < scalar( @{ $RunParameters->{'vcfPaths'} } ) ; ++$i ) {
            $RunParameters->{'vcfPath'} = $RunParameters->{'vcfPaths'}->[$i];
            unless ( $RunParameters->{'vcfPath'} =~ /\/$/ ) {
                $RunParameters->{'vcfPath'} = $RunParameters->{'vcfPath'} . "/";
            }
            print "\t+Adding vcfPath: $RunParameters->{'vcfPath'}\n";
#            @list = `find $RunParameters->{'vcfPath'} -name "*.vcf"`;
            if (exists ( $RunParameters->{'mtime'})) {
               @list = `find $RunParameters->{'vcfPath'} -iname "*.vcf" -mtime $RunParameters->{'mtime'}`;
            }else{
                @list = `find $RunParameters->{'vcfPath'} -iname "*.vcf"`;
           }
          LOOP: while ( $filepath = shift(@list) ) {    ## Limiting to add Only
                $filepath = File::Spec->rel2abs($filepath);
                chomp($filepath);
                if ( exists $RunParameters->{'AddOnly'} ) {
                    my $val = $RunParameters->{'AddOnly'};
                    if ( $filepath !~ /$val/ ) { next }
                }             
                if ( !( -e $filepath ) )   { next }
                if ( $filepath =~ /\/\./ ) { next }
                $filepathHash{$filepath} = stat($filepath)->mtime;
                $vcfPathHash{$filepath}  = $RunParameters->{'vcfPath'};
                push( @{ $RunParameters->{'filepaths'} },    $filepath );
                push( @{ $RunParameters->{'vcfPathsHash'} }, $RunParameters->{'vcfPath'} );
                push( @{ $RunParameters->{'filedates'} },    stat($filepath)->mtime );
            }
        }
        checkFiles($RunParameters);
    }
    else {
        die "\t+No vcfPaths to look\n";
    }
    return $RunParameters;
}

sub checkFiles {
    my $RunParameters = shift;
    my (%filepath2filedate) = ();
    my (%vcfpath2filepath)  = ();
    my (%filepath2vcfpath)  = ();
    my (%filepathHash)      = ();
    my $curFiledate;
    my $curVcfpath;
    my $curFilepath;
    while ( my $curFilepath = shift( @{ $RunParameters->{'filepaths'} } ) ) {
        $curFiledate                     = shift( @{ $RunParameters->{'filedates'} } );
        $curVcfpath                      = shift( @{ $RunParameters->{'vcfPathsHash'} } );
        $filepathHash{$curFilepath}      = $curFilepath;
        $filepath2filedate{$curFilepath} = $curFiledate;
        $vcfpath2filepath{$curVcfpath}   = $curFilepath;
        $filepath2vcfpath{$curFilepath}  = $curVcfpath;
    }
    if (exists($RunParameters->{'sync'})) {
       if ($RunParameters->{'sync'}==1) {
         $RunParameters->{'enableScrub'}=1;
       }
    } else {
       $RunParameters->{'sync'}=0;
    }
    if (exists($RunParameters->{'clean'})) {
       if ($RunParameters->{'clean'}==1) {
         $RunParameters->{'eraseOld'}=1;
       }
    }
    my $cursor = $RunParameters->{'filesCollection'}->find();
    $cursor->immortal(1);
    # Remove Variants If Files Are Deleted
  LOOP2: while ( my $doc = $cursor->next ) {
        if ( exists( $doc->{'filepath'} ) && exists( $doc->{'filedate'} ) ) {
            my $mongoFilepath = $doc->{'filepath'};
            my $mongoFiledate = $doc->{'filedate'};
            my $mongoVcf      = $doc->{'vcfPath'};
            if ( exists( $filepathHash{$mongoFilepath} ) ) {
                if ( $filepath2filedate{$mongoFilepath} eq $mongoFiledate ) {
                    if ( $doc->{'status'} == 1 && !( exists( $RunParameters->{'CreateOnPath'} ) ) ) {
                        delete $filepathHash{$mongoFilepath};
                        next LOOP2;
                    }
                    elsif ( exists( $RunParameters->{'clean'} ) ) {
                        print "\t!Unexpected status in clean mode ($doc->{'status'}) Scrubbing:$mongoFilepath\n";
                        InsertVar->scrubFilepath( $mongoFilepath, $doc->{'collectionDb'}, $RunParameters );
                    } 
                    else {
                        delete $filepathHash{$mongoFilepath};
                        next LOOP2;
                    }
                }
            }
            else {
                unless ( exists( $RunParameters->{'AddOnly'} ) ) {
                    if ( exists( $vcfpath2filepath{ $doc->{'vcfPath'} } ) && !( exists( $filepathHash{ $doc->{'filepath'} } ) ) ) {
                        if (exists($RunParameters->{'eraseOld'})) {
                            if  ($RunParameters->{'eraseOld'}==1) {
                                print "\t!Databased file $doc->{'filepath'} has been deleted, scrubbing\n";  
                                if (InsertVar->scrubFilepath( $mongoFilepath, $doc->{'collectionDb'}, $RunParameters )==0) {
                                   delete( $filepathHash{$doc->{'filepath'}}); 
                                   next LOOP2;
                                }
                            }                                                  
                        } else {
                #@            print "\t\t+Warning - Database files exist without matched VCF in vcfPath: $doc->{'filepath'}\n";
                            if ( exists( $RunParameters->{'eraseOld'} ) ) {
                               InsertVar->scrubFilepath( $mongoFilepath, $doc->{'collectionDb'}, $RunParameters );
                            } else {
                               delete( $filepathHash{$doc->{'filepath'}});                             
                            }
                        }   
                    }
                    else {
                        #print "\t\t+This will be added: $doc->{'filepath'}  $doc->{'vcfPath'}\n";
                    }
                }
            }
        }
    }
  LOOP3: foreach my $curFilepath ( keys %filepathHash ) {
        my $cursorFile  = $RunParameters->{'filesCollection'}->find( { 'filepath' => $curFilepath } );
        my $curFiledate = $filepath2filedate{$curFilepath};
        my $curVcfPath  = $filepath2vcfpath{$curFilepath};
        if ( my $doc = $cursorFile->next ) {
            if ( exists( $doc->{'filepath'} ) && exists( $doc->{'filedate'} ) ) {
                if ( $doc->{'filedate'} eq $curFiledate ) {
                    if ( exists( $RunParameters->{'CreateOnPath'} ) ) {
                        print "\t+Scrubbing due to CreateOnPath\n";
                        InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters );
                    }
                    else {
                        next LOOP3;
                    }
                }
                else {
                    if ( exists( $doc->{'filepath'} ) && exists( $doc->{'collectionDb'} ) ) {
                        print "\t!Set ($doc->{'filepath'})) to Scrub because dates are not equal: $doc->{'filedate'}  vs. $curFiledate\n";
                        if (InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters ) ==0) {
                            next LOOP3;
                        }
                    }
                    else {
                        print "strange\n";
                    }
                }
            }
        }
        push( @{ $RunParameters->{'filepaths'} },    $curFilepath );
        push( @{ $RunParameters->{'filedates'} },    $curFiledate );
        push( @{ $RunParameters->{'vcfPathsHash'} }, $curVcfPath );
    }
    return $RunParameters;
}

sub forceVariantsType {
    my $function      = shift;
    my $biomarkerRef    = shift;
    my $RunParameters = shift;
    foreach my $key ( keys %{ $RunParameters->{'IncludeInfoFields'} } ) {
        if ( exists( $RunParameters->{'IncludeInfoFields'}->{$key} ) ) {
            if ( $RunParameters->{'IncludeInfoFields'}->{$key} == 2 ) {
                if ( exists( $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                    unless ( looks_like_number(  $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                         $biomarkerRef->{'variants'}->[0]->{$key} = undef;
                    } else {
                        MongoDB::force_double(  $biomarkerRef->{'variants'}->[0]->{$key}* 1.00 );
                    }
                }
            } elsif ( $RunParameters->{'IncludeInfoFields'}->{$key} == 3 ) {
                if ( exists(  $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                    unless ( looks_like_number(  $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                         $biomarkerRef->{'variants'}->[0]->{$key} = undef;
                    } else {
                         $biomarkerRef->{'variants'}->[0]->{$key} = int( sprintf( "%.0f",  $biomarkerRef->{'variants'}->[0]->{$key} ) );
                    }
                }
            }elsif ( $RunParameters->{'IncludeInfoFields'}->{$key} == 1 ) {
                if ( exists(  $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                     $biomarkerRef->{'variants'}->[0]->{$key} = "$biomarkerRef->{'variants'}->[0]->{$key}";
                }
            }elsif ( $RunParameters->{'IncludeInfoFields'}->{$key} == 0 ) {
                if ( exists(  $biomarkerRef->{'variants'}->[0]->{$key} ) ) {
                    undef(  $biomarkerRef->{'variants'}->[0]->{$key} );
                    delete(  $biomarkerRef->{'variants'}->[0]->{$key} );
                }
            }
        }
    }
}

sub scrubFilepath {
    my $function            = shift;
    my $mongoFilepath       = shift;
    my $variantCollectionDb = shift;
    my $RunParameters       = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    if ( exists( $RunParameters->{'pidRunning'} ) ) { 
        if ($RunParameters->{'pidRunning'}==1) {
           print "\t\t+Another instance of adding variants, aborting scrub\n";
           return 0; 
        }
    }
    if ( exists( $RunParameters->{'sync'} ) ) { 
        if ($RunParameters->{'sync'}==1) {
            $RunParameters->{'scrubEnabled'}=1;
        }
    }
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();    
    if (exists($RunParameters->{'scrubEnabled'} ) ) { 
        if ($RunParameters->{'scrubEnabled'} eq "on" || $RunParameters->{'scrubEnabled'} eq "1" || $RunParameters->{'scrubEnabled'} eq "true") {
            print "\t\t\t+Scrubbing:$mongoFilepath in $variantCollectionDb....";    
            $RunParameters->{'filesCollection'}->update( { 'filepath' => $mongoFilepath }, { '$set' => { 'status' => int(2) } }, { safe => 1 } );
            $RunParameters->{$variantCollectionDb}->ensure_index( { 'variants.filepath' => 1 } );
            $RunParameters->{$variantCollectionDb}->remove( { 'variants.filepath' => $mongoFilepath, 'info.varArray' => 1 }, { safe => 1, 'multiple' => 1 } );
            my $cursor=$RunParameters->{$variantCollectionDb}->find({'variants.filepath'=>$mongoFilepath});
            while (my $doc=$cursor->next) {
                $RunParameters->{$variantCollectionDb} ->update(  {'_id'=>$doc->{'_id'}}, { '$inc' => { 'info.varArray' => -1 } }, { safe =>1 });
                $RunParameters->{$variantCollectionDb}->update( {'_id'=>$doc->{'_id'}}, { '$pull' => { 'variants' => { 'filepath' => $mongoFilepath } } }, { 'safe' => 11  } );
                $RunParameters->{$variantCollectionDb}->remove( { '_id'=>$doc->{'_id'},'variants' => { '$size' => 0 } }, { 'safe' => 1 } );                
            }            
            $RunParameters->{'filesCollection'}->remove( { 'filepath' => $mongoFilepath }, { safe => 1, 'multiple' => 1 } );
            $RunParameters->{'uidCollection'}->update( {}, { '$pull' => { 'filepaths' => { 'filepath' => $mongoFilepath } } }, { 'multiple' => 1, 'safe' => 1 } );
            $RunParameters->{'uidCollection'}->remove( { 'filepaths' => { '$size' => 0 } }, { 'multiple' => 1, 'safe' => 1 } );
            print "Scrubbing empty records...";
            print "Scrub Complete.\n";
            my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
            print "___________Completed scrubbing ($diff seconds)_______  \n";        
        } else {
                print "\t\t+Scrubbing not turned on, use scrubEnabled=1\n";

        }
    } else {
        print "\t\t+Scrubbing not turned on, use scrub=on\n";
    }
    return 1;
}    # scrubFilepath (mongoFilepath)=1

sub openVcf {
    ##################################
    # openVCF
    ##################################
    my ( $RunParameters, $fileInfoRef, $filepath, $curVcfPath ) = @_;
    my $line = "";
    open( FILE, $fileInfoRef->{'filepath'} ) or warn "Can't open file in openVCF: $fileInfoRef->{'filepath'}!\n";
    $fileInfoRef->{'vcfParameters'} = {};
  HEADER: while ( $line = <FILE> ) {
        chomp($line);
        if ( !( $line =~ /^\#/ ) )     { last HEADER }
        if ( $line =~ /##fileformat/ ) { next HEADER }
        if ( $line =~ /##INFO=<ID=(.+?),/ ) {
            my $i = $1;
            $i =~ s/\./\_/g;
            my $id = { 'flag' => 'INFO' };
            if ( $line =~ /Number=(.+?),/ ) {
                $id->{'Number'} = $1;
            }
            if ( $line =~ /Type=(.+?),/ ) {
                $id->{'Type'} = $1;
            }
            if ( $line =~ /Description=(.+?)>/ ) {
                $id->{'Description'} = $1;
                $id->{'Description'} =~s/\"//g;
            }
            $fileInfoRef->{'vcfParameters'}->{'INFO'}->{$i} = $id;
        }
        elsif ( $line =~ /##FORMAT=<ID=(.+?),/ ) {
            my $i = $1;
            $i =~ s/\./\_/g;
            my $id = { 'flag' => 'FORMAT' };
            if ( $line =~ /Number=(.+?),/ ) {
                $id->{'Number'} = $1;
            }
            if ( $line =~ /Type=(.+?),/ ) {
                $id->{'Type'} = $1;
            }
            if ( $line =~ /Description=(.+?)>/ ) {
                $id->{'Description'} = $1;
                $id->{'Description'} =~s/\"//g;                
            }
            $fileInfoRef->{'vcfParameters'}->{'FORMAT'}->{$i} = $id;
        }
        elsif ( $line =~ /\#\#ALT\<ID=(.+?),/ ) {
            my $i = $1;
            $i =~ s/\./\_/g;
            my $id = { 'flag' => 'ALT' };
            if ( $line =~ /Number=(.+?),/ ) {
                $id->{'Number'} = $1;
            }
            if ( $line =~ /Type=(.+?),/ ) {
                $id->{'Type'} = $1;
            }
            if ( $line =~ /Description=(.+?)>/ ) {
                $id->{'Description'} = $1;
                $id->{'Description'} =~s/\"//g;
            }
            $fileInfoRef->{'vcfParameters'}->{'ALT'}->{$i} = $id;
        }
        elsif ( $line =~ /\#\#DLR=(.*)/ ) {
            my $dlr = $1;
            unless ( $dlr =~ /[[:alpha:]]/ ) {
                $fileInfoRef->{'DLR'} = $dlr;
            }
        }
        elsif ( $line =~ /##/ ) {
            push( @{ $fileInfoRef->{'vcfParameters'}->{'OTHER'} }, $line );
        }
        elsif ( $line =~ /^#CHR/ ) {
            my @temp = split( /\t/, $line );
            $fileInfoRef->{'vcfParameters'}->{'chrLine'} = $line;
            for ( my $i = 9 ; $i <= $#temp ; ++$i ) {
                push( @{ $fileInfoRef->{'sampleNames'} }, $temp[$i] );
            }
            next HEADER;
        }
        push( @{ $fileInfoRef->{'vcfParameters'}->{'HEADER'} }, $line );
    }
    $RunParameters->{'filesCollection'}->update( { 'filepath' => $fileInfoRef->{'filepath'} }, $fileInfoRef, {  safe => 1 } );
    return ( $fileInfoRef, $line );
}    #    openVcf($RunParameters,$fileInfoRef,$filepath)=(fileInfoRef,line)

sub CreateOnPath {
    my $function      = shift;
    my $RunParameters = shift;
    print "\t+Create On Path - Removing old records with CreateOnPath=$RunParameters->{'CreateOnPath'}, and Adding new ones\n";
    unless ( exists( $RunParameters->{ $RunParameters->{'uidName'} } ) ) { die "+\t! With CreateOnPath, specify $RunParameters->{'uidName'}\n"; }
    if ( exists( $RunParameters->{'CreateOnPath'} ) ) {
        unless ( $RunParameters->{'CreateOnPath'} =~ /\/$/ ) {
            $RunParameters->{'CreateOnPath'} = $RunParameters->{'CreateOnPath'} . "/";
        }
        $RunParameters->{'vcfPath'} = $RunParameters->{'CreateOnPath'};
        $RunParameters->{'AddOnly'} = $RunParameters->{'CreateOnPath'};
        for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
            my $collectionName=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};   
           $RunParameters->{$collectionName}->remove( { 'CreateOnPath' => $RunParameters->{'CreateOnPath'} }, { safe => 1, 'multiple' => 1 } );
           $RunParameters->{$collectionName}->remove( { 'CreateOnPath' => $RunParameters->{'CreateOnPath'} }, { safe => 1, 'multiple' => 1 } );    
        }            

        #$RunParameters->{'AddOnly'}=$RunParameters->{'CreateOnPath'};
        for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
            my $variantCollectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};   
            my $conn = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($variantCollectionDb);
            my @list = ();
            if (exists ( $RunParameters->{'mtime'})) {
               @list = `find $RunParameters->{'vcfPath'} -iname "*.vcf" -mtime $RunParameters->{'mtime'}`;
            }else{
                @list = `find $RunParameters->{'vcfPath'} -iname "*.vcf"`;
           }
            foreach my $file (@list) {
                $file = File::Spec->rel2abs($file);
                chomp($file);
                $conn->remove( { $RunParameters->{ $RunParameters->{'uidName'} } => $RunParameters->{'uidName'}, 'varArray' => 1 },
                    { safe => 1, 'multiple' => 1 } );
                $conn->update(
                    { 'variants' => { '$elemMatch' => { $RunParameters->{ $RunParameters->{'uidName'} } => $RunParameters->{'uidName'} } } },
                    { '$pull'    => { 'variants'   => { $RunParameters->{ $RunParameters->{'uidName'} } => $RunParameters->{'uidName'} } } },
                    { 'multiple' => 1, 'safe' => 1 }
                );
                $conn->update(
                    { 'variants' => { '$elemMatch' => { $RunParameters->{ $RunParameters->{'uidName'} } => $RunParameters->{'uidName'} } } },
                    { '$inc'     => { 'varArray'   => -1 } },
                    { safe => 1, 'multiple' => 1 }
                );
                $RunParameters->{'filesCollection'}->remove( { 'filepath' => $file }, { safe => 1, 'multiple' => 1 } );
                $RunParameters->{'uidCollection'}
                  ->remove( { $RunParameters->{'uidName'} => $RunParameters->{ $RunParameters->{'uidName'} } }, { safe => 1, 'multiple' => 1 } );
            }
        }
        $RunParameters->{'filesCollection'}
          ->remove( { $RunParameters->{'uidName'} => $RunParameters->{ $RunParameters->{'uidName'} } }, { 'safe' => 1, 'multiple' => 1 } );
        $RunParameters->{'uidCollection'}
          ->remove( { $RunParameters->{'uidName'} => $RunParameters->{ $RunParameters->{'uidName'} } }, { safe => 1, 'multiple' => 1 } );
        $RunParameters->{'filesCollection'}->remove( { 'CreateOnPath' => $RunParameters->{'CreateOnPath'} }, { safe => 1, 'multiple' => 1 } );
        $RunParameters->{'uidCollection'}->remove( { 'vcfPath' => $RunParameters->{'CreateOnPath'} }, { 'multiple' => 1, 'safe' => 1 } );
    }
}

sub DeleteOnPath {
    my $function      = shift;
    my $RunParameters = shift;
    unless ( exists( $RunParameters->{ $RunParameters->{'uidName'} } ) ) { die "+\t! With CreateOnPath, specify $RunParameters->{'uidName'}\n"; }
    if ( exists( $RunParameters->{'DeleteOnPath'} ) ) {
        unless ( $RunParameters->{'DeleteOnPath'} =~ /\/$/ ) {
            $RunParameters->{'DeleteOnPath'} = $RunParameters->{'DeleteOnPath'} . "/";
        }
        print "\t!Deleting on $RunParameters->{'DeleteOnPath'}\n";

    }
}


############## Startup And Misc. ##########################
sub countGenotypes {
    my $variantRef = shift;
    my $het        = int(0);
    my $hom        = int(0);
    my $oth        = int(0);
    my $all        = int(0);
    my $ref        = int(0);
    foreach my $variant ( keys %{ $variantRef->{'genotype'} } ) {
        if ( exists( $variantRef->{'genotype'}->{$variant}->{'gt'} ) ) {
            if ( $variantRef->{'genotype'}->{$variant}->{'gt'} eq '0/1' ) {
                ++$het;
            }
            elsif ( $variantRef->{'genotype'}->{$variant}->{'gt'} eq '1/1' ) {
                ++$hom;
            }
            elsif ( $variantRef->{'genotype'}->{$variant}->{'gt'} eq "0/0" ) {
                ++$ref;
            }
            elsif ( $variantRef->{'genotype'}->{$variant}->{'gt'} ne '0/0' ) {
                ++$oth;
            }
            ++$all;
        }
        else { ++$oth; ++$all }
    }
    return ( $het, $hom, $oth, $all, $ref );
}

sub swapGT {
    my ($variantRef) = @_;
    my @swap_ar     = split( /\|/, $variantRef->{'swap'} );
    my %swap_change = ();
    my %swap_willbe = ();
    my $set         = "";
    foreach $set (@swap_ar) {
        if ( $set =~ /(.*)->(.*)/ ) {
            $swap_change{$1} = $2;                   #gt1->gt2
            $swap_willbe{$1} = $variantRef->{$2};    #val for gt1 is now val for gt2
        }
    }
    foreach $set ( keys %swap_change ) {
        $variantRef->{$set} = $swap_willbe{$set};
    }
    return $variantRef;
}    #     swapGT(variantRef)=variantRef

sub mongoParse {
    my ( $line, $fileInfoRef ) = @_;
    my @temp = split( /\s+/, $line );
    my $fileInfo = $_[1];
  LOOP: foreach my $set (@temp) {
        if ( $set =~ /=/ ) {
            my ( $key, $val ) = split( "=", $set );
            $key = lc($key);
            if ( $key eq "patient" || $key eq "study" ) { next LOOP }
            $fileInfoRef->{$key} = $val;
        }
    }
    return $fileInfoRef;
}    #    mongoParse(line,fileInfoRef)=fileInfoRef

sub addInHashRef {
    my ( $main_hash_ref, $sub_hash_ref ) = @_;
    foreach my $key ( keys %$sub_hash_ref ) {
        unless ( exists( $main_hash_ref->{$key} ) ) {
            $main_hash_ref->{$key} = $sub_hash_ref->{$key};
        }
    }
    return $main_hash_ref;
}    #    addInHashRef (mainHashRef,subHashRef)=mainHashRef

sub printDefaults {
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    print "\n________________________ System Defaults__________________ \n";
    $Data::Dumper::Indent = 0;
    $Data::Dumper::Pair   = " : ";
    my $str = Dumper $RunParameters;
    $str =~ s/\$VAR/\$RunParameters/g;
    $str =~ s/\s//g;
    $str =~ s/\,/, /g;
    print "sysRef: $str\n______________________________________________ \n";
}


############## Filtering ##########################
sub FilterVariantPipeline {
    my $function      = shift;
    my $variantRef    = shift;
    my $RunParameters = shift;
    unless ( exists( $variantRef->{'filetype'} ) ) { return }
    my $fileFilterType    = $variantRef->{'filetype'};
    my $newDatabaseFilter = $variantRef->{'varDatabase'};
    my $newVariantFilter  = $RunParameters->{'InitVariantStatus'};
    if ( exists( $RunParameters->{'FilterVariantPipeline'} ) ) {
      LOOP: for ( my $i = 0 ; $i < scalar( @{ $RunParameters->{'FilterVariantPipeline'} } ) ; ++$i ) {
            my $state = 'none';
            foreach my $message ( sort keys %{ $RunParameters->{'FilterVariantPipeline'}->[$i] } ) {
                if ( $message =~ /SetVariantTo(.*)/ ) {
                    $state = $1;
                    ( $newDatabaseFilter, $newVariantFilter ) = InsertVar->SetVariantTo( $state, $variantRef, $RunParameters, $i, $newDatabaseFilter, $newVariantFilter );
                } else {
                    print "\t\t+Warning, we shouldn't be here.  Did you set up FilterVariantPipeline with SetVariantToPASS?\n";
                }
            }
        }
    }
    $variantRef->{'varDatabase'}=$newDatabaseFilter;
    $variantRef->{'filter'}=$newVariantFilter ;
    
}

sub SetVariantTo {
    my $function          = shift;
    my $state             = shift;
    my $variantRef        = shift;
    my $RunParameters     = shift;
    my $i                 = shift;
    my $newDatabaseFilter = shift;
    my $newVariantFilter  = shift;
    if ( length($state) == 0 ) { return ( $newDatabaseFilter, $newVariantFilter ) }
    my $change = "SetVariantTo" . $state;
    unless ( exists( $variantRef->{'filetype'} ) ) {
        print "\t\t+Warning, I don't think you should ever be here\n";
        return ( $newDatabaseFilter, $newVariantFilter );
    }
    my $fileFilterType = $variantRef->{'filetype'};
    #my $curFilter = $variantRef->{'filter'};
    unless ( exists( $RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{$fileFilterType} ) ) {
        return ( $newDatabaseFilter, $newVariantFilter );
    }
    unless (exists($RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{'DatabaseIfFilter'})) {
        $RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{'DatabaseIfFilter'}="PASS";
    }
    my $dbState = $RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{'DatabaseIfFilter'};    
    for ( my $j = 0 ; $j < scalar( @{ $RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{$fileFilterType} } ) ; ++$j ) {
         my $ref=$RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{$fileFilterType}->[$j];
         my $filter=$RunParameters->{'FilterVariantPipeline'}->[$i]->{$change}->{$fileFilterType}->[$j]->{'field'};
         if ( exists( $variantRef->{$filter} )) { 
              if ( $newVariantFilter ne $state ) { #success,fail
                  $newVariantFilter = Maintain->checkConditions( $variantRef, $filter, $ref, $state, $newVariantFilter );
              }
         }  elsif ($ref->{'oper'} eq 'UNDEF') {
              $newVariantFilter=$state;
         }
         if ( $newVariantFilter eq $state ) {
             $newVariantFilter  = $state;
             $newDatabaseFilter = $dbState;
         } else {
             $newVariantFilter  = $newVariantFilter;
             $newDatabaseFilter = $newDatabaseFilter;
         }
    }
    return ( $newDatabaseFilter, $newVariantFilter );
}

sub excludeNonReportableVariants {
    my $RunParameters = shift;
    unless ( exists( $RunParameters->{'nonReportableVariants'} ) ) {
        return 0;
    }
    print "..Excluding variants in $RunParameters->{'nonReportableVariants'}.";
    my $count = 0;
    open( FILE, $RunParameters->{'nonReportableVariants'} ) or die "Can't open file in openVC: $RunParameters->{'nonReportableVariants'}!\n";
  LOOP: while (<FILE>) {
        my $line = $_;
        chomp($line);
        if ( $line =~ /^#/ ) { next LOOP }
        my @vcfFields = split( /\t+/, $line );
        if ( $#vcfFields < 4 ) { next LOOP }
        $vcfFields[0] =~ s/chr//g;
        my $coord = "chr" . $vcfFields[0] . ":" . $vcfFields[1] . ":" . $vcfFields[3] . ":" . $vcfFields[4];
        $excludedVariants->{$coord} = 1;
        ++$count;
    }
    close(FILE);
    print "\t\t+Found $count variants to exclude.\n";
}
1;
