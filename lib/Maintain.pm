#!/usr/bin/perl -w
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
# * Author(s):
#    David Craig
#    Release 9/1/15
#/
######################################################################################
package Maintain;
use strict;
use MongoDB;
use MongoDB::OID;
use Time::localtime;
use Data::Dumper;
use Storable;
use Sys::Hostname;
use Scalar::Util qw(looks_like_number);
use Digest::MD5 qw(md5 md5_hex md5_base64);


$| = 1;

sub compactCollection {
    my $function      = shift;
    my $RunParameters = shift;
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        my $holder = $RunParameters->{'markersDb'}->run_command( { compact => $collectionDb, paddingBytes => int(3) } );
    }
    return 1;
}

sub validStartCheck {
    my $function      = shift;
    my $RunParameters = shift;
    print "\t\t+Checking parameters file and validating RunParameters Object\n";
    if ( exists( $RunParameters->{'vcfPath'} ) ) {
        unless ( $RunParameters->{'vcfPath'} =~ /\/$/ ) {
            $RunParameters->{'vcfPath'} = $RunParameters->{'vcfPath'} . "/";
        }
    }
    unless ( exists( $RunParameters->{'varGlobalMax'} ) ) {
        $RunParameters->{'varGlobalMax'} = int(10000);
    }
    if (exists( $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} }) &&
        exists( $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'CreateOnPath'} })) {
        $RunParameters->{'CreateOnUid'} = $RunParameters->{'StartupRunParameters'}->{ $RunParameters->{'uidName'} };
    }
    if ( exists( $RunParameters->{'CreateOnUid'} ) || exists( $RunParameters->{'CreateOnPath'} ) ) {
        unless ( $RunParameters->{'CreateOnPath'} =~ /\/$/ ) {
            $RunParameters->{'CreateOnPath'} = $RunParameters->{'CreateOnPath'} . "/";
        }
        unless ( exists( $RunParameters->{'CreateOnUid'} ) || exists( $RunParameters->{'CreateOnPath'} ) ) {
            die "\t!!!! Please specify both CreateOnUid/$RunParameters->{'uidName'} with CreateOnPath on commandline";
        }
        $RunParameters->{'vcfPath'} = $RunParameters->{'CreateOnPath'};
    }
    unless ( exists ($RunParameters->{'Collections'})) {
       die "\t!!!!! No Collections Array Specified\n";
    }
    for (my $i=0; $i< scalar(@{$RunParameters->{'Collections'}});++$i) {
       if (exists($RunParameters->{'Collections'}->[$i]->{'CollectionName'}) &&
           exists($RunParameters->{'Collections'}->[$i]->{'GroupField1'}) &&
           exists($RunParameters->{'Collections'}->[$i]->{'GroupField2'})
          ) {
              print "\t\t+ Valid Collection: $RunParameters->{'Collections'}->[$i]->{'CollectionName'}.\n";
              $RunParameters->{$RunParameters->{'Collections'}->[$i]->{'CollectionName'} . '-index'}=$i;
              $RunParameters->{'CollectionCount'}=$i;
          } else {
              die "\t  !Missing parameter of either CollectionName, GroupField1, GroupField2 in Collections Run Parameter\n";
          }
    }
    $RunParameters->{loadLocations} = 1;
    foreach my $parameter (
                              'uidName','MongodbHost','AnnotateHost','MongodbHost',
                              'DatabaseName','Snpeff4Path','SnpeffParameter',
                              ,'maxInsertsPerRun', 'gemGeneLocPath', 'reportableGenes',
                              'runDir', 'StartupRunParameters'
                            ) {
        unless ( exists( $RunParameters->{$parameter} ) ) {
            die "!!!! MISSING MANDATORY PARAMETER ($parameter !!!!!!!)\n";
        }
    }
}

sub writeVCF {
    my $function      = shift;
    my $RunParameters = shift;
    my $filepath = shift;
    my $dir      = shift;
    my $fileConn = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection('files');
    my $cursor   = $fileConn->find( { 'filepath' => $filepath } );
    while ( my $doc = $cursor->next ) {
        my $patient        = $doc->{'patient'};
        my $study          = $doc->{'study'};
        my $assay          = $doc->{'assay'};
        my $tissue         = $doc->{'tissue'};
        my $collectionDb   = $doc->{'collectionDb'};
        my $collectionDate = $doc->{'collectionDate'};
        my $projectRun     = $doc->{ $RunParameters->{'uidName'} };
        my $name           = "$doc->{'projectRun'}.$doc->{'study'}";
        my @sampleNames    = ();
        my $numNames       = 0;

        foreach my $sampleName ( @{ $doc->{'sampleNames'} } ) {
            push( @sampleNames, $sampleName );
            ++$numNames;
        }
        if ( exists( $doc->{'patient'} ) )      { $name .= "_$doc->{'patient'}" }
        if ( exists( $doc->{'assay'} ) )        { $name .= "_$doc->{'assay'}" }
        if ( exists( $doc->{'tissue'} ) )       { $name .= "_$doc->{'tissue'}" }
        if ( exists( $doc->{'collectionDb'} ) ) { $name .= "_$doc->{'collectionDb'}" }
        $name .= ".$doc->{'filetype'}.GCF";
        unless ( -e $dir ) { system("mkdir -p $dir") }
        open( VCF, ">$dir/$name" ) or die "Can't open final VCF file\n";
        print VCF "##fileformat=VCFv4.1\n";
        print VCF "##source=GemDb\n";

        foreach my $field ( keys %{ $doc->{'vcfParameters'}->{'INFO'} } ) {
            print VCF "##INFO=<" . $field;
            foreach my $property ( keys %{ $doc->{'vcfParameters'}->{'INFO'}->{$field} } ) {
                print VCF ",$property=$doc->{'vcfParameters'}->{'INFO'}->{$field}->{$property}";
            }
            print VCF ">\n";
        }
        foreach my $field ( keys %{ $doc->{'vcfParameters'}->{"FORMAT"} } ) {
            print VCF "##FORMAT=<" . $field;
            foreach my $property ( keys %{ $doc->{'vcfParameters'}->{'FORMAT'}->{$field} } ) {
                print VCF ",$property=$doc->{'vcfParameters'}->{'FORMAT'}->{$field}->{$property}";
            }
            print VCF ">\n";
        }
        foreach my $line ( @{ $doc->{'vcfParameters'}->{'OTHER'} } ) {
            if ( $line =~ /SnpEffVersion/ ) {
                $line = "SnpEffVersion=4.0\n";
            }
            elsif ( $line =~ /SnpEffCmd/ ) {
                $line = "SnpEffCmd=$RunParameters->{'Snpeff4Path'} $RunParameters->{'SnpeffParameter'}\n";
            }
            print VCF "$line\n";
        }
        print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" . join( "\t", @sampleNames ) . "\n";
        my $conn = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
        my $cursor1;
        my $cursor2 = $conn->aggregate( [ { '$match' => { "variants.projectRun" => $projectRun } }, { '$unwind' => '$variants' }, { '$match' => { "variants.filepath" => $filepath } }, { '$limit' => 50 } ], { 'cursor' => { batchSize => 0 } }, $cursor1 );
        while ( my $vDoc = $cursor2->next ) {
            my $infoField    = "GemDb=$vDoc->{'info'}->{'version'}";
            my $gtFieldsText = "";
            my $gtValsText   = "";
            foreach my $infoFieldKey ( 'biomarker', 'gene', 'coord', 'study', 'effect', 'collectionDb', 'group', 'vcfPath', 'type', 'molecule', 'dbFreqCat', 'dbFreq', 'biomarkerCount', 'drug_rule_matched' ) {
                if ( exists( $vDoc->{$infoFieldKey} ) ) {
                    if ( defined( $vDoc->{$infoFieldKey} ) ) {
                        $infoField .= ";$infoFieldKey=$vDoc->{$infoFieldKey}";
                    }
                }
            }
            foreach my $variantKey ( keys %{ $vDoc->{'variants'} } ) {
                if ( defined($variantKey) && $variantKey ne 'genotype' ) {
                    $infoField .= ";$variantKey=$vDoc->{'variants'}->{$variantKey}";
                }
                elsif ( $variantKey eq "genotype" ) {
                    my @gtFields = ();
                    my @gtVals   = ();
                    foreach my $gtNum ( sort keys %{ $vDoc->{'variants'}->{'genotype'} } ) {
                        foreach my $gtKey ( sort keys %{ $vDoc->{'variants'}->{'genotype'}->{$gtNum} } ) {
                            if ( $gtKey ne "sampleName" ) {
                                push( @gtFields, $gtKey );
                                push( @gtVals,   $vDoc->{'variants'}->{'genotype'}->{$gtNum}->{$gtKey} );
                            }
                        }
                    }
                    if ( scalar( @sampleNames > 0 ) ) {
                        $gtFieldsText = join( ":", @gtFields );
                        $gtValsText .= "\t" . join( ":", @gtVals );
                    }
                }
            }
            foreach my $annotateKey ( 'dbsnp', 'frequency', 'frequencyCategory', 'phyloP', '1000gAmrMaf', 'fathmm', '1000gAsnMaf' ) {
                if ( defined($annotateKey) && exists( $vDoc->{'annotate'}->{$annotateKey} ) ) {
                    $infoField .= ";annotate_$annotateKey=$vDoc->{'annotate'}->{$annotateKey}";
                }
            }
            foreach my $snpeffInfoKey ( 'hgsv', 'AminoAcidChange', 'functionalImpact', 'functionalClass', 'transcriptId', 'functionalEffect' ) {
                if ( defined($snpeffInfoKey) && exists( $vDoc->{'snpeff'}->{$snpeffInfoKey} ) ) {
                    $infoField .= ";snpeff_$snpeffInfoKey=$vDoc->{'snpeff'}->{$snpeffInfoKey}";
                }
            }
            foreach my $aberrationInfoKey ( keys %{ $vDoc->{'aberration'} } ) {
                if ( defined($aberrationInfoKey) && exists( $vDoc->{'aberration'}->{$aberrationInfoKey} ) ) {
                    $infoField .= ";aberration_$aberrationInfoKey=$vDoc->{'aberration'}->{$aberrationInfoKey}";
                }
            }
            print VCF join( "\t",
                $vDoc->{'variants'}->{'chr'},  $vDoc->{'variants'}->{'hg19pos'},             $vDoc->{'variants'}->{'varName'}, $vDoc->{'variants'}->{'ref'}, $vDoc->{'variants'}->{'alt'},
                $vDoc->{'variants'}->{'qual'}, $vDoc->{'aberration'}->{'aberration_filter'}, $infoField,                       $gtFieldsText,                $gtValsText )
              . "\n";
        }
        close(VCF);
    }
    print "\t\t+Finished Writing VCF for $filepath\n";

}

sub replaceInHash {
    my ( $main_hash_ref, $sub_hash_ref ) = @_;
    foreach my $key ( keys %$sub_hash_ref ) {
        $main_hash_ref->{$key} = $sub_hash_ref->{$key};
    }
    return $main_hash_ref;
}    #    replaceInHashRef (mainHashRef,subHashRef)=mainHashRef

sub assignGroup {
    my $function      = shift;
    my $biomarkerRef  = shift;
    my $RunParameters = shift;
    my $collectionIndex=$RunParameters->{$biomarkerRef->{'collectionDb'} . '-index'};
    if (exists($biomarkerRef->{'variants'}->[0]->{$RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'}})) {
        $biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'}}=
           $biomarkerRef->{'variants'}->[0]->{$RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'}};
    } else {
        $biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'}}='none';
    }
    if (exists($biomarkerRef->{'variants'}->[0]->{$RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'}})) {
        $biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'}}=
           $biomarkerRef->{'variants'}->[0]->{$RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'}};
    } else {
        $biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'}}='none';
    }
    $biomarkerRef->{'info'}->{'groupField'} = $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'} .
                                    ":" . $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'} .
                                    ":" . 'vcfPath';
    $biomarkerRef->{'info'}->{'group'}      = md5_hex($biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField1'}  } .
                                    ":"  . $biomarkerRef->{ $RunParameters->{'Collections'}->[$collectionIndex]->{'GroupField2'} } .
                                    ":" . $biomarkerRef->{'vcfPath'});
}

sub calculateDbFreq {
    my $function      = shift;
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    print "--------------------- Counting biomarkers ----\n";
    print "\t+Counting inserts in new projects\n";
    my $start      = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $cursor     = $RunParameters->{'uidCollection'}->find( { 'inserted' => int(-1) }, { safe => 1 } );
    my $newInserts = 0;
    $cursor->immortal(1);
    ## Count per projectRun and files
    while ( my $doc = $cursor->next ) {
        $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $doc->{ $RunParameters->{'uidName'} } },{'$set'=>{'inserted'=>int(0)}},{safe=>0});
               $RunParameters->{'germline'}->ensure_index( { 'variants.filepath' => 1 } );
               $RunParameters->{'tumor'}->ensure_index( { 'variants.filepath' => 1 } );
               $RunParameters->{'tumor'}->ensure_index( { 'variants.projectRun' => 1 } );
               $RunParameters->{'germline'}->ensure_index( { 'variants.projectRun' => 1 } );
               $RunParameters->{'tumor'}->ensure_index( { 'aberration.aberration_filter' => 1 } );
               $RunParameters->{'germline'}->ensure_index( { 'aberration.aberration_filter' => 1 } );
               $RunParameters->{'tumor'}->ensure_index( Tie::IxHash->new('variants.projectRun' => 1, 'aberration.aberration_filter' => 1) );
               $RunParameters->{'germline'}->ensure_index( Tie::IxHash->new('variants.projectRun' => 1, 'aberration.aberration_filter' => 1) );
               $RunParameters->{'tumor'}->ensure_index( { 'type' => 1 } );
               $RunParameters->{'germline'}->ensure_index( { 'type' => 1 } );
                my $fileInsertCount = $RunParameters->{ 'tumor' }->count(  Tie::IxHash->new('variants.projectRun' =>  $doc->{$RunParameters->{'uidName'} },'aberration.aberration_filter'=>'PASS') )
                                       +$RunParameters->{ 'germline' }->count( Tie::IxHash->new( 'variants.projectRun' =>  $doc->{$RunParameters->{'uidName'} },'aberration.aberration_filter'=>'PASS' ));
                my $countSNV=  $RunParameters->{ 'tumor' }->count( { 'type'=>'SNV','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countINS=  $RunParameters->{ 'tumor' }->count( { 'type'=>'smallInsertion','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countDEL=  $RunParameters->{ 'tumor' }->count( { 'type'=>'smallDeletion','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countSV=  $RunParameters->{ 'tumor' }->count( { 'type'=>'StructuralVariant','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countFocalGain=  $RunParameters->{ 'tumor' }->count( { 'type'=>'FocalGain','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countFocalDeletion=  $RunParameters->{ 'tumor' }->count( { 'type'=>'FocalDeletion','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countFusions=  $RunParameters->{ 'tumor' }->count( { 'type'=>'Fusions','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countSNPs=  $RunParameters->{ 'germline' }->count( { 'type'=>'SNV','variants.projectRun' =>  $doc->{$RunParameters->{'uidName'}},'aberration.aberration_filter'=>'PASS'});
                my $countPointMutations=$countSNV+$countINS+$countDEL;
                my $countCNA=$countFocalGain+$countFocalDeletion;
                $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $doc->{ $RunParameters->{'uidName'} } },
                                                           { '$set'=>{'countSNV'=>$countSNV,'countINS'=>$countINS,'countDel'=>$countDEL,'countSV'=>$countSV,
                                                                      'countFocalDeletion'=>$countFocalDeletion,'countPointMutations'=>$countPointMutations,'inserted'=>$fileInsertCount,
                                                                      'countCNA'=>$countCNA ,'countSNPs'=>$countSNPs,'countFusions'=>$countFusions}},
                                                           { safe => 1 ,'multiple'=>1}
                                                         );
                    my $cursor3=$RunParameters->{ 'tumor' }->find({'variants.projectRun' =>$doc->{ $RunParameters->{'uidName'}},'type'=>'SNV'});
                    my @array=();
                    while (my $doc3=$cursor3->next) {
                        if (exists($doc->{'variants'}->[0]->{'AR2'})) {
                            push(@array,$doc3->{'variants'}->[0]->{'AR2'});
                        } elsif (exists($doc3->{'variants'}->[0]->{'SEURAT_AR_TUMOR'})) {
                            push(@array,$doc3->{'variants'}->[0]->{'SEURAT_AR_TUMOR'});
                        } elsif (exists($doc3->{'variants'}->[0]->{'DNA_ALT_ALLELE_TOTAL_FRACTION'})) {
                            push(@array,$doc3->{'variants'}->[0]->{'DNA_ALT_ALLELE_TOTAL_FRACTION'});
                        }
                    }
                    if ($#array>0) {
                          my (@data) = sort { $a <=> $b } @array;
                          my $med=0;
                           if ( scalar(@data) % 2 ) {
                               $med=$data[ @data / 2 ]
                           } else {
                               my ( $upper, $lower );
                               $lower = $data[ @data / 2 ];
                               $upper = $data[ @data / 2 - 1 ];
                               $med=($upper+$lower)/2;
                           }
                       MongoDB::force_double($med);
                       $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $doc->{ $RunParameters->{'uidName'} } },
                                                                  { '$set' => { 'medianMaf' => $med } },
                                                                  { safe => 1 }
                                                                );
                    }   
                foreach my $filepath ( @{ $doc->{'filepaths'} } ) {
                    my $cursor2 = $RunParameters->{'filesCollection'}->find( { 'filepath' => $filepath }, { safe => 1 } );
                    LOOP: while ( my $doc2 = $cursor2->next ) {
                        unless (exists($doc2->{'collectionDb'})) {next LOOP}
                       $RunParameters->{$doc2->{'collectionDb'}}->ensure_index( { 'variants.filepath' => 1 } );
                        my $fileInsertCount = $RunParameters->{ $doc2->{'collectionDb'} }->count( { 'variants.filepath' => $filepath } );
                        $RunParameters->{'filesCollection'}->update( { '_id' => $doc2->{'_id'} }, { '$set' => { 'inserted' => $fileInsertCount } }, { safe => 1 ,multi=>1} );
                        if (exists($doc2->{'filename'})) {
                            print "\t\t\t+Updated projectRun/$RunParameters->{'uidName'}: $doc2->{$RunParameters->{'uidName'}} filename:$doc2->{'filename'} $fileInsertCount inserted\n";
                        } else {print "\t\t\tmissing projectRun/$RunParameters->{'uidName'}\n";}
                    }
                }
            $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $doc->{ $RunParameters->{'uidName'} } },
                                                        { '$set' => { 'status' => int(1) ,'reportable'=>int(1)} },
                                                        { safe => 1 }
                                                      );
    }
    print "\t\t+Total of $newInserts\n";
    print "\t+Counting inserts going through files\n";
    $newInserts = 0;
    $cursor = $RunParameters->{'filesCollection'}->find( { 'inserted' => int(-1),'status'=>int(1) }, { safe => 1 } );
    LOOP: while ( my $doc = $cursor->next ) {
        my $filepath = $doc->{'filepath'};
       unless (exists($doc->{'collectionDb'})) {next LOOP}
        my $fileInsertCount = $RunParameters->{ $doc->{'collectionDb'} }->count( { 'variants.filepath' => $filepath } );
        $RunParameters->{'filesCollection'}->update( { '_id' => $doc->{'_id'} }, { '$set' => { 'inserted' => $fileInsertCount } }, { safe => 1 ,multi=>1} );
        $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $doc->{ $RunParameters->{'uidName'} } },
                                                   { '$inc' => { 'inserted' => $fileInsertCount } },
                                                   { safe => 1 }
                                                 );
        print "\t\t\t+Updated files/$RunParameters->{'uidName'}, $doc->{'filename'} $fileInsertCount inserted\n";
        $newInserts += $fileInsertCount;
    }
    print "\t\t+Total of $newInserts\n";
    print "\t+Run dbFreq $function...\n";
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $vcfCollection=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        my $fileConn = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection('files');
        $RunParameters->{$vcfCollection}->ensure_index( { 'variants.filepath' => 1 } );
#        $RunParameters->{$vcfCollection}->ensure_index( { 'variants.coord' => 1 } );
        $RunParameters->{$vcfCollection}->ensure_index( { 'biomarker' => 1 } );
        my $cursor = $fileConn->find( { collectionDb => $vcfCollection } );
        my $patientNum = 0;
        print "\t+Counting number of samples\n";
      PATIENT: while ( my $doc = $cursor->next ) {    #last PATIENT;
            if ( exists( $doc->{'sampleNames'} ) ) {
                my $val = scalar( @{ $doc->{'sampleNames'} } );
                $patientNum = $patientNum + $val;
            }
            else { ++$patientNum }
        }
        print "\t+Counting Markers in $vcfCollection\n";
        my $markerConn = $RunParameters->{$vcfCollection};
        $markerConn->ensure_index({'info.checkCount'=>1});
        if (exists($RunParameters->{'resetCount'})) {
            if ($RunParameters->{'resetCount'}==1) {
                $markerConn->update( { }, { '$set' => { 'info.checkCount' => int(1) } }, { 'safe' => 1,'multiple'=>1 } );
            }
        }
        my %counted = ();
        my $superCount=0;
       RESET:
        $cursor = $markerConn->find({'info.checkCount'=>1});
        $cursor->immortal(1);
        my $counter = 0;
        my $totalCount=0;
        my $totalCountPrev=$totalCount;
        my $reportStep=10000;
        my $prev=localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
      REQ: while ( my $doc = $cursor->next ) {
            my $biomarkerAlleleCount = 0;
            my $biomarkerCount       = 0;
            ++$counter;
            ++$superCount;
            if ( $counter % $reportStep == 0 ) {
                 my $diff = sprintf("%2.1f",($totalCount-$totalCountPrev)/(localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $prev+0.001));
                 print "\t\t\t+Db Freq Count is at $counter, $totalCount with $diff per s.\n";
                 $prev=localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
                 $totalCountPrev=$totalCount;
            }
            if ( exists( $counted{ $doc->{'biomarker'} } ) ) { next REQ }
            my $biomarker = $doc->{'biomarker'};
            my $out=$markerConn->aggregate([{'$match'=>{'biomarker' => $biomarker }},{'$group'=>{'_id'=>'biomarker','count'=>{'$sum'=>'$info.varArray'}}}]);
            $biomarkerAlleleCount=int($out->[0]->{'count'});
            $totalCount+=$biomarkerAlleleCount;
            my $thresh = 0.03;   
            my $cat    = 'unknown';                     
            if ( exists( $RunParameters->{'rareFreqThresh'} ) ) { $thresh = $RunParameters->{'rareFreqThresh'}; }                        
            my $freq = $thresh;            
            if (exists( $doc->{'annotate'}->{'maxPopFreq'})) {
#                if ($doc->{'annotate'}->{'maxPopFreq'} > $freq) {
                    $freq=$doc->{'annotate'}->{'maxPopFreq'};
#                }elsif ($freq>$doc->{'annotate'}->{'maxPopFreq'}) {
            }            
            if ( $patientNum > 20 ) {
                $freq = sprintf( "%0.5f", $biomarkerAlleleCount / $patientNum );
            } else {
                $freq =$thresh;
            }
            MongoDB::force_double($freq);
            if ( $freq <= $thresh ) { $cat = "rare";} else { $cat = "common"; }
            $markerConn->update( { 'biomarker' => $biomarker },
                                 { '$set' => { 'dbFreq' => $freq,
                                               'dbFreqCat' => $cat,
                                               'biomarkerCount' => int($biomarkerAlleleCount),
                                               'info.checkCount'=>int(0)
                                            }
                                 },
                                 { safe => 'true', multiple => 1 }
                               );
            $counted{$biomarker} = 1;
            if ( $totalCount > 1000000 && $counter> 10*$reportStep ) {
                        my $remainingToCount = $markerConn->count({'info.checkCount'=>1});
                        print "\t\t+Reset Count (total counted=$superCount, remaining $remainingToCount)\n";
                        goto RESET
            }
        }
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "--------------------- Done Counting (Db Freq Completed in $diff seconds) ----\n";
    1;
}

sub cleanOld {
    my $RunParameters = shift;
    print "\n------Checking and Scrubbing Old Files----\n";
    my $cursor3 = $RunParameters->{'markersDb'}->get_collection('files')->find( { 'status' => 3 } );
    while ( my $doc = $cursor3->next ) {
        print "\t\t- Found status 3: $doc->{'filepath'}\n";
        if ( exists( $doc->{'collectionDb'} ) ) {
            InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters );
        }
        else {
            print "Could not found  collectionDb in $doc->{'filepath'} \n";
        }
    }
    my $cursor2 = $RunParameters->{'markersDb'}->get_collection('files')->find( { 'status' => 2 } );
    while ( my $doc = $cursor2->next ) {
        print "\t\t- Found status 2: $doc->{'filepath'}\n";
        if ( exists( $doc->{'collectionDb'} ) ) {
            InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters );
        }
        else {
            print "Could not found  collectionDb in $doc->{'filepath'} \n";
        }
    }
    my $cursor0 = $RunParameters->{'markersDb'}->get_collection('files')->find( { 'status' => 0 } );
    $cursor0->immortal(1);
    while ( my $doc = $cursor0->next ) {
        print "\t\t- Found status 0: $doc->{'filepath'}\n";
        if ( exists( $doc->{'collectionDb'} ) ) {
            InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters );
        }
        else {
            $RunParameters->{'filesCollection'}->remove( { 'filepath' => $doc->{'filepath'} }, { safe => 1, 'multiple' => 1 } );
        }
    }
    my $cursor4 = $RunParameters->{'markersDb'}->get_collection('files')->find( { 'status' => 4 } );
    $cursor4->immortal(1);
    while ( my $doc = $cursor4->next ) {
        print "\t\t- Found status 4: $doc->{'filepath'}\n";
        InsertVar->scrubFilepath( $doc->{'filepath'}, $doc->{'collectionDb'}, $RunParameters );
    }
}

sub purgePartials {
    my $function      = shift;
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my $timeOut = 360000000;
    my $start   = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    if (exists($RunParameters->{'scrubEnabled'} ) ) {
        if ($RunParameters->{'scrubEnabled'} eq "on" || $RunParameters->{'scrubEnabled'} eq "1" || $RunParameters->{'scrubEnabled'} eq "true") {
            cleanOld($RunParameters);
        }
    }
    for (my $i=1;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        print "\n--------------Purge Partials ($collectionDb)-------------- \n";
        my $connMarkers = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
        if (exists($RunParameters->{'largeRun'})) {
            if ($RunParameters->{'largeRun'}==1) {
                #print "\t+Reindexing $collectionDb\n";
                #$connMarkers->drop_indexes;
            }
        }
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
  #  if ( glob("$RunParameters->{'runDir'}/snpeff.*.vcf") ) { `rm $RunParameters->{'runDir'}/snpeff.*.vcf`; }
  #  if ( glob("$RunParameters->{'runDir'}/err_file*") )    { `rm $RunParameters->{'runDir'}/err_file*`; }
    print "---- Purge Completed in $diff seconds-------- \n";
}

sub cleanAndMerge {
    my $function      = shift;
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my $timeOut = 360000000;
    my $start   = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $end=1;
    if (exists($RunParameters->{'scrubEnabled'} ) ) {
        if ($RunParameters->{'scrubEnabled'} eq "on" || $RunParameters->{'scrubEnabled'} eq "1" || $RunParameters->{'scrubEnabled'} eq "true") {
            cleanOld($RunParameters);
        }
    }
    for (my $i=1;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        print "\n--------------Cleaning and Merging ($collectionDb)-------------- \n";
        my $connMarkers = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
        if (exists($RunParameters->{'largeRun'})) {
            if ($RunParameters->{'largeRun'}==1) {
                #print "\t+Reindexing $collectionDb\n";
                #$connMarkers->drop_indexes;
            }
        }
        $connMarkers->ensure_index( { 'info.group' => 1 } );
        $connMarkers->ensure_index( { 'info.compact' => 1 } );
        $connMarkers->ensure_index(Tie::IxHash->new(
                                    'info.group' => 1,
                                    'info.toDelete'=>1,
                                    'info.compact'  => 1
                                   ));
        $connMarkers->ensure_index(Tie::IxHash->new(
                                    'info.group' => 1,
                                    'info.toDelete'=>1,
                                    'info.varFilled'  => 1
                                   ));
        print "\n\n+Starting merging for $collectionDb\n";
        my %foundGroup  = ();
        if (exists($RunParameters->{'resetCompact'})) {
            if ($RunParameters->{'resetCompact'}==1) {
                $connMarkers->update( { }, { '$set' => { 'info.compact' => int(0) } }, { 'safe' => 1,'multiple'=>1 } );
            }
        }
        my $cursor = $connMarkers->find({'info.compact'=>0});
        $cursor->immortal(1);
        my $count  = 1;
        my $totalCount=1;
        my $prev   = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
        my $prev2=1;
        my $firstprev   = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
      OUTERLOOP: while ( my $doc = $cursor->next ) {
            my $id = $doc->{'_id'};
            if ( exists( $foundGroup{ $doc->{'info'}->{'group'} } ) ) { next OUTERLOOP }
            unless ($doc->{'info'}->{'toDelete'} == 0) {next OUTERLOOP}
            if ($doc->{'info'}->{'compact'} == 1) {next OUTERLOOP}
            my $interval=1000;
            my $temp=localtime->hour * 3600 + localtime->min() * 60 + localtime->sec()-$firstprev;
            if ($temp<0) { $firstprev=$firstprev-24*3600;}
            if ($temp>8*3600) {print "\t\t+Ending merging early ($temp - $prev - $firstprev) $count \n"; $end=2; last OUTERLOOP;}
            if ( $count % $interval == 0 ) {
                my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $prev;
                my $diffperinterval=sprintf("%3.2f",$interval/($diff+0.001));
                my $diff2perinterval=sprintf("%3.2f",($totalCount-$prev2)/($diff+0.0001));
                print "\t\t+Merging Count is at $count in " . $diff . "s or ($diffperinterval per s)  with total of $totalCount marked ($diff2perinterval per s), over past $temp seconds.\n";
                $prev = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
                $prev2=$totalCount;
            }
            my $holder = $RunParameters->{'markersDb'}->run_command(
                {
                    findAndModify => $collectionDb,
                    query         => { '_id' => $id },
                    update        => { '$set' => { 'info.toDelete' => $$ } },
                    'new'         => 1
                }
            );
            unless ( $holder->{'value'}->{'info'}->{'toDelete'} eq $$ ) {
                next OUTERLOOP;
            }
            ++$count;
            ++$totalCount;
            my $idx=Tie::IxHash->new(
                                        'info.group' => $doc->{'info'}->{'group'},
                                        'info.toDelete' => int(0),
                                        'info.varFilled'  => int(0)
                                       );
            my $out=$connMarkers->aggregate([{
                                                '$match'=>$idx
                                            },{'$limit'=>$RunParameters->{'varGlobalMax'}+3}]);
            if (scalar(@{$out}) <1) {
                $connMarkers->update( { '_id' => $doc->{'_id'} }, { '$set' => { 'info.compact' => int(1),'info.toDelete'=>int(0) } }, { 'safe' => 1 } );
                next OUTERLOOP;
            }
            my $biomarkerRef=$doc;
        #    $connMarkers->update( { '_id' => $doc->{'_id'} }, { '$set' => { 'info.toDelete' => int($$) } }, { 'safe' => 1 } );
            my $varArrayCount= scalar( @{ $biomarkerRef->{'variants'} });
            my $moreFound=0;
            my $elements=scalar(@{$out});
            MIDLOOP: for ($i=0;$i<$elements;++$i) {
                my $docFirst=$out->[$i];
                unless ($docFirst->{'info'}->{'toDelete'} ==0) {next MIDLOOP}
                $biomarkerRef->{'info'}->{'varMax'}    = $RunParameters->{'varGlobalMax'} ;
                $biomarkerRef->{'info'}->{'varFilled'} = int(0);
                 my $holder = $RunParameters->{'markersDb'}->run_command(
                     {
                         findAndModify => $collectionDb,
                         query         => { '_id' => $docFirst->{'_id'} },
                         update        => { '$set' => { 'info.toDelete' => $$ } },
                         'new'         => 1
                     }
                 );
                 unless ( $holder->{'value'}->{'info'}->{'toDelete'} eq $$ ) {
                     next MIDLOOP;
                 }
#                $connMarkers->update( { '_id' => $docFirst->{'_id'} }, { '$set' => { 'info.toDelete' => int($$) } }, { 'safe' => 1 } );
                foreach my $variantRef ( @{ $docFirst->{'variants'} } ) {
                    ++$totalCount;
                    ++$varArrayCount;
                    push( @{ $biomarkerRef->{'variants'} }, $variantRef );
                    push( @{ $biomarkerRef->{$RunParameters->{'UidName'} } }, $variantRef->{ $RunParameters->{'uidName'} }  );
                }
                $biomarkerRef->{'count'}->{'all'}= $biomarkerRef->{'count'}->{'all'} + $docFirst->{'count'}->{'all'};
                $biomarkerRef->{'count'}->{'hets'}= $biomarkerRef->{'count'}->{'hets'} + $docFirst->{'count'}->{'hets'};
                $biomarkerRef->{'count'}->{'homs'}= $biomarkerRef->{'count'}->{'homs'} + $docFirst->{'count'}->{'homs'};
                $biomarkerRef->{'count'}->{'oths'}= $biomarkerRef->{'count'}->{'oths'} + $docFirst->{'count'}->{'oths'};
                $biomarkerRef->{'count'}->{'refs'}= $biomarkerRef->{'count'}->{'refs'} + $docFirst->{'count'}->{'refs'};
                if ($varArrayCount >= $RunParameters->{'varGlobalMax'} ) {
                    $biomarkerRef->{'info'}->{'varFilled'} = int(1);
                    $moreFound=1;
                    last MIDLOOP;
                }
            }
            $biomarkerRef->{'info'}->{'varArray'}  = scalar( @{ $biomarkerRef->{'variants'} } );
            $biomarkerRef->{'info'}->{'compact'}   = int(1);
            $biomarkerRef->{'info'}->{'toDelete'}   = int(0);
            $biomarkerRef->{'info'}->{'partialInsert'}   = int($$);
            if ( $biomarkerRef->{'info'}->{'varArray'} >= $RunParameters->{'varGlobalMax'} ) {
                $biomarkerRef->{'info'}->{'varMax'}    = $RunParameters->{'varGlobalMax'};
                $biomarkerRef->{'info'}->{'varFilled'} = int(1);
                $biomarkerRef->{'info'}->{'compact'} = int(1);
            }
            if ($moreFound==0) {
                $foundGroup{ $biomarkerRef->{'info'}->{'group'} } = 1;
            }
            my $oldId=$biomarkerRef->{'_id'};
            undef $biomarkerRef->{'_id'};
            delete $biomarkerRef->{'_id'};
            my $newId = $connMarkers->insert( $biomarkerRef, { safe => 1 } );
            $connMarkers->remove( { 'info.toDelete' => int($$) }, { 'safe' => 1 });
            $connMarkers->update( { '_id' => $newId }, { '$unset' => { 'info.partialInsert' => "" } }, { 'safe' => 1 } );
        }
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "---- Cleaning and Merging Completed in $diff seconds-------- \n";
    return $end;
}

sub mongoConnect {
    my $function      = shift;
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my $timeOut = 360000000;
    $MongoDB::Cursor::timeout = $timeOut;
    $RunParameters->{'mongoConnection'} = MongoDB::MongoClient->new( host => "mongodb://$RunParameters->{'MongodbHost'}",
                                                                    query_timeout => $timeOut, timeout => $timeOut, wtimeout => $timeOut
                                                                   ) or die "Could not connect to $RunParameters->{'MongodbHost'}";
    $RunParameters->{'annotateByConn'} = MongoDB::MongoClient->new( host => "mongodb://$RunParameters->{'AnnotateHost'}",
                                                                    query_timeout => $timeOut, timeout => $timeOut, wtimeout => $timeOut
                                                                  )
                                                                   or die "Could not connect to $RunParameters->{'AnnotateHost'}";
    $RunParameters->{'markersDb'} = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} );
    $RunParameters->{'filesCollection'} = $RunParameters->{'markersDb'}->get_collection('files');
    $RunParameters->{'uidCollection'}   = $RunParameters->{'markersDb'}->get_collection( $RunParameters->{'uidName'} );
    $RunParameters->{'locksCollection'} = $RunParameters->{'markersDb'}->get_collection('locks');
    $RunParameters->{'locksCollection'}->ensure_index( { 'bufferIndex' => 1 } );
    $RunParameters->{'filesCollection'}->ensure_index( { 'filepath'    => 1 } );
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        $RunParameters->{$collectionDb}= $RunParameters->{'markersDb'}->get_collection($collectionDb);
    }
    $RunParameters->{'annotateByCoordCollection'} = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->get_collection('coord');
    $RunParameters->{'annotateByRulesCollection'} = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->get_collection('rules');
    $RunParameters->{'annotateByGeneCollection'}  = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->get_collection('gene');
    $RunParameters->{'annotateByChrPosCollection'}  = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->get_collection('chrPos');
    $RunParameters->{'annotateByGeneCodonCollection'}  = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->get_collection('geneCodon');

    if ( exists( $RunParameters->{'BufferHostList'} ) ) {
        undef $RunParameters->{'bufferHost'};
        delete $RunParameters->{'bufferHost'};
        for ( my $index = 0 ; $index <= $#{ $RunParameters->{'BufferHostList'} } ; ++$index ) {
            for ( my $i = 1 ; $i <= $RunParameters->{'BufferHostList'}->[$index]->{'connections'} ; ++$i ) {
                push( @{ $RunParameters->{'bufferHost'} }, $RunParameters->{'BufferHostList'}->[$index]->{'host'} );
            }
        }
    }
    if ( exists( $RunParameters->{'bufferHost'} ) ) {
        for ( my $index = 0 ; $index <= $#{ $RunParameters->{'bufferHost'} } ; ++$index ) {
            if ( $RunParameters->{'locksCollection'}->count( { 'bufferIndex' => $index } ) == 0 ) {
                $RunParameters->{'locksCollection'}->insert( { 'bufferHost' => $RunParameters->{'bufferHost'}->[$index],
                                                               'status' => int(1),
                                                               'bufferIndex' => $index,
                                                               'hostname' => 'none',
                                                               'pid' => '0'
                                                             },
                                                             { safe => 1 }
                                                           );
            }
        }
    }
    else {
        die "!!!!!Died: No bufferHost given in configuration file!!!!!!\n";
    }
    return $RunParameters;
}

sub loadConf {
    my $function      = shift;
    my $confFile      = shift;
    my $RunParameters = {};
    open( CONF, $confFile ) or die "Can't find $confFile\n";
    while (<CONF>) {
        if (/^#/) { next }
        chomp;
        if (/(.*)=(.*)/) { $RunParameters->{$1} = $2 }
        print "Conf: $1=$2\n";
    }
    close(CONF);
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    return $RunParameters;
}

sub removeCollection {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    $RunParameters->{'collectionConn'}->remove();
    $RunParameters->{'filesCollection'}->remove();
    return 1;

}

sub dropDb {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    print "\n--------------Dropping Markers-------------- \n";
    $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->drop();
    Maintain->closeAllConnections($RunParameters);
    Maintain->mongoConnect($RunParameters);
    print "\t\t-Dropped Markers and Reset Db - re-establishing connections\n";
    return 1;
    die;
}

sub loadPad {
    my $function = shift;
    my $mul      = shift;
    ++$mul;
    my $pad = "";
    for ( my $i = 0 ; $i <= $mul * 1000 ; ++$i ) {
        $pad = $pad . "t";
    }
    return $pad;
}

sub buildConnection {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    my $in      = shift;
    my $timeOut = 1000000;
  LOOP: for ( my $i = 0 ; $i <= $#{ $RunParameters->{'bufferHost'} } ; ++$i ) {
        my $bufferName = "";
        my $bufferDb   = "";
        if ( $RunParameters->{'bufferHost'}->[$i] =~ /(.+?)\/(.*)/ ) {
            $bufferName = $1;
            $bufferDb   = $2;
        }
        unless ( length($bufferName) > 0 && length($bufferDb) > 0 ) { next LOOP }
        my $cursor = $RunParameters->{'markersDb'}->run_command(
            {
                findAndModify => 'locks',
                query         => { 'status' => int(1), "bufferIndex" => int($i) },
                update        => { '$set' => { 'status' => 0, 'in' => $in, "hostname" => hostname, "pid" => "$$" } }
            }
        );
        if ( $cursor->{'ok'} eq "1" && exists( $cursor->{'value'}->{'bufferIndex'} ) ) {
            my $doc = $cursor->{'value'};
            $RunParameters->{'bufferIndex'}->{$in} = $doc->{'bufferIndex'};
            $RunParameters->{'bufferConn'}->{$in}  = MongoDB::MongoClient->new( host => "mongodb://$bufferName", query_timeout => $timeOut, timeout => $timeOut, wtimeout => $timeOut );
            $RunParameters->{'bufferDb'}->{$in}    = $RunParameters->{'bufferConn'}->{$in}->get_database($bufferDb);
            my $DatabaseName = $in . $doc->{'bufferIndex'} . hostname . $$;
            $DatabaseName =~ s/-//g;
            $RunParameters->{'vcfCollection'}->{$in} = $RunParameters->{'bufferConn'}->{$in}->get_database($bufferDb)->get_collection($DatabaseName);
            return 1;
        }
    }
    die "Sleeping: Everything locked!!\n";
    sleep(60);
    goto LOOP;
}

sub closeConnection {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    my $in = shift;
    $RunParameters->{'vcfCollection'}->{$in}->drop();
    $RunParameters->{'locksCollection'}->update( { 'bufferIndex' => int( $RunParameters->{'bufferIndex'}->{$in} ) }, { '$set' => { 'status' => int(1), 'in' => "None", "hostname" => "none", "pid" => "None" } }, { safe => 1 } );
    undef $RunParameters->{'vcfCollection'}->{$in};
    delete $RunParameters->{'vcfCollection'}->{$in};
    undef $RunParameters->{'bufferConn'}->{$in};
    delete $RunParameters->{'bufferConn'}->{$in};
    undef $RunParameters->{'bufferDb'}->{$in};
    delete $RunParameters->{'bufferDb'}->{$in};
    undef $RunParameters->{'bufferIndex'}->{$in};
    delete $RunParameters->{'bufferIndex'}->{$in};
}

sub closeAllConnections {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    my $timeOut = 1000000;
    print "\t+Brutally Closing All Connections & Deleting Buffers\n";
    undef $RunParameters->{'vcfCollection'};
    delete $RunParameters->{'vcfCollection'};
    $RunParameters->{'locksCollection'}->remove( {}, { safe => 1 } );
    foreach my $bufferHost ( @{ $RunParameters->{'bufferHost'} } ) {
        my $bufferName = "";
        my $bufferDb   = "";
        if ( $bufferHost =~ /(.+?)\/(.*)/ ) {
            $bufferName = $1;
            $bufferDb   = $2;
        }
        my $conn = MongoDB::MongoClient->new( host => "mongodb://$bufferName", query_timeout => $timeOut, timeout => $timeOut, wtimeout => $timeOut );
        $conn->get_database($bufferDb)->drop();
    }
    undef $RunParameters->{'vcfCollection'};
    delete $RunParameters->{'vcfCollection'};
    undef $RunParameters->{'bufferConn'};
    delete $RunParameters->{'bufferConn'};
    undef $RunParameters->{'bufferDb'};
    delete $RunParameters->{'bufferDb'};
    undef $RunParameters->{'bufferIndex'};
    delete $RunParameters->{'bufferIndex'};
    Maintain->mongoConnect($RunParameters);
    print "\t\t-Closed all connections and reconnected\n";
}

sub buildBufferConnections {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        Maintain->buildConnection( $RunParameters, 'nakedInsert_' . $collectionDb );
        Maintain->buildConnection( $RunParameters, 'geneInsert_' . $collectionDb );
        Maintain->buildConnection( $RunParameters, 'staged_' . $collectionDb );
    }
    return 1;
}

sub closeBufferConnections {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    print "\t+Gently Closing Buffers\n";
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        Maintain->closeConnection( $RunParameters, 'nakedInsert_' . $collectionDb);
        Maintain->closeConnection( $RunParameters, 'geneInsert_' . $collectionDb);
        Maintain->closeConnection( $RunParameters, 'staged_' . $collectionDb);
    }
    return 1;
}

sub checkConditions {
    my $function           = shift;
    my $variantRef        = shift;
    my $filter            = shift;
    my $ref               = shift;
    my $success           = shift;
    my $fail              = shift;
    if ($ref->{'oper'} eq "EQUALS"      ) {
        if ( $variantRef->{$filter} eq $ref->{'val1'}) {return $success;}
    } elsif ($ref->{'oper'} eq "NOT"    ) {
        if ( $variantRef->{$filter} ne $ref->{'val1'}) {return $success;}
    } elsif ($ref->{'oper'} eq "CONDITIONAL") {
         if (exists( $variantRef->{$ref->{'fieldIfTrue'}} ) && exists( $variantRef->{$ref->{'fieldIfFalse'}} )) {
             my $condRef={};
             $condRef->{'field'}=$ref->{'field'};
             $condRef->{'oper'}=$ref->{'conditOperator'};
             $condRef->{'val1'}=$ref->{'conditVal1'};
             if ($condRef->{'oper'} eq "BETWEEN") {$condRef->{'val2'}=$condRef->{'conditVal2'}};
             if (checkConditions( $variantRef, $condRef->{'field'}, $condRef,  'pass', 'fail' ) eq 'pass') {
                my $trueRef={};
                $trueRef->{'field'}=$ref->{'fieldIfTrue'};
                $trueRef->{'oper'}=$ref->{'operIfTrue'};
                $trueRef->{'val1'}=$ref->{'val1IfTrue'};
                if ($trueRef->{'oper'} eq "BETWEEN") {$trueRef->{'val2'}=$ref->{'val2IfTrue'}};
                if (checkConditions( $variantRef, $trueRef->{'field'}, $trueRef, 'pass', 'fail') eq 'pass') {
                    return $success;
                } else {return $fail}
              } else {
                my $falseRef={};
                $falseRef->{'field'}=$ref->{'fieldIfFalse'};
                $falseRef->{'oper'}=$ref->{'operIfFalse'};
                $falseRef->{'val1'}=$ref->{'val1IfFalse'};
                if ($falseRef->{'oper'} eq "BETWEEN") {$falseRef->{'val2'}=$ref->{'val2IfFalse'}};
                if (Maintain->checkConditions( $variantRef, $falseRef->{'field'}, $falseRef, 'pass', 'fail' ) eq 'pass') {
                    return $success;
                }    else {return $fail}
             }
         }
    } elsif ($ref->{'oper'} eq "AND") {
         if (exists( $variantRef->{$ref->{'field2'}})) {
              my $andRef={};
              $andRef->{'field'}=$ref->{'field'};
              $andRef->{'oper'}=$ref->{'operForField'};
              $andRef->{'val1'}=$ref->{'val1ForField'};
              if ($andRef->{'oper'} eq "BETWEEN") {$andRef->{'val2'}=$ref->{'val2ForField'}};
              if (Maintain->checkConditions( $variantRef, $andRef->{'field'}, $andRef, 'pass', 'fail' ) eq 'pass') {
                  my $andRef2={};
                  $andRef2->{'field'}=$ref->{'field2'};
                  $andRef2->{'oper'}=$ref->{'operForField2'};
                  $andRef2->{'val1'}=$ref->{'val1ForField2'};
                  if ($andRef2->{'oper'} eq "BETWEEN") {$andRef2->{'val2'}=$ref->{'val2ForField2'}};
                  if (Maintain->checkConditions( $variantRef, $andRef->{'field'}, $andRef, 'pass', 'fail') eq 'pass') {
                      return $success;
                  }    else {return $fail}
              }    else {return $fail}
         }
    } elsif ($ref->{'oper'} eq "REGEX") {
        if (exists($ref->{'val1'})) {
           if ( $variantRef->{$filter} =~/$ref->{'val1'}/) {return $success;}
        }
    } elsif (looks_like_number($variantRef->{$filter} )) {
        if ($ref->{'oper'} eq "<"      ) {
           if ( $variantRef->{$filter} < $ref->{'val1'}) {return $success;}
        } elsif ($ref->{'oper'} eq "<="     ) {
            if ( $variantRef->{$filter} <= $ref->{'val1'}) {return $success;}
        } elsif ($ref->{'oper'} eq ">"      ) {
            if ( $variantRef->{$filter} > $ref->{'val1'}) {return $success;}
        } elsif ($ref->{'oper'} eq ">="     ) {
            if ( $variantRef->{$filter} >= $ref->{'val1'}) {return $success;}
        } elsif ($ref->{'oper'} eq "=="     ) {
            if ( $variantRef->{$filter} == $ref->{'val1'}) {return $success;}
        } elsif ($ref->{'oper'} eq "BETWEEN") {
            if (exists($ref->{'val2'})) {
                if (
                     ($variantRef->{$filter} > $ref->{'val1'} && $variantRef->{$filter} < $ref->{'val2'}) ||
                     ($variantRef->{$filter} > $ref->{'val2'} && $variantRef->{$filter} < $ref->{'val1'})
                   ) {return $success;}
            }
        }
   } elsif ($ref->{'oper'} eq "DEF") {
           if ( exists($variantRef->{$filter})) {return $success;}
   }
   return $fail;
}

sub printSystemDefaults {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    print "\n________________________ System Defaults_______________________________________ \n";
    $Data::Dumper::Indent = 0;
    $Data::Dumper::Pair   = " : ";
    my $str = Dumper $RunParameters;
    $str =~ s/\$VAR/\$RunParameters/g;
    $str =~ s/\s//g;
    $str =~ s/\,/, /g;
    $str =~ s/\'//g;

    for ( my $index = 0 ; $index <= $#{ $RunParameters->{'BufferHostList'} } ; ++$index ) {
        my $val = $RunParameters->{'BufferHostList'}->[$index]->{'host'};
        $str =~ s/, $val//g;
    }
    print "sysRef: $str\n_________________________________________________________________________________________ \n\n";
    return $RunParameters;
}

sub reloadStartupRunParameters {
    my $function      = shift;
    my $RunParameters = shift;
    foreach my $runParameter ( keys %{ $RunParameters->{'StartupRunParameters'} } ) {
        $RunParameters->{$runParameter} = $RunParameters->{'StartupRunParameters'}->{$runParameter};
    }    
    if ( exists( $RunParameters->{'CreateOnPath'} ) ) {
        unless ( $RunParameters->{'CreateOnPath'} =~ /\/$/ ) {
            $RunParameters->{'CreateOnPath'} = $RunParameters->{'CreateOnPath'} . "/";
        }
    }
    if ( exists( $RunParameters->{'vcfPath'} ) ) {
        unless ( $RunParameters->{'vcfPath'} =~ /\/$/ ) {
            $RunParameters->{'vcfPath'} = $RunParameters->{'vcfPath'} . "/";
        }
    }
    return $RunParameters;
}

sub loadLocations {
    my $function      = shift;
    my $RunParameters = shift;
    open( DB, "$RunParameters->{'gemGeneLocPath'}" ) or warn "!Can't find $RunParameters->{'gemGeneLocPath'}\n";
    print "\t\t+Loading Locations\n";
    my $isGene = retrieve("$RunParameters->{'gemGeneLocPath'}.isGene.HASH");
    my $isExon = retrieve("$RunParameters->{'gemGeneLocPath'}.isExon.HASH");
    if (exists($RunParameters->{'customAnnotation'})) {
        foreach my $entry (keys %{$RunParameters->{'customAnnotation'}}) {
            my $file=$RunParameters->{'customAnnotation'}->{$entry};
            open (FILE, "$file") or die "\t!!!!Can't open file: $file!!!\n";
            LOOP: while (<FILE>) { 
               my $line=$_;
               chomp($line);
               if ($line=~/^#/) {next LOOP};
               my @fields=split("\t",$line);
               if ($#fields<7) { die "\t!! $RunParameters->{'customAnnotation'}->{$entry} has the wrong number of fields, should be 7 (see $line) becomes @fields\n";}
               my $gene=$entry;
               my $chr=$fields[0];
               $chr=~s/chr//g;
               if ($line=~/[;|\t]GENE=(.*?)\;/) {
                   $gene=$1;
               } elsif ($line=~/[;|\t]gene=(.*?)\;/) {
                   $gene=$1;
               }
               if ($line=~/[;|\t]END=(.*?)\;/) {
                  my $enSt=$1;
                  my $st=int($fields[1]/100)-1;
                  my $en=int($enSt/100)+1;
                  for ( my $i = $st ; $i <= $en ; ++$i ) {
                      my $pos = $i * 100;
                      if ( exists( $isGene->{"chr$chr:$pos"} ) ) {
                          unless ( $isGene->{"chr$chr:$pos"} =~ /;$gene;/ ) {
                              $isGene->{"chr$chr:$pos"} = $isGene->{"chr$chr:$pos"} . "$gene;";
                          } else {
                              $isGene->{"chr$chr:$pos"} = ";$gene;";                         
                          }
                      }
                      else {
                         $isGene->{"chr$chr:$pos"} = ";$gene;";                                        
                      }

                   }
                   $st=int($fields[1]/10)-1;
                   $en=int($enSt/10)+1;
                   for ( my $i = $st ; $i <= $en ; ++$i ) {
                       my $pos = $i * 10;
                       $isExon->{"chr$chr:$pos"} = int(1);
                   }
               } else  {
                   my $st=int($fields[1]/100)-1;
                   my $en=int($fields[1]/100)+1;
                   for ( my $i = $st ; $i <= $en ; ++$i ) {
                       my $pos = $i * 100;
                       if ( exists( $isGene->{"chr$chr:$pos"} ) ) {
                           unless ( $isGene->{"chr$chr:$pos"} =~ /;$gene;/ ) {
                               $isGene->{"chr$chr:$pos"} = $isGene->{"chr$chr:$pos"} . "$gene;";
                           }else {
                              $isGene->{"chr$chr:$pos"} = ";$gene;";                         
                          }
                       }
                       else {
                           $isGene->{"chr$chr:$pos"} = ";$gene;";
                       }
                   }
                   $st=int($fields[1]/10)-1;
                   $en=int($fields[1]/10)+1;
                   for ( my $i = $st ; $i <= $en ; ++$i ) {
                       my $pos = $i * 10;
                       $isExon->{"chr$chr:$pos"} = int(1);
                   }
               }
            }
        }
    }
    return ( $isExon, $isGene );
    #my $isExon = {};
    #my $isGene = {};
    while (<DB>) {
        if (/^#/) { next }
        chomp;
        my @fields = split( /\t/, $_ );
        if ( $fields[1] ne "protein_coding" ) { next }
        my $chr = $fields[0];

        if ( $fields[2] eq "gene" ) {
            my ( $st, $en ) = 0;
            my $geneName = 1;
            if ( $fields[4] > $fields[3] ) {
                $st = int( $fields[3] / 100 ) - 1;
                $en = int( $fields[4] / 100 ) + 1;
            }
            else {
                $st = int( $fields[4] / 100 ) - 1;
                $en = int( $fields[3] / 100 ) + 1;
            }
            if ( $fields[8] =~ /gene_name \"(.*?)\"\; gene_source/ ) {
                $geneName = $1;
            }
            for ( my $i = $st ; $i <= $en ; ++$i ) {
                my $pos = $i * 100;
                if ( exists( $isGene->{"chr$chr:$pos"} ) ) {
                    unless ( $isGene->{"chr$chr:$pos"} =~ /;$geneName;/ ) {
                        $isGene->{"chr$chr:$pos"} = $isGene->{"chr$chr:$pos"} . "$geneName;";
                    }
                }
                else {
                    $isGene->{"chr$chr:$pos"} = ";$geneName;";
                }
            }
        }
        elsif ( $fields[2] eq "exon" || $fields[2] eq "UTR") {
            my ( $st, $en ) = 0;
            if ( $fields[4] > $fields[3] ) {
                $st = int( $fields[3] / 10 ) - 2;
                $en = int( $fields[4] / 10 ) + 2;
            }
            else {
                $st = int( $fields[4] / 10 ) - 2;
                $en = int( $fields[3] / 10 ) + 2;
            }
            for ( my $i = $st ; $i <= $en ; ++$i ) {
                my $pos = $i * 10;
                $isExon->{"chr$chr:$pos"} = int(1);
            }
        }
    }

#    store $isExon, "$RunParameters->{'gemGeneLocPath'}.isExon.HASH";
#    store $isGene, "$RunParameters->{'gemGeneLocPath'}.isGene.HASH";
#die;
    return ( $isExon, $isGene );
}
sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}
1;
