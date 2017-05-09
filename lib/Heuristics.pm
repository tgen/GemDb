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
# * Author (s):
#    David Craig
#    Release 9/1/15
#
#  'status' 1: valid and completed, 2: set to scrub 3:  4: Moved through InsertVar
#/
######################################################################################
package Heuristics;
use strict;
use Time::localtime;
use File::stat;
use Scalar::Util qw(looks_like_number);
use MongoDB;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);


$| = 1;
sub readPed {
   # my $function      = shift;
    my $RunParameters = shift;
    my $ped  = shift;
    my $fields={'ped'=>$ped,'validPed'=>'0'};
    if ($ped=~/(^\w*?)_(\w*?)_(\w*?)_(\w*?)_(\w*?)_F(....)G(.)S(.)(M.*?)(P.*?)_(\w*)/) {
         $fields->{'validPed'}=1;
         print "\t\t+Valid Ped Format: $ped\n";
         if (length($1)>0) {
            $fields->{'studyFam'}=$1;
         } else {die "\tDie Unknown $1\n"}
         if (length($2)>0) {
            $fields->{'patientId'}=$2;
         }else {die "\tDie Unknown $2\n"}
         if (length($3)>0) {
            $fields->{'visitSampleNumber'}=$3;
         }  else {die "\tDie Unknown $3\n"}    
         if (length($4)>0) {
            $fields->{'tissue'}=$4;
         } else {die "\tDie Unknown $4\n"}  
         if (length($6)>0) {
            $fields->{'family'}=$6;
         }else {die "\tDie Unknown $6\n"}           
         if (length($7)>0) {
            $fields->{'gender'}=$7;
            if ($fields->{'gender'}eq "1") {$fields->{'gender'}="Male"}else{$fields->{'gender'}="Female"}
         }else {die "\tDie Unknown $7\n"}
         if (length($8)>0) {
            $fields->{'affectedStatus'}=$8;
            if ($fields->{'affectedStatus'}eq "1") {
                $fields->{'affectedStatus'}="Unaffected";
             }elsif($fields->{'affectedStatus'}eq "2") {
                 $fields->{'affectedStatus'}="Affected";
             }               
         } else {die "\tDie Unknown $1\n"}
         if (length($9)>0) {
            $fields->{'momId'}=$9;
         } else {die "\tDie Unknown $9\n"}
         if (length($10)>0) {
            $fields->{'dadId'}=$10;
            $fields->{'dadId'}=~s/P/D/;
         }  else {die "\tDie Unknown $10\n"}
         if (length($11)>0) {
            $fields->{'assay'}=$11;
         }    else {die "\tDie Unknown $11\n"}    
         if ( $fields->{'affectedStatus'} eq "Affected") {
            $fields->{'proband'}='S'.$fields->{'studyFam'};                                                                                                     
         } elsif ($fields->{'affectedStatus'} eq "Unaffected") {
         }
     }
     return $fields;
}

sub joinVCF {
    my $function      = shift;
    my $RunParameters = shift;
    print "--------------------- joining VCFs $function----\n";
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();    
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    print "\t\t+Verifying indexed studyPatientTissue/gene\n";      
    $RunParameters->{'markersDb'}->get_collection('tumor')->ensure_index( Tie::IxHash->new('variants.filepath' => 1, 'gene' => 1) );        
    my $cursor = $RunParameters->{'filesCollection'}->find({'joinGeneInMarkers'=>1});
    $cursor->immortal(1);
    
    LOOP: while (my $doc=$cursor->next) {
        print "\t+Joining file: $doc->{'filepath'}\n";
        my $filepath=$doc->{'filepath'};
        my $studyPatientTissue=$doc->{'studyPatientTissue'};
        my $cursor2 = $RunParameters->{'filesCollection'}->find({'studyPatientTissue'=>$studyPatientTissue});
        my @filepaths=();
        my $collectionDb=$doc->{'collectionDb'};         
        while (my $doc2=$cursor2->next) {
           if ($doc2->{'filepath'} ne $doc->{'filepath'}) {
                  push (@filepaths,$doc2->{'filepath'});   
            }                                  
        }
        my $joinInName="";
        if (exists($doc->{'joinInName'})) {
           $joinInName=$doc->{'joinInName'}; 
        }             
        my $markerConn1 = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
        my $cursor3=$markerConn1->find({'variants.studyPatientTissue'=>$studyPatientTissue});
        my %genesLookup=();
        while (my $doc3=$cursor3->next) {
           $genesLookup{$doc3->{'gene'}}=1;
        }
        
        my $infoFields={};
        foreach my $field (keys %{$doc}) {
            if ($field=~/joinInfo(.*)/) {               
               if ($doc->{$field}==1 || $doc->{$field}==2 || $doc->{$field}==3) {
                  $infoFields->{$1}=$doc->{$field};
               }
            }
        }
        open (FILE,$filepath) or die "\t\t!Can't open $filepath\n";
        LOOP2: while (<FILE>) {
            my $line=$_;
            if ($line=~/^#/) {next LOOP2}
            chomp($line);
            my @fields=split(/\t/,$line);
            my $info="";
            if (defined($fields[7])) {
              $info=$fields[7];
            }
            my $alt=$fields[4];
            my $qual=$fields[5];
            my $filter=$fields[6];
            my @infoFields=split(/\;/,$info);
            my @genes=();
            if ($info=~/GENE=(.*?)\;/) {
               my $infoGene=$1;
               foreach my $gene (split(/\,/,$infoGene)) {
                  if (length($gene)>0) { 
                      push(@genes,$gene)
                  }
               }
            } elsif ($info=~/GENE=(.*?)/) {
               my $infoGene=$1;
               foreach my $gene (split(/\,/,$infoGene)) {
                  if (length($gene)>0) {
                      push(@genes,$gene)
                  }
               }
            }
            foreach my $infoField (@infoFields) {
                if ($infoField=~/(.*?)=(.*)/) {
                    my $key=$1;
                    my $val=$2;
                    LOOP3: foreach my $gene (@genes) {
                        unless (exists($genesLookup{$gene})) {next LOOP3};
                         my $insert=0;#change
                         if (exists($infoFields->{$key})) {
                            if ($infoFields->{$key}==3) {
                               if ( looks_like_number($val)) {
                                  $val=int($val);
                                  $insert=1;
                               }
                            } elsif ($infoFields->{$key}==2) {
                               if ( looks_like_number($val)) {
                                  MongoDB::force_double($val);
                                  $insert=1;
                               }
                     
                            } elsif ($infoFields->{$key}==1) {
                                  $insert=1;                          
                            }
                         }
                         if ($insert==1) {                             
                             foreach my $filepathToAdd (@filepaths) {
                                     $markerConn1->update( Tie::IxHash->new('variants.filepath' => $filepathToAdd , 'gene' => $gene),
                                                           { '$set' => { 'variants.$.' . $joinInName . $key => $val } }, 
                                                           { safe => 1, multiple=>1} );                       
                              }
                         }
                    }
                }
            }
        }
        close (FILE);
        $RunParameters->{'filesCollection'}->update({'filepath'=>$filepath},{'$set'=>{'joinGeneInMarkers'=>0}},{'safe'=>1,'multiple'=>1});
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t+Joined in $diff Seconds\n";    
}

sub indexCollections {
    my $function      = shift;
    my $RunParameters = shift;
    $RunParameters = Maintain->reloadStartupRunParameters($RunParameters);
    print "\t+Building Indexes...";
    my @vals = (
        'biomarker', 'gene', 'variants.filepath','impact','report',
        'info.checkCount', "variants.$RunParameters->{uidName}",
        'dbFreqCat',  'info.group','variants.studyPatientTissue','studyPatient','variants.studyPatient',
        'info.varFilled',  'info.toDelete','info.md5sum','geneInfo.chr','geneInfo.startPos','geneInfo.endPos',
        'aberration.aberration_name','aberration.aberration_filter','variants.filter'
    );
    $RunParameters->{'markersDb'}->get_collection('tumor')->ensure_index( Tie::IxHash->new('info.group' => 1, 'info.varFilled' => 1) );        
    $RunParameters->{'markersDb'}->get_collection('germline')->ensure_index( Tie::IxHash->new('info.group' => 1, 'info.varFilled' => 1) );        
    
    $RunParameters->{'markersDb'}->get_collection('germline')->ensure_index( Tie::IxHash->new('variants.projectRun' => 1, 'impact' => 1, 'dbFreqCat'=>1) );        
    #$RunParameters->{'markersDb'}->get_collection('tumor')->ensure_index( Tie::IxHash->new('variants.projectRun' => 1, 'impact' => 1, 'dbFreqCat'=>1) );        
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    foreach my $val (@vals) {
        for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
            my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};    
            $RunParameters->{'markersDb'}->get_collection($collectionDb)->ensure_index( { $val => 1 } );
        }
    }
#    my @valsGerm = ( 'annotate.maxPopFreq' );
#    foreach my $val (@valsGerm) {
#        $RunParameters->{'markersDb'}->get_collection('germline')->ensure_index( { $val => 1 } );
#    }
    my @valsTum = ('rules_checked_array','rules_checked');
    foreach my $val (@valsTum) {
        $RunParameters->{'markersDb'}->get_collection('tumor')->ensure_index( { $val => 1 } );
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "Indexed in $diff Seconds\n";
    return 1;
}

sub discoverFileFields {
##################################
    # discoverFileFields
##################################
    my ( $function, $RunParameters, $filepath, $filedate, $curVcfPath, $fileInfoRef ) = @_;

    $fileInfoRef->{'style'}                        = 'unknown';
    if ( $fileInfoRef->{'filepath'} =~ /reports\/(.*?)\/(.*?)\/(.*?)\/(.*?)\/(.*?)\/analysisResults\/(.*?\.vcf)$/ ) {
        $fileInfoRef->{'study'}                           = $1;
        $fileInfoRef->{'patient'}                         = $2;
        $fileInfoRef->{'collectionDate'}                  = $3;
        $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} = $4;
        $fileInfoRef->{'assay'}                           = $5;
        $fileInfoRef->{'filename'}                        = $6;
        $fileInfoRef->{'style'}                        = 'reports';
        $fileInfoRef->{ $RunParameters->{'uidName'} } =
          join( "_", $fileInfoRef->{'study'}, $fileInfoRef->{'patient'}, $fileInfoRef->{'collectionDate'}, $fileInfoRef->{'sampleSubgroupVisitSampleNumber'}, $fileInfoRef->{'assay'} );
        if ( $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} =~ /(.)(.)/ ) {
            $fileInfoRef->{'sampleSubgroup'}    = $1;
            $fileInfoRef->{'VisitSampleNumber'} = $2;
        }
        if ( $fileInfoRef->{'patient'} =~ /(.*?)s(.*?)p(.*)/ ) {
            $fileInfoRef->{'referenceSample'} = $1;
            $fileInfoRef->{'library'}         = $2;
            $fileInfoRef->{'pool'}            = $3;
        }
    }
    elsif ( $fileInfoRef->{'filepath'} =~ /reports\/(.*?)\/(.*?)\/(.*?)\/(.*?)\/(.*?)\/(.*?\.vcf)$/ ) {
        $fileInfoRef->{'study'}                           = $1;
        $fileInfoRef->{'patient'}                         = $2;
        $fileInfoRef->{'collectionDate'}                  = $3;
        $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} = $4;
        $fileInfoRef->{'assay'}                           = $5;
        $fileInfoRef->{'filename'}                        = $6;
        $fileInfoRef->{'style'}                        = 'reports';
        $fileInfoRef->{ $RunParameters->{'uidName'} } = join( "_", $fileInfoRef->{'study'}, $fileInfoRef->{'patient'}, $fileInfoRef->{'collectionDate'}, $fileInfoRef->{'sampleSubgroupVisitSampleNumber'}, $fileInfoRef->{'assay'} );
        if ( $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} =~ /(.)(.)/ ) {
            $fileInfoRef->{'sampleSubgroup'}    = $1;
            $fileInfoRef->{'VisitSampleNumber'} = $2;
        }
        if ( $fileInfoRef->{'patient'} =~ /(.*?)s(.*?)p(.*)/ ) {
            $fileInfoRef->{'referenceSample'} = $1;
            $fileInfoRef->{'library'}         = $2;
            $fileInfoRef->{'pool'}            = $3;
        }
    }
    elsif ( $fileInfoRef->{'filepath'} =~ /reports\/(.*?)\/(.*?)\/(.*?)\/(.*?)\/(.*?\.vcf)$/ ) {
        $fileInfoRef->{'study'}                           = $1;
        $fileInfoRef->{'patient'}                         = $2;
        $fileInfoRef->{'collectionDate'}                  = "none";
        $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} = $3;
        $fileInfoRef->{'assay'}                           = $4;
        $fileInfoRef->{'filename'}                        = $5;
        $fileInfoRef->{'style'}                        = 'reports';
        $fileInfoRef->{ $RunParameters->{'uidName'} } = join( "_", $fileInfoRef->{'study'}, $fileInfoRef->{'patient'}, $fileInfoRef->{'collectionDate'}, $fileInfoRef->{'sampleSubgroupVisitSampleNumber'}, $fileInfoRef->{'assay'} );
        if ( $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} =~ /(.)(.)/ ) {
            $fileInfoRef->{'sampleSubgroup'}    = $1;
            $fileInfoRef->{'VisitSampleNumber'} = $2;
        }
        if ( $fileInfoRef->{'patient'} =~ /(.*?)s(.*?)p(.*)/ ) {
            $fileInfoRef->{'referenceSample'} = $1;
            $fileInfoRef->{'library'}         = $2;
            $fileInfoRef->{'pool'}            = $3;
        }
    }
    elsif ( $fileInfoRef->{'filepath'} =~ /reports\/(.*?)\/(.*?)\/(.*?)\/(.*?\.vcf)$/ ) {
        $fileInfoRef->{'study'}          = $1;
        $fileInfoRef->{'style'}                        = 'reports';
        $fileInfoRef->{'patient'}        = $2;
        $fileInfoRef->{'collectionDate'} = "none";
        $fileInfoRef->{'tissue'}         = $3;
        $fileInfoRef->{'filename'}       = $4;
        if ( $fileInfoRef->{'tissue'} eq "germline" ) {
            $fileInfoRef->{'collectionDb'} = "germline";
            $fileInfoRef->{'assay'}        = "germline";
        }
        else {
            $fileInfoRef->{'collectionDb'} = "tumor";
            $fileInfoRef->{'assay'}        = "tumor";
        }
        $fileInfoRef->{ $RunParameters->{'uidName'} } = join( "_", $fileInfoRef->{'study'}, $fileInfoRef->{'patient'} );
    }
    elsif ( $fileInfoRef->{'filepath'} =~ /reports\/(.*)\/(.*?\.vcf)$/ ) {
        $fileInfoRef->{ $RunParameters->{'uidName'} } = $1;
        $fileInfoRef->{'style'}                        = 'reports';
        $fileInfoRef->{ $RunParameters->{'uidName'} } =~ s/.*\///;
        $fileInfoRef->{'filename'} = $2;
        if ( $fileInfoRef->{ $RunParameters->{'uidName'} } =~ /^(.+?)_(.+?)_(.+?)_(.+?)_(.*)/ ) {
            $fileInfoRef->{'study'}                           = $1;
            $fileInfoRef->{'patient'}                         = $2;
            $fileInfoRef->{'recordId'}                        = $3;
            $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} = $4;
            if ( exists( $RunParameters->{'reportableGenes'}->{$5} ) ) {
                $fileInfoRef->{'assay'} = $5;
            }
            else {
                if (exists($RunParameters->{'defaultAssay'})) {
                    $fileInfoRef->{'assay'} =  $RunParameters->{'defaultAssay'};
                } else {
                    $fileInfoRef->{'assay'}='tumor';
                }
            }
        }
    } else {
        $fileInfoRef->{'study'}   = "NoStudy";
        $fileInfoRef->{'patient'} = "NoPatient";
        $fileInfoRef->{'style'}                        = 'VCF';
        if ( $fileInfoRef->{'filepath'} =~ /VCF\/(.*?)\/(.*?)\/(.*?)\/(.*.vcf)$/ ) {
            $fileInfoRef->{'study'}    = $1;
            $fileInfoRef->{'patient'}  = $2;
            $fileInfoRef->{'sample'}   = $3;
            $fileInfoRef->{'tissue'}   = $3;
            $fileInfoRef->{'filename'} = $4;
            if ( $fileInfoRef->{'sample'} eq "germline" ) { $fileInfoRef->{'collectionDb'} = 'germline'; $fileInfoRef->{'assay'} = 'germline';} else {$fileInfoRef->{'assay'} = 'tumor';}
            $fileInfoRef->{'valid'}        = int(1);                        
        }
        elsif ( $fileInfoRef->{'filepath'} =~ /VCF\/(.*?)\/(.*?)\/(.*.vcf)$/ ) {
            $fileInfoRef->{'study'}    = $1;
            $fileInfoRef->{'patient'}  = "$1.$2";
            $fileInfoRef->{'filename'} = $3;
            $fileInfoRef->{'valid'}        = int(1);            
        }
        elsif ( $fileInfoRef->{'filepath'} =~ /VCF\/(.*?)\/(.*.vcf)$/ ) {
            $fileInfoRef->{'study'}    = $1;
            $fileInfoRef->{'patient'}  = "$1.$2";
            $fileInfoRef->{'filename'} = $2;
            $fileInfoRef->{'valid'}        = int(1);                        
        }
        elsif ( $fileInfoRef->{'filepath'} =~ /VCF\/(.*.vcf)$/ ) {
            $fileInfoRef->{'study'}    = "no_study";
            $fileInfoRef->{'patient'}  = $1;
            $fileInfoRef->{'filename'} = $1;            
        }
        elsif ( $fileInfoRef->{'filepath'} =~ /VCF\/.*\/(.*?.vcf)$/ ) {
            $fileInfoRef->{'filename'} = $1;
            if ( $fileInfoRef->{'filename'} =~ /.*\/(.*?)\.vcf$/ ) {
                $fileInfoRef->{ $RunParameters->{'uidName'} } = $1;
                if ( $fileInfoRef->{ $RunParameters->{'uidName'} } =~ /^(.*?)_(.*?)_(.*?)_(.*?)_(.*?)/ ) {
                    $fileInfoRef->{'study'}   = $1;
                    $fileInfoRef->{'patient'} = "$1.$2";
                }
                else {
                    $fileInfoRef->{'study'}   = "no_study";
                    $fileInfoRef->{'patient'} = $fileInfoRef->{ $RunParameters->{'uidName'} };
                }
            }
            else {
                $fileInfoRef->{'study'}   = "no_study";
                $fileInfoRef->{'patient'} = $fileInfoRef->{'filename'};
            }
            $fileInfoRef->{'sample'} = "tumor";
            if ( !exists( $fileInfoRef->{ $RunParameters->{'uidName'} } ) ) {
                die "\t!Please provide a $RunParameters->{'uidName'}, I can't figure it out.\n";
            }
        }
        unless ( exists( $fileInfoRef->{'sample'} ) ) {
            $fileInfoRef->{'sample'} = 'tumor';
            $fileInfoRef->{'assay'}       = "tumor";
        }
        if ( $fileInfoRef->{'sample'} =~ /germ/i && !( exists( $fileInfoRef->{'collectionDb'} ) ) ) {
            $fileInfoRef->{'sample'}       = "germline";
            $fileInfoRef->{'assay'}       = "germline";
            $fileInfoRef->{'collectionDb'} = "germline";
        }
        $fileInfoRef->{'collectionDate'}                  = '01012000';
        $fileInfoRef->{'sampleSubgroupVisitSampleNumber'} = $fileInfoRef->{'sample'};
        $fileInfoRef->{'tissue'}                          = $fileInfoRef->{'sample'};
    }
    $fileInfoRef->{'studyPatient'} = $fileInfoRef->{'study'} . "." . $fileInfoRef->{'patient'};
    unless ( exists( $fileInfoRef->{ $RunParameters->{'uidName'} } ) ) {
        $fileInfoRef->{ $RunParameters->{'uidName'} } = $fileInfoRef->{'studyPatient'};
    }
    return $fileInfoRef;

}    #    validateFile(sysRef,filepath)=fileInfoRef

sub runInheritance {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb  = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my %projectRunPs   = ();
    my %toSkip   = ();    
    my $projectRunPNum = 0;
    my $fileConn   = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection('files');
#    my $cursor     = $fileConn->find( { 'filepath' => '/IVY/VCF/crdc/crdc0193/germline/C4RCD0193_C1_2_1438343-C4RCD0193_M1_2_1438344-C4RCD0193_D1_1_1438345.HC_All.snpEff.vcf', 'inheritanceCheck' => { '$ne' => 'yes' } } );    
    my $cursor     = $fileConn->find( { 'study' => 'crdc' , 'inheritanceCheck' => { '$ne' => 'yes' }} );    
    print "\t+Run Inheritance $function\n";    
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
  PATIENT: while ( my $doc = $cursor->next ) { 
        if ( exists( $doc->{'sampleNames'} ) ) {
            my $val = scalar( @{ $doc->{'sampleNames'} } );
            $projectRunPNum = $projectRunPNum + $val;
        }
        $projectRunPs{$doc->{'projectRun'}}{inheritance}=1;
        if ( exists( $doc->{'sampleNames'}->[2] ) )  {
            my $mom   = 0;
            my $momId = "";            
            my $dad   = 0;
            my $dadId = "";
            my $child = 0;
            my $validPedCount=0;
            my $ln= scalar(@{$doc->{'sampleNames'}});
            my $ped=readPed($RunParameters,$doc->{'sampleNames'}->[0]); 
            if ($ped->{'validPed'} eq "1") {   
                my $match=-1;        
                LOOP:for ( my $i = 0 ; $i <$ln; ++$i ) {
                    my $pos= $i+1;  
                    my $mem=readPed($RunParameters,$doc->{'sampleNames'}->[$i]);
                    print "patientId: $mem->{'patientId'} $mem->{'affectedStatus'} mom: $mem->{'momId'}  dad: $mem->{'dadId'}\n";
                    if (!exists(   $ped->{'patientId'}) || !exists($mem->{'patientId'}) ||!exists($mem->{'affectedStatus'}) || $mem->{'validPed'} ne "1" ) {
                        print "\t\tWarning, missing $i $doc->{'sampleNames'}->[$i] with $ped->{'ped'} $doc->{'projectRun'}\n";
                        $match=-1;
                        last LOOP;
                    } elsif ( $mem->{'affectedStatus'}eq "Affected") {
                        if (length($mem->{'momId'})>0  && length($mem->{'dadId'})>0){
                            $match=$i;
                            $projectRunPs{ $doc->{'projectRun'} }{'child'} = 's' . $pos; 
                            ++$child; 
                            $momId=     $mem->{'momId'};
                            $dadId=     $mem->{'dadId'};                        
                        }
                    }
                }
                if ($match>-1) {
                LOOP1:for ( my $i = 0 ; $i <$ln; ++$i ) {    
                        my $pos= $i+1;
                        my $mem=readPed($RunParameters,$doc->{'sampleNames'}->[$i]);                 
                        unless ($i==$match) {            
                            if ("M" .$mem->{'patientId'} eq  $momId) {
                                $projectRunPs{ $doc->{'projectRun'} }{'mom'} = 's' . $pos; 
                                  ++$mom;   
                            } elsif ("D" .$mem->{'patientId'} eq $dadId) {
                                $projectRunPs{ $doc->{'projectRun'} }{'dad'} = 's' . $pos; 
                                  ++$dad;       
                            }            
                        }
                    }
                }
            } else {
                for ( my $i = 0 ; $i <= 2 ; ++$i ) {
                    my $pos= $i+1;
                    if ( $doc->{'sampleNames'}->[$i] =~ /C1/ ) { 
                        $projectRunPs{ $doc->{'projectRun'} }{'child'} = 's' . $pos; 
                        ++$child; 
                    }
                    if ( $doc->{'sampleNames'}->[$i] =~ /D1/ || $doc->{'sampleNames'}->[$i] =~ /F1/ ) { 
                        $projectRunPs{ $doc->{'projectRun'} }{'dad'} = 's' . $pos; 
                        ++$dad;  
                    }
                    if ( $doc->{'sampleNames'}->[$i] =~ /M1/ ) { 
                        $projectRunPs{ $doc->{'projectRun'} }{'mom'} = 's' . $pos; 
                        ++$mom; 
                    }
                }
            }
            unless ($mom>=1 && $dad>=1 && $child>=1) {
                print "\t\tMissing some of trio : mom=$mom, dad=$dad,child=$child\n";                                                                              
                $toSkip{$doc->{'projectRun'}}=1;
            }
        }
    }
    print "\t\t-Loaded Patients: $projectRunPNum\n";
    my $markerConn = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
    $markerConn->ensure_index( Tie::IxHash->new('variants.projectRun' => 1, 'impact' => 1, 'dbFreqCat'=>1) );        
   # $markerConn->ensure_index( Tie::IxHash->new('gene'=>1,'dbFreqCat' => 1, 'impact' => 1,'variants.genotype.s1.gt'=>1) );          
    my $markerConn1 = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
    my $markerConn2 = $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection($collectionDb);
    foreach my $projectRunP ( reverse sort keys %projectRunPs ) {
        print "\t\t-Patient: $projectRunP\n";
        my $cursor = $markerConn1->find(  Tie::IxHash->new(
                                           'variants' => { '$elemMatch' => { 'projectRun'=> $projectRunP } } ,
                                           'impact'    => 'significant',
                                           'dbFreqCat' => 'rare')
                                        );
        $cursor->fields( { biomarker => 1, 'variants.$' => 1, snpeff => 1, dbFreq => 1,impact=>1,dbFreqCat=>1 } );
        my $fp = "";
        $cursor->immortal(1);
        VARIANT: while ( my $doc = $cursor->next ) {
            my $inheritance = "NA";
            my $id          = $doc->{'_id'};
            if ( my $gene = $doc->{snpeff}->{gene} ) {
                my $biomarker = $doc->{biomarker};
                my $dbFreq  = $doc->{dbFreq};
                #my $geneCount = $markerConn2->count(Tie::IxHash->new(
                #        'gene'                    => $gene,
                #        'dbFreqCat'               => 'rare',
                #        'impact'           => 'significant',
                #        'variants'                  => { '$elemMatch' => { projectRun => $projectRunP, "variants.projectRun.genotype.s1.gt" => { '$ne' => "0/0" } } } )
                #);                      
                #$markerConn1->update( { '_id' => $id, "variants.projectRun" => $projectRunP }, 
                #                      { '$set' => { 'variants.$.rareGeneCount' => int($geneCount) } }, 
                #                      { safe => 1} );
                if (exists($toSkip{$projectRunP})) {next VARIANT}                
                unless ( exists( $projectRunPs{$projectRunP}{'child'} ) ) {
                    next VARIANT;
                }
                unless ( exists( $doc->{variants}->[0]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} ) ) {
                    next VARIANT;
                }
                my $geneCount =     int(0);     
                for (my $m=0; $m<scalar(@{$doc->{variants}});++$m) {
                    unless ( exists( $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} ) 
                             && exists( $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} ) 
                           ) { next VARIANT }
                    if ( $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} eq "0/0" ) { 
                         $inheritance = "Reference";                 
                    } elsif ( $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} eq "0/0" && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} eq "0/0" ) {
                        $geneCount = $markerConn2->count( {  gene => $gene, 
                                                  dbFreqCat => 'rare', 
                                               #   impact   => 'significant',
                                                  variants => { '$elemMatch' => { projectRun => $projectRunP, "genotype.$projectRunPs{$projectRunP}{child}.gt" => { '$ne' => "0/0" } } } } );
                        if ( $doc->{variants}->[$m]->{filter} eq "PASS" ) {
                            $inheritance = "De Novo";
                        }    else {
                            $inheritance = "LQ De Novo";
                        }
                    }    elsif ($doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} eq "1/1"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} eq "0/0"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} eq "0/1" )
                    {
                        $inheritance = "X-Linked";
                    }    elsif ($doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} eq "1/1"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} eq "0/1"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} eq "0/1" )
                    {
                        $inheritance = "Homozygous";
                    }   elsif ($doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} eq "0/1"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} eq "0/1"
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} eq "0/1" ) {
                        $inheritance = "CommonInFamily";                    
                    } elsif ($doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'child'} }->{gt} eq $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt}
                        && $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'dad'} }->{gt} eq $doc->{variants}->[$m]->{genotype}->{ $projectRunPs{$projectRunP}{'mom'} }->{gt} ) {
                        $inheritance = "CommonInFamily";             
                    } else {
                        my $momPassCount = $markerConn2->count(
                            {
                                'gene'                    => $gene,
                                dbFreqCat               => 'rare',
                                'impact'           => 'significant',
                                variants => { '$elemMatch' => { projectRun => $projectRunP, "genotype.$projectRunPs{$projectRunP}{mom}.gt" => "0/1", "genotype.$projectRunPs{$projectRunP}{dad}.gt" => "0/0", "genotype.$projectRunPs{$projectRunP}{child}.gt" => "0/1", filter => "PASS" } }
                            }
                        );
                        my $dadPassCount = $markerConn2->count(
                            {
                                'gene'                    => $gene,
                                dbFreqCat               => 'rare',
                                'impact'           => 'significant',
                                variants => { '$elemMatch' => { projectRun => $projectRunP, "genotype.$projectRunPs{$projectRunP}{dad}.gt" => "0/1", "genotype.$projectRunPs{$projectRunP}{mom}.gt" => "0/0", "genotype.$projectRunPs{$projectRunP}{child}.gt" => "0/1", filter => "PASS" } }
                            }
                        );
                        my $momCount = $markerConn2->count(
                            {
                                'gene'                    => $gene,
                                dbFreqCat               => 'rare',
                            'impact'           => 'significant',
                                variants                  => { '$elemMatch' => { projectRun => $projectRunP, "genotype.$projectRunPs{$projectRunP}{mom}.gt" => "0/1", "genotype.$projectRunPs{$projectRunP}{dad}.gt" => "0/0", "genotype.$projectRunPs{$projectRunP}{child}.gt" => "0/1" } }
                            }
                        );
                        my $dadCount = $markerConn2->count(
                            {
                                'gene'                    => $gene,
                                dbFreqCat               => 'rare',
                                'impact'           => 'significant',
                                variants                  => { '$elemMatch' => { projectRun => $projectRunP, "genotype.$projectRunPs{$projectRunP}{dad}.gt" => "0/1", "genotype.$projectRunPs{$projectRunP}{mom}.gt" => "0/0", "genotype.$projectRunPs{$projectRunP}{child}.gt" => "0/1" } }
                            }
                        );
                        if ( $momCount > 0     && $dadCount > 0 )     { $inheritance = "LQ Phased Compount Het"; }
                        if ( $momPassCount > 0 && $dadPassCount > 0 ) { $inheritance = "Phased Compount Het"; }
                        $fp = $doc->{variants}->[$m]->{'filepath'};
                    }
                    $markerConn1->update( { '_id' => $id, "variants.projectRun" => $projectRunP }, 
                                          { '$set' => { 'variants.$.inheritance' => $inheritance, 'variants.$.rareGeneCount' => int($geneCount) } }, 
                                          { safe => 1 } 
                                        );
                }
            }
        }
        $fileConn->update( {  'projectRun'=> $projectRunP }, { '$set' => { 'inheritanceCheck' => "yes" } }, { safe => 1 ,multiple=>1} );
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t-runInheritance: $diff seconds\n";
    1;
}

sub addRules {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    unless ( exists( $RunParameters->{'RulesFiles'} ) ) { die "\t! I can't find a RulesFiles Parameter\n" }
    my %rulesSaw = ();
    print "--------------------- Adding Rules ----\n";
    for ( my $i = 0 ; $i <= $#{ $RunParameters->{'RulesFiles'} } ; ++$i ) {
        my $rulesVersion = $RunParameters->{'RulesFiles'}->[$i];
        if ( $rulesVersion eq "*" ) {
            my $holder = $RunParameters->{'annotateByConn'}->get_database( $RunParameters->{'AnnotateDbName'} )->run_command( [ 'distinct' => 'rules', 'key' => 'RulesVersion', 'query' => {} ] );
            $RunParameters->{'dynamicRuleUpdating'}=1;
            foreach my $rule ( @{ $holder->{'values'} } ) {
                $rulesSaw{$rule} = 1;
                print "\t\t+Will run rules: $rule\n";
            }
        }
        else {
            $rulesSaw{$rulesVersion} = 1;
            print "\t\t+Will run rules: $rulesVersion\n";
        }
        print "\t+Loading Rules ($rulesVersion) to new files\n";
        my %lastFilepath = ();
        my %lastUid      = ();
        $RunParameters->{$collectionDb}->ensure_index( { 'rules_checked_array' => 1 } );
        $RunParameters->{$collectionDb}->ensure_index( { 'rules_checked_array' => 1, 'CreateOnPath' => 1 } );
        $RunParameters->{$collectionDb}->ensure_index( { 'rules_checked'       => 1 } );
        my $rulesConn = $RunParameters->{'annotateByRulesCollection'};
        foreach my $newRules ( keys %rulesSaw ) {
            my $biomarkerCursor;
            ## Determine how to look up rules
            if ( exists( $RunParameters->{'CreateOnPath'} ) ) {
                $biomarkerCursor = $RunParameters->{$collectionDb}->find( { 'rules_checked_array' => { '$nin' => [$newRules] }, 'CreateOnPath' => $RunParameters->{'CreateOnPath'} } );
            }
            elsif ( exists( $RunParameters->{'dynamicRuleUpdating'} ) ) {
                if ( $RunParameters->{'dynamicRuleUpdating'} == 1 ) {
                    $biomarkerCursor = $RunParameters->{$collectionDb}->find( { 'rules_checked_array' => { '$nin' => [$newRules] }});
                }
                else {
                    $biomarkerCursor = $RunParameters->{$collectionDb}->find( { 'rules_checked' => { '$exists' => 0 } } );
                }
            }
            else {
                $biomarkerCursor = $RunParameters->{$collectionDb}->find( { 'rules_checked' => { '$exists' => 0 } } );
            }
            my $c = 0;
            print "\t\t+Adding $newRules rules to biomarkers in $collectionDb\n";
          MARKER: while ( my $biomarker = $biomarkerCursor->next ) {
                ++$c;
                if ( $c % 10000 == 0 ) { print "\t\t- Count at $c\n"; }
                my $id           = $biomarker->{'_id'};
                my @matchingRule = ();
                if ( exists( $biomarker->{'variants'}->['0']->{'filepath'} ) ) {
                    $lastFilepath{ $biomarker->{'variants'}->['0']->{'filepath'} } = 1;
                }
                else {
                    print "\t\t\tSlight warning: $biomarker->{'biomarker'} seems to be missing expected filepath\n";
                    next MARKER;
                }
                $lastUid{ $biomarker->{ $RunParameters->{'uidName'} } } = 1;
                unless ( exists( $biomarker->{'aberration'}->{'lookup'} ) ) { next MARKER }
                my $lookup = uc( $biomarker->{'aberration'}->{'lookup'} );
                unless ( exists( $biomarker->{'aberration'}->{'gene'} ) && exists( $biomarker->{'aberration'}->{'aberration_type'} ) ) { next MARKER }
                my $changeLookup = uc( $biomarker->{'aberration'}->{'gene'} . "_" . $biomarker->{'aberration'}->{'aberration_type'} . "_CHANGE" );
                my $wildcard     = 0;
                my $noWildcard   = 0;
                my $snv          = 0;
                my $fusion       = 0;
                my $cursor1      = $rulesConn->find( { 'aberration_lookup' => { '$in' => [ $lookup, $changeLookup ] }, 'drug_rules_version'=>$newRules } );
                while ( my $rule = $cursor1->next ) {
                    my $finalRule = Annotate->joinOn( '_id', $rule->{'_id'}, $RunParameters->{'annotateByRulesCollection'} );
                    if ( exists( $finalRule->{'drug_rule_snv_wild_card_match'} ) || exists( $finalRule->{'drug_rule_fusion_wild_card_match'} ) ) {
                        $wildcard = 1;
                        if ( exists( $finalRule->{'drug_rule_snv_wild_card_match'} ) )    { $snv    = 1 }
                        if ( exists( $finalRule->{'drug_rule_fusion_wild_card_match'} ) ) { $fusion = 1 }
                    }
                    push( @matchingRule, $finalRule );
                }
                my @finalMatchingRule = ();
                foreach my $rule (@matchingRule) {
                    if ( $wildcard == 1 ) {
                        if ( exists( $rule->{'drug_rule_snv_wild_card_match'} ) || exists( $rule->{'drug_rule_fusion_wild_card_match'} ) ) {
                            push( @finalMatchingRule, $rule );
                        }
                    }
                    else {
                        push( @finalMatchingRule, $rule );
                    }
                }
                $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set' => { 'drug_rule_matched_flag' => int(0) } }, { 'safe' => 1 } );
                my $count = scalar(@finalMatchingRule);    #print "count:$count @finalMatchingRule\n";
                if ( scalar(@finalMatchingRule) > 0 ) {
                    if ( exists( $finalMatchingRule[0]->{'aberration_lookup'} ) ) {
                        $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set' => { 'matching_rule' => \@finalMatchingRule } }, { 'safe' => 1 } );
                        if ( $wildcard == 1 && $snv == 1 )    { 
                             $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set' => { 'drug_rule_snv_wild_card_match'    => int(1) } }, { 'safe' => 1 } ) 
                        }
                        if ( $wildcard == 1 && $fusion == 1 ) { 
                              $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set' => { 'drug_rule_fusion_wild_card_match' => int(1) } }, { 'safe' => 1 } ) 
                        }
                        $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set' => { 'drug_rule_matched_flag' => int(1) } }, { 'safe' => 1 } );
                    }
                }
                $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$set'  => { 'rules_checked'       => int(1) } },    { 'safe' => 1 } );
                $RunParameters->{$collectionDb}->update( { '_id' => $id }, { '$push' => { 'rules_checked_array' => $newRules } }, { 'safe' => 1 } );       
            }     
        }
        foreach my $newRules ( keys %rulesSaw ) {
            foreach my $filepath ( keys %lastFilepath ) {
                $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection('files') ->update( 
                                                              { 'filepath' => $filepath, 'rules_checked_array' => { '$nin' => [$newRules] } }, 
                                                              { '$addToSet' => { 'rules_checked_array' => $newRules } }, { 'safe' => 1, 'multiple' => 1 } 
                                                            );
            }
            foreach my $uid ( keys %lastUid ) {
                $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $uid, 'rules_checked_array' => { '$nin' => [$newRules] } }, 
                                                           { '$addToSet' => { 'rules_checked_array' => $newRules } }, 
                                                           { 'safe' => 1, 'multiple' => 1 } 
                                                         );
            }
        }
    }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\n------  Add Rules Completed for ($diff seconds)----  \n";
    return 1;
}


1;
