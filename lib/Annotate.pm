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
# *
# * Author(s):
#    David Craig
#    Release 9/1/15
#/
######################################################################################

package Annotate;
use strict;
use File::stat;
use MongoDB;
use MongoDB::OID;
use Time::localtime;
use Digest::MD5 qw(md5 md5_hex md5_base64);

$| = 1;

#### MAIN #################
sub annotate {
    my $function      = shift;
    my $RunParameters = shift;
    Maintain->reloadStartupRunParameters($RunParameters);
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    for (my $i=0;$i<=$RunParameters->{'CollectionCount'};++$i) {
        my $collectionDb=$RunParameters->{'Collections'}->[$i]->{'CollectionName'};
        print "\n\n------Annotation for $collectionDb ------ \n";
        Annotate->annotateCustom( $RunParameters, $collectionDb );
        Annotate->runSnpeff( $RunParameters, $collectionDb );
        Annotate->annotateSnpeff( $RunParameters, $collectionDb );
        Annotate->annotateGene( $RunParameters, $collectionDb );
        if ($collectionDb eq 'tumor') {
            Annotate->addRulesByFilename( $RunParameters, $collectionDb );
        }
        Annotate->finalize( $RunParameters, $collectionDb );
        my $stat->{'elapsedInsertTime'} = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $RunParameters->{'GlobalStart'};
        if ( exists( $RunParameters->{'addedFilepaths'} ) ) {
            print "\t+Completed files:\n";
            for ( my $i = 0 ; $i <= $#{ $RunParameters->{'addedFilepaths'} } ; ++$i ) {
                my $uid     = $RunParameters->{'addedUids'}->[$i];
                my $fp      = $RunParameters->{'addedFilepaths'}->[$i];
                #$RunParameters->{'filesCollection'}->update( { 'filepath' => $fp },{ '$set' => { 'status' => int(1)}}, { safe => 1, 'multiple' => 1 } );
                my $cursor2 = $RunParameters->{'filesCollection'}->find( { 'filepath' => $fp, 'collectionDb' => $collectionDb } );
                if ( my $doc2 = $cursor2->next ) {
                    if ( $doc2->{'valid'} eq "1" ) {
                        my $stat->{'End_ps'} = '';#`ps -o cmd,%cpu,%mem`;
                        $stat->{'elapsedInsertTime'} = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $RunParameters->{'GlobalStart'};
                        $RunParameters->{'uidCollection'}->update(
                            { $RunParameters->{'uidName'} => $uid },
                            {'$set'=>{
                                $RunParameters->{'uidName'} => $uid,
                                'inserted'                  => int(-1),
                                'countSNV'                  => int(-1),
                                'countSV'                   => int(-1),
                                'countCNA'                  => int(-1),
                                'countPointMutations'       => int(-1),
                                'countDifferentialExpressedGenes'                  => int(-1),
                                'countFusions'                  => int(-1),
                                'version'                   => $RunParameters->{'version'},
                                'study'                     => $doc2->{'study'},
                                'patient'                   => $doc2->{'patient'},
                                'studyPatient'              => $doc2->{'studyPatient'},
                                'studyPatientTissue'                     => $doc2->{'studyPatientTissue'},
                                'vcfPath'                     => $doc2->{'vcfPath'},
                                'assay'                     => $doc2->{'assay'},
                                'style'                     => $doc2->{'style'},
                                'status'                     => int(0),
                                'reportable'                     => int(0),
                                'compute' => { 'end' => $stat }
                                }
                            },
                            { safe => 1, upsert => 1 }
                        );
                        my $fileInsertCount  = int(-1);
                        my $myUidInsertCount = int(-1);
                        $RunParameters->{'filesCollection'}->update( { 'filepath' => $doc2->{'filepath'} },
                                                                     { '$set' => { 'status' => int(1), 'compute.end' => $stat, 'inserted' => $fileInsertCount } },
                                                                     { safe => 1, 'multiple' => 1 }
                                                                    );
                        if ( exists( $RunParameters->{'writeVcfTo'} ) ) {
                            Maintain->writeVCF( $RunParameters, $fp, $RunParameters->{'writeVcfTo'} );
                        }
                        $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $uid },
                                                                   { '$addToSet' => { 'filepaths' => $fp } },
                                                                   { safe => 1 }
                                                                );
                        print "\t+Successful Added file: $fp\n";
                    } else {
                        $RunParameters->{'filesCollection'}->update( { 'filepath' => $doc2->{'filepath'} },
                                                                     { '$set' => { 'status' => int(1), 'compute.end' => $stat, 'inserted' => int(0) } },
                                                                     { safe => 1, 'multiple' => 1 }
                                                                    );
                    }
                }
            }
        }
        else {
            print "\t+No Added files\n";
        }
        my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
        print "\n------  Annotation Completed for $collectionDb ($diff seconds)----  \n";
    }
}

sub addRulesByFilename {
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
            foreach my $rule ( @{ $holder->{'values'} } ) {
                $rulesSaw{$rule} = 1;
                print "\t+Will run rules: $rule\n";
            }
        }
        else {
            $rulesSaw{$rulesVersion} = 1;
            print "\t+Will run rules: $rulesVersion\n";
        }
        my %lastFilepath = ();
        for ( my $i = 0 ; $i <= $#{ $RunParameters->{'addedFilepaths'} } ; ++$i ) {
            $lastFilepath{$RunParameters->{'addedFilepaths'}->[$i]}=1;
        }
        my %lastUid      = ();
        my $conn  = $RunParameters->{'vcfCollection'}->{ 'staged_' . $collectionDb };
        my $rulesConn = $RunParameters->{'annotateByRulesCollection'};
        foreach my $newRules ( keys %rulesSaw ) {
            my $biomarkerCursor = $conn->find();
            my $c = 0;
            MARKER: while ( my $biomarker = $biomarkerCursor->next ) {
                ++$c;
                my $id           = $biomarker->{'_id'};
                my @matchingRule = ();
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
                        unless ( exists( $rule->{'drug_rule_snv_wild_card_match'} ) && exists( $rule->{'drug_rule_fusion_wild_card_match'} ) ) {
                            push( @finalMatchingRule, $rule );
                        }
                    }
                    else {
                        push( @finalMatchingRule, $rule );
                    }
                }
                unless (exists ($biomarker->{ 'drug_rule_matched_flag'})) {
                    $conn->update( { '_id' => $id }, { '$set' => { 'drug_rule_matched_flag' => int(0) } }, { 'safe' => 1 } );
                }
                my $count = scalar(@finalMatchingRule);    #print "count:$count @finalMatchingRule\n";
                if ( scalar(@finalMatchingRule) > 0 ) {
                    for (my $i=0;$i< $count;++$i) {
                        if ( exists( $finalMatchingRule[$i]->{'aberration_lookup'} ) ) {
                            $conn->update( { '_id' => $id }, { '$push' => { 'matching_rule' => $finalMatchingRule[$i] } }, { 'safe' => 1 } );
                        }
                    }
                    if ( $wildcard == 1 && $snv == 1 )    {
                         $conn->update( { '_id' => $id }, { '$set' => { 'drug_rule_snv_wild_card_match'    => int(1) } }, { 'safe' => 1 } )
                    }
                    if ( $wildcard == 1 && $fusion == 1 ) {
                          $conn->update( { '_id' => $id }, { '$set' => { 'drug_rule_fusion_wild_card_match' => int(1) } }, { 'safe' => 1 } )
                    }
                    $conn->update( { '_id' => $id }, { '$set' => { 'drug_rule_matched_flag' => int(1) } }, { 'safe' => 1 } );
                }
                $conn->update( { '_id' => $id }, { '$set'  => { 'rules_checked'       => int(1) } },    { 'safe' => 1 } );
                $conn->update( { '_id' => $id }, { '$push' => { 'rules_checked_array' => $newRules } }, { 'safe' => 1 } );
            }
        }
        foreach my $newRules ( keys %rulesSaw ) {
            foreach my $filepath ( keys %lastFilepath ) {
                $RunParameters->{'mongoConnection'}->get_database( $RunParameters->{'DatabaseName'} )->get_collection('files') ->update(
                                                              { 'filepath' => $filepath },
                                                              { '$addToSet' => { 'rules_checked_array' => $newRules } }, { 'safe' => 1, 'multiple' => 1 }
                                                            );
            }
            foreach my $uid ( keys %lastUid ) {
                $RunParameters->{'uidCollection'}->update( { $RunParameters->{'uidName'} => $uid },
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
sub annotateCustom {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    if (exists($RunParameters->{'customAnnotation'})) {
        my $start        = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
        print "\t+Annotating/Assigning Custom $collectionDb\n";
        my $conn      = $RunParameters->{'vcfCollection'}->{ 'nakedInsert_' . $collectionDb };
        my $connStage = $RunParameters->{'vcfCollection'}->{ 'staged_' . $collectionDb };
        $conn->ensure_index( { 'biomarker' => 1 } );
        my $custAnnRef          = {};
        my $counted=0;
        foreach my $entry (keys %{$RunParameters->{'customAnnotation'}}) {
             my $file=$RunParameters->{'customAnnotation'}->{$entry};
             open (FILE, "$file") or die "\t!!!!Can't open file: $file!!!\n";
             LOOP: while (<FILE>) {
                my $line=$_;
                chomp($line);
                if ($line=~/^#/) {next LOOP};
                my @fields=split("\t",$line);
                my $chr=$fields[0];
                $chr=~s/chr//g;
                my $pos=$fields[1];
                my $id=$fields[2];
                my $ref=$fields[3];
                my $alt=$fields[4];
                my $qual=$fields[4];
                my $filter=$fields[6];
                my $info=$fields[7];
                my $gene=$entry;
                my $end=$pos;
                if ($line=~/[;|\t]END=(.*?)\;/) {
                    $end=$1;
                }
                for (my $i=$pos;$i<=$end;++$i) {
                    if ($ref eq ".") {
                        $ref="A,T,C,G";
                    }
                    if ($alt eq ".") {
                        $alt="A,T,C,G";
                    }
                    my @refs=split(",",$ref);
                    my @alts=split(",",$alt);
                    foreach my $refPos (@refs) {
                       foreach my $altPos (@alts) {
                           my $coord="chr$chr:$i:$refPos:$altPos";
                           my @infoFields=split(";",$info);
                           foreach my $infoField (@infoFields) {
                               if ($infoField=~/(.*?)=(.*)/) {
                                   my $key=$1; 
                                   my $val=$2;
                                   if (length($key)>0 && length($val)>0) {
                                       if ($key eq "gene" || $key eq "GENE") {
                                          $custAnnRef->{$coord}->{'gene'}=$val;
                                       }elsif($key eq "effect"){
                                          $custAnnRef->{$coord}->{'effect'}=$val;
                                       }elsif($key eq "report"){
                                          $custAnnRef->{$coord}->{'report'}=$val;
                                       }elsif($key eq "impact"){
                                          $custAnnRef->{$coord}->{'impact'}=$val;
                                       }elsif($key eq "type") {
                                          $custAnnRef->{$coord}->{'type'}=$val;
                                       }
                                       if ($key=~/aberration_(.*)/) {
                                          my $aberr=$1;
                                          if (length($aberr)>0) {
                                             $custAnnRef->{$coord}->{'aberration'}->{$aberr}=$val;
                                          }
                                       }
                                   }
                               }
                           }
                           unless (exists($custAnnRef->{$coord}->{'gene'})) {
                               $custAnnRef->{$coord}->{'gene'}="none";
                           };
                           unless (exists($custAnnRef->{$coord}->{'aberration'}->{'custom'})) {$custAnnRef->{$coord}->{'aberration'}->{'aberration_custom'}=$entry};
                           unless (exists($custAnnRef->{$coord}->{'aberration'}->{'name'})) {$custAnnRef->{$coord}->{'aberration'}->{'aberration_name'}='Custom Annotation'};
                           unless (exists($custAnnRef->{$coord}->{'aberration'}->{'gene'})) {$custAnnRef->{$coord}->{'gene'}=$entry};
                           unless (exists($custAnnRef->{$coord}->{'effect'})) {$custAnnRef->{$coord}->{'effect'}="none"};
                           unless (exists($custAnnRef->{$coord}->{'impact'})) {$custAnnRef->{$coord}->{'impact'}="low"};
                           unless (exists($custAnnRef->{$coord}->{'report'})) {$custAnnRef->{$coord}->{'report'}="comprehensive"};
                           unless (exists($custAnnRef->{$coord}->{'type'})) {$custAnnRef->{$coord}->{'type'}="SNV"};
                           my $cursor = $conn->find( { 'biomarker' => $coord } );
                           $cursor->immortal(1);
                          LOOP: while ( my $biomarkerRef = $cursor->next ) {
                                 my $gene       = $custAnnRef->{$coord}->{'gene'};
                                 my @ProjectRun = ();
                                 my @variants   = ();
                                 my $id=$biomarkerRef->{'_id'};
                                # my $biomarkerRef = dclone($doc);
#                                my $biomarkerRef=$doc;
                                 foreach my $variantRef ( @{ $biomarkerRef->{'variants'} } ) {
                                     $variantRef->{'gene'}=$gene;
                                     push( @variants, $variantRef );
                                     push (@ProjectRun,  $variantRef->{ $RunParameters->{'uidName'} });
                                 }
                                $biomarkerRef->{'aberration'}       = $custAnnRef->{$coord}->{'aberration'};
                                unless (exists($biomarkerRef->{'aberration'}->{'aberration_filter'})) {
                                    $biomarkerRef->{'aberration'}->{'aberration_filter'}="PASS";
                                }
                                $biomarkerRef->{'gene'}         = $gene;
                                $biomarkerRef->{'report'}         = $custAnnRef->{$coord}->{'report'};
                                $biomarkerRef->{'impact'}         = $custAnnRef->{$coord}->{'impact'};
                                $biomarkerRef->{'significance'} = "None";
                                $biomarkerRef->{'type'} = "SNV";
                                if (exists($RunParameters->{'AnnotateByCoordFields'})) {
                                     if ( $RunParameters->{'AnnotateByCoordFields'}->{'ALL'}==0) {
                                         $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                                  $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                                  $RunParameters->{'annotateByCoordCollection'},
                                                                                                  $RunParameters->{'AnnotateByCoordFields'}
                                                                                                 );
                                     } else {
                                         $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                                  $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                                  $RunParameters->{'annotateByCoordCollection'}
                                                                                                 );
                                     }
                                 } else {
                                         $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                                  $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                                  $RunParameters->{'annotateByCoordCollection'}
                                                                                                 );
                                 }
                                $biomarkerRef->{'toAnnotateByCoord'} = int(0);
                                if (exists($RunParameters->{'AnnotateByGeneFields'})) {
                                    if ( $RunParameters->{'AnnotateByGeneFields'}->{'ALL'}==0) {
                                         $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'},$RunParameters->{'AnnotateByGeneFields'}  );
                                     } else {
                                         $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'}  );
                                     }
                                } else {
                                         $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'} );
                                }
                               $biomarkerRef->{'toAnnotateByGene'}            = int(0);
                               $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRun;
                               $biomarkerRef->{'info'}->{'varMax'}                      = int(1);
                               if ( ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) > $biomarkerRef->{'info'}->{'varMax'} && ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) < $RunParameters->{'varGlobalMax'} ) {
                                   $biomarkerRef->{'info'}->{'varMax'} = int( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 4 );
                               }
                               elsif ( ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) >= $RunParameters->{'varGlobalMax'} ) {
                                   $biomarkerRef->{'info'}->{'varMax'} = int( $RunParameters->{'varGlobalMax'} );
                               }
                               if ( $biomarkerRef->{'info'}->{'varArray'} >= $biomarkerRef->{'info'}->{'varMax'} ) {
                                   $biomarkerRef->{'info'}->{'varFilled'} = int(1);
                               } else {
                                   $biomarkerRef->{'info'}->{'varFilled'} = int(0);
                               }
                               $biomarkerRef->{'variants'}                    = \@variants;
                               $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRun;
                               $biomarkerRef                                  = Annotate->isCommon($biomarkerRef);
                               delete( $biomarkerRef->{'_id'} );
                               $connStage->insert( $biomarkerRef, { safe => 1 } );
                               ++$counted;print "$coord $gene\n";
                               $conn->remove({ '_id' => $id} );
                          }
                        }
                    }
                }
             }
         }
         my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;    
         print "\t\t-Done with Custom Annotation for $diff seconds, adding $counted\n";
    }
    
}
sub runSnpeff {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    print "\t+runSnpeff for $collectionDb \n";
    my ( $chr, $pos, $ref, $alt ) = "";
    my %done  = ();
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $conn  = $RunParameters->{'vcfCollection'}->{ 'nakedInsert_' . $collectionDb };
    $conn->ensure_index( { 'biomarker' => 1 } );
    open( SNPEFF, ">$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf" ) or die "Can't open $RunParameters->{'runDir'}/snpeff.$collectionDb.tmp.vcf";
    print SNPEFF "##fileformat=VCFv4.1\n";
    print SNPEFF "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO\n";
    my $cursor = $conn->find();
    $cursor->immortal(1);
    my $count = 0;
  LOOP: while ( my $doc = $cursor->next ) {
        ++$count;
        unless ( $chr = $doc->{'variants'}->[0]->{'chr'} )     { next LOOP }
        unless ( $pos = $doc->{'variants'}->[0]->{'hg19pos'} ) { next LOOP }
        unless ( $ref = $doc->{'variants'}->[0]->{'ref'} )     { next LOOP }
        unless ( $alt = $doc->{'variants'}->[0]->{'alt'} )     { next LOOP }
        my $coord = "$chr:$pos:$ref:$alt";
        if ( $alt =~ /\</ ) { next LOOP }
        if ( exists( $done{ $doc->{'variants'}->[0]->{'coord'} } ) ) { next LOOP }
        $done{ $doc->{'variants'}->[0]->{'coord'} } = 1;
        my $line = "";

        if ( exists( $doc->{'collectionDb'} ) ) {
            $line = "$chr\t$pos\t.\t$ref\t$alt\t40\tPASS\tversion=original;variantCollection=$doc->{'collectionDb'};final=1\n";
            print SNPEFF $line;
        }
    }
    close SNPEFF;
    if ( $count > 0 ) {
        SLEEP: my $numRunning=`ps -ef | grep "$RunParameters->{'Snpeff4Path'}" | grep snpEff\.jar | wc -l `;
        if ($numRunning >1) {
                sleep(60);
                print "\t\t+Sleeping because these are running: \n";
                goto SLEEP;
        }
        my $err = `java -Xmx32G -jar $RunParameters->{'Snpeff4Path'}/snpEff.jar eff -c $RunParameters->{'Snpeff4Path'}/snpEff.config $RunParameters->{'SnpeffParameter'} -onlyProtein -quiet -noStats -noLog -o vcf -t $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf > $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf 2> $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf.err`;
        if (stat("$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf.err")->size > 5) {
           print "\t\t++Snpeff with $err $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf.err\n";
           `rm -f $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf.err`;
        }  else {
           `rm -f $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf.err`;
        }
    }
    `rm -f $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf`;
    print "\t\t+rm -f $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf\n";
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t-Running SnpEff Completed for $collectionDb: $diff seconds\n";
}

sub annotateSnpeff {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    my %ids          = ();
    my $start        = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    print "\t+Annotating/Assigning Snpeff $collectionDb\n";
    if ( -e "$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf" ) {
        open( SNPEFF, "$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf" ) or die "Can't find $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf\n";
    }
    else {
        print "\t\t-Nothing run for annotateSnpeff ($RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf)\n";
        return 0;
    }
    my $conn      = $RunParameters->{'vcfCollection'}->{ 'nakedInsert_' . $collectionDb };
    my $connStage = $RunParameters->{'vcfCollection'}->{ 'staged_' . $collectionDb };
    $conn->ensure_index( { 'biomarker' => 1 } );
    my $geneListRef = {};
    foreach my $libraryToPrune ( keys %{ $RunParameters->{'reportableGenes'} } ) {
        my $filename = $RunParameters->{'reportableGenes'}->{$libraryToPrune};
        Annotate->loadLibrary( $filename, $geneListRef, $libraryToPrune, $RunParameters );
        if (exists($RunParameters->{'reported_aberration_names'}->{'ALL'})) {
            if ($RunParameters->{'reported_aberration_names'}->{'ALL'}==1) {
                $geneListRef->{$libraryToPrune}->{'none'}=1;
                $geneListRef->{$libraryToPrune}->{'Unknown'}=1;
                $geneListRef->{$libraryToPrune}->{'ALL'}=1;                
                                
            }
        }
    }
    my $last    = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $counter = 0;
  VCFLINE: while (<SNPEFF>) {
        ++$counter;
        if ( $counter % 100000 == 0 ) {
            print "\t\t+SNPeff Counter at: $counter " . sprintf( "%0.1f", 100000 / ( time - $last ) ) . " per 1 sec\n";
            $last = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
        }
        my $checkInsert = 0;
        my $line        = $_;
        if ( $line =~ /^#/ ) {
            next VCFLINE;
        }
        chomp($line);
        my @temp = split( /\t/, $line );
        my $insert_check = 0;
        if ( $#temp < 4 ) { next VCFLINE }
        my $biomarkerRef     = {};
        my $snpeffVariantRef = { 'toJoinSnpeff' => int(0) };
        my $snpeffMultiplesVariantRef = [];
        my $coord            = join( ":", $temp[0], $temp[1], $temp[3], $temp[4] );
        if ( $#temp == 7 ) {
            my @fields = split( /\;+/, $temp[7] );
            foreach my $info (@fields) {
                if ( $info =~ /(.*?)=(.*)/ ) {
                    my $key = $1;
                    my $val = $2;
                    if ( $key eq "ANN" ) {
                        ($snpeffVariantRef,$snpeffMultiplesVariantRef) = Annotate->addSnpeffInfo( $val, $snpeffVariantRef,$snpeffMultiplesVariantRef, $RunParameters ,$geneListRef);
                        # push(@{$snpeffMultiplesVariantRef},$val);
                    } elsif ($key eq "LOF") {
                        if ($val=~ /\(.*?\|.*?\|.*?\|(.*?)\)/) {
                            $snpeffVariantRef->{'LOF'}=$1;
                        }
                    }elsif ($key eq "NMD") {
                        if ($val =~ /\(.*?\|.*?\|.*?\|(.*?)\)/) {
                             $snpeffVariantRef->{'NMD'}=$1;
                        }
                    } else {
                        $biomarkerRef->{$key} = $val;
                    }
                }
            }
        }
        if (!exists($snpeffVariantRef->{'functionalEffect'} ) && exists($RunParameters->{'snpeffValidAnnotations'}->{'ALL'})) {
            if ($RunParameters->{'snpeffValidAnnotations'}->{'ALL'}==1 && !(exists($snpeffVariantRef->{'functionalEffect'}))) {
                   $snpeffVariantRef->{'functionalEffect'}='Unknown';            
                   $snpeffVariantRef->{'toJoinSnpeff'} = 1;
                   $snpeffVariantRef->{'report'}='comprehensive';
                   $snpeffVariantRef->{'impact'}='low';              
            }
        }
        if ( $snpeffVariantRef->{'toJoinSnpeff'} == 1 && exists( $snpeffVariantRef->{'functionalEffect'} ) ) {
            my $cursor = $conn->find( { 'biomarker' => $coord } );
            $cursor->immortal(1);
            ### First Record
          LOOP: while ( my $doc = $cursor->next ) {
                unless(exists($snpeffVariantRef->{'gene'})) {
                    $snpeffVariantRef->{'gene'}='none';
                }
                my $gene       = $snpeffVariantRef->{'gene'};
                my @ProjectRun = ();
                my @variants   = ();
                if ( exists( $ids{ $doc->{'_id'} } ) ) { next LOOP }
#                $biomarkerRef = dclone($doc);
				$biomarkerRef = $doc;
                $ids{ $doc->{'_id'} } = 1;
                foreach my $variantRef ( @{ $biomarkerRef->{'variants'} } ) {
                    $variantRef->{'gene'}=$gene;
                    if ( exists( $variantRef->{'assay'} ) ) {
                        if ( exists( $geneListRef->{ $variantRef->{'assay'} }->{$gene} ) || exists($geneListRef->{ $variantRef->{'assay'} }->{'ALL'})) {
                            cleanVariants($variantRef);
                            push( @variants, $variantRef );
                            push (@ProjectRun,  $variantRef->{ $RunParameters->{'uidName'} });
                        }
                    }
                }
                $biomarkerRef->{'snpeff'}       = $snpeffVariantRef;
                $biomarkerRef->{'gene'}         = $snpeffVariantRef->{'gene'};
                $biomarkerRef->{'significance'} = "None";
                $biomarkerRef = Annotate->assignAberrationDoc( $biomarkerRef, $RunParameters );
                if ( exists( $biomarkerRef->{'aberration'}->{'effect'} ) ) {
                    $biomarkerRef->{'effect'} = $biomarkerRef->{'aberration'}->{'effect'};
                }
                if ( ( $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "PASS" ||
                       $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "LOWQC" ) && scalar(@variants) > 0 ) {
                    if (exists($RunParameters->{'AnnotateByCoordFields'})) {
                        if ( $RunParameters->{'AnnotateByCoordFields'}->{'ALL'}==0) {
                            $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                     $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                     $RunParameters->{'annotateByCoordCollection'},
                                                                                     $RunParameters->{'AnnotateByCoordFields'}
                                                                                    );
                        } else {
                            $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                     $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                     $RunParameters->{'annotateByCoordCollection'}
                                                                                    );
                        }
                    } else {
                            $biomarkerRef->{'annotate'}          = Annotate->joinOn( 'coord',
                                                                                     $biomarkerRef->{'variants'}->[0]->{'coord'},
                                                                                     $RunParameters->{'annotateByCoordCollection'}
                                                                                    );
                    }
                    my $tempAnnotate;
                    if (exists($RunParameters->{'AnnotateByGeneCodon'})) {
                        if ($RunParameters->{'AnnotateByGeneCodon'}==1 && exists ($biomarkerRef->{'snpeff'}->{'geneCodon'})) {
                            $tempAnnotate          = Annotate->joinOn( 'geneCodon',
                                                                                     $biomarkerRef->{'snpeff'}->{'geneCodon'},
                                                                                     $RunParameters->{'annotateByGeneCodonCollection'}
                                                                                    );
                             if (exists($tempAnnotate->{'geneCodon'}) && exists($tempAnnotate->{'description'})) {
                                  $biomarkerRef->{'annotate'}->{'geneCodon'}=$tempAnnotate->{'geneCodon'};
                                  $biomarkerRef->{'annotate'}->{'codonInfo'}=$tempAnnotate->{'description'};
                             }
                        }
                        $biomarkerRef->{'toAnnotateByGeneCodon'} = int(0);
                    }
                    $biomarkerRef->{'toAnnotateByCoord'} = int(0);
                    if (exists($RunParameters->{'AnnotateByGeneFields'})) {
                        if ( $RunParameters->{'AnnotateByGeneFields'}->{'ALL'}==0) {
                             $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'},$RunParameters->{'AnnotateByGeneFields'}  );
                         } else {
                             $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'}  );
                         }
                    } else {
                             $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'} );
                    }
                    if ( exists( $biomarkerRef->{'geneInfo'}->{'effect'} ) ) {
                        foreach my $effect ( @{ $biomarkerRef->{'geneInfo'}->{'effect'} } ) {
                            if ( $effect->{'variant'} eq $biomarkerRef->{'effect'} ) {
                            }
                        }
                    }
                    if ( exists( $biomarkerRef->{'annotate'}->{'clinvarDiseaseKnownVariant'} ) ) {
                        $biomarkerRef->{'significance'} =~s/^None//;
                        $biomarkerRef->{'significance'} = $biomarkerRef->{'significance'} . ";" . $biomarkerRef->{'annotate'}->{'clinvarDiseaseKnownVariant'};
                    }
                    $biomarkerRef->{'toAnnotateByGene'}            = int(0);
                    $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRun;
                    $biomarkerRef->{'info'}->{'varMax'}                      = int(1);
                    if (exists($RunParameters->{'allSnpeffTranscripts'})) {
                        if ($RunParameters->{'allSnpeffTranscripts'} ==1) {
                            if (scalar(@{$snpeffMultiplesVariantRef})>0) {
                                $biomarkerRef->{'transcripts'}=$snpeffMultiplesVariantRef;
                            }
                        }
                    }
                    if ( ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) > $biomarkerRef->{'info'}->{'varMax'} && ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) < $RunParameters->{'varGlobalMax'} ) {
                        $biomarkerRef->{'info'}->{'varMax'} = int( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 4 );
                    }
                    elsif ( ( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) >= $RunParameters->{'varGlobalMax'} ) {
                        $biomarkerRef->{'info'}->{'varMax'} = int( $RunParameters->{'varGlobalMax'} );
                    }
                    if ( $biomarkerRef->{'info'}->{'varArray'} >= $biomarkerRef->{'info'}->{'varMax'} ) {
                        $biomarkerRef->{'info'}->{'varFilled'} = int(1);
                    } else {
                        $biomarkerRef->{'info'}->{'varFilled'} = int(0);
                    }
                    $biomarkerRef->{'variants'}                    = \@variants;
                    $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRun;
                    $biomarkerRef                                  = Annotate->isCommon($biomarkerRef);
                    if (exists($biomarkerRef->{'snpeff'}->{'report'})) {
                        $biomarkerRef->{'report'}=$biomarkerRef->{'snpeff'}->{'report'};
                        $biomarkerRef->{'impact'}=$biomarkerRef->{'snpeff'}->{'impact'};
                    } else {
                        $biomarkerRef->{'report'}='comprehensive';
                        $biomarkerRef->{'impact'}='low';
                    }
                    delete( $biomarkerRef->{'_id'} );
                    $connStage->insert( $biomarkerRef, { safe => 1 } );
                }
            }
        }
    }
    $conn->remove();
    if ( -e ("$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf") ) { `rm $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.an.tmp.vcf` }
    if ( -e ("$RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf") ) { `rm $RunParameters->{'runDir'}/snpeff.$collectionDb.pid$$.tmp.vcf` }
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t-Completed annotateSnpeff: $diff seconds\n";
    return 1;
}

sub annotateGene {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    print "\t+Annotating Genes\n";
    my $isGene    = $RunParameters->{'isGene'};
    my $start     = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $conn      = $RunParameters->{'vcfCollection'}->{ 'geneInsert_' . $collectionDb };
    my $connStage = $RunParameters->{'vcfCollection'}->{ 'staged_' . $collectionDb };
    #$conn->ensure_index( { 'biomarker' => 1 } );
    #$conn->ensure_index( { 'group' => 1 } );
    my $geneListRef = {};
    foreach my $libraryToPrune ( keys %{ $RunParameters->{'reportableGenes'} } ) {
        my $filename = $RunParameters->{'reportableGenes'}->{$libraryToPrune};
        Annotate->loadLibrary( $filename, $geneListRef, $libraryToPrune, $RunParameters );
        if (exists($RunParameters->{'reported_aberration_names'}->{'ALL'})) {
            if ($RunParameters->{'reported_aberration_names'}->{'ALL'}==1) {
                $geneListRef->{$libraryToPrune}->{'none'}=1;
                $geneListRef->{$libraryToPrune}->{'ALL'}=1;
                $geneListRef->{$libraryToPrune}->{'Unknown'}=1;
                
            }
        }
    }
    my $cursor = $conn->find();
    $cursor->immortal(1);
    my $last            = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $counter         = 0;
    my %biomarkersFound = ();
  LOOP: while ( my $doc = $cursor->next ) {
        ++$counter;
        if ( $counter % 100000 == 0 ) { print "\t\t-Gene Annotate Counter at: $counter " . sprintf( "%0.3f", 100 / ( time - $last ) ) . " 1sec\n"; $last = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec(); }
        my $id           = $doc->{'_id'};
        my $biomarkerRef = $doc;
        my @genes        = ();
        my @biomarkers   = ();
        my @ProjectRuns  = ();
        my @coords=();
        my $ProjectRun   = $biomarkerRef->{'variants'}->[0]->{ $RunParameters->{'uidName'} };
        push (@ProjectRuns,$ProjectRun);
        my $uid = "";
        unless (exists($biomarkerRef->{'report'})) {
            $biomarkerRef->{'report'}='standard';
        } elsif ($biomarkerRef->{'report'} eq "comprehensive") {
            $biomarkerRef->{'report'}='standard';
        }
        ##  Structural Variant ###
        if ( $biomarkerRef->{'type'} eq 'StructuralVariant' ) {
            my $chr1 = $biomarkerRef->{'variants'}->[0]->{'chr'};
            $chr1 =~ s/chr//;
            my $pos1    = 100 * int( $biomarkerRef->{'variants'}->[0]->{'hg19pos'} / 100 );
            my $chr2    = "None";
            my $pos2    = "None";
             my @coordSet=();
            my @geneSet = ();
            if ( exists( $biomarkerRef->{'variants'}->[0]->{'END'} ) ) {
                if ( $biomarkerRef->{'variants'}->[0]->{'END'} =~ /(.*?):(.*)/ ) {
                    $chr2 = $1;
                    $pos2 = $2;
                    $chr2 =~ s/chr//;
                    $pos2 = 100 * int( $2 / 100 );
                }elsif (exists($biomarkerRef->{'variants'}->[0]->{'END'})) {
                    $chr2 = $biomarkerRef->{'variants'}->[0]->{'CHR2'};
                    $pos2 = 100 * int( $biomarkerRef->{'variants'}->[0]->{'END'} / 100 );
                }
                else {
                    $chr2 = $chr1;
                    $pos2 = 100 * int( $biomarkerRef->{'variants'}->[0]->{'END'} / 100 );
                }
            }
            elsif ( exists( $biomarkerRef->{'variants'}->[0]->{'SV'} ) ) {
                my $tVar = $biomarkerRef->{'variants'}->[0]->{'SV'};
                if ( $tVar =~ /.:(.*?):(.*?)\|.:(.*?):(.*)/ ) {
                    $chr1 = $1;
                    $pos1 = 100 * int( $2 / 100 );
                    $chr2 = $3;
                    $pos2 = 100 * int( $4 / 100 );
                }
            }
            my $sv = "$chr1:$pos1-$chr2:$pos2";
            foreach my $val ( "chr$chr1:$pos1", "chr$chr2:$pos2" ) {
                if (exists($RunParameters->{'reported_aberration_names'}->{'ALL'})) {
                    if ($RunParameters->{'reported_aberration_names'}->{'ALL'}==1 && !(exists ($isGene->{$val} )) ) {
                        $isGene->{$val}='none';
                    }
                }
                if ( exists( $isGene->{$val} ) ) {
                    my @temp = split( /;+/, $isGene->{$val} );
                    foreach my $gene (@temp) {
                        if ( (length($gene) > 0 && exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{$gene} ))
                            || exists($geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{'ALL'}) ) {
                            my $ext="break";
                            if ($gene eq "none") {
                                my $chr=$biomarkerRef->{'variants'}->[0]->{'chrPos'};
                               $chr=~s/:.*//g;
                               my $pos=$biomarkerRef->{'variants'}->[0]->{'chrPos'};
                               $pos=~s/.*://g;
                               $pos=int($pos/1000000)*1000000;
                               $ext="break($chr:$pos)";
                            }
                            $uid = $gene . "_" . $ext . "-" .  $biomarkerRef->{'variants'}->[0]->{'filepath'} . "-" . $biomarkerRef->{'variants'}->[0]->{'varName'};
                            unless ( exists( $biomarkersFound{$uid} ) ) {
                               cleanVariants($biomarkerRef->{'variants'}->[0]);
                                push @genes,      $gene;
                                push @biomarkers, $gene . "_" . $ext;
                                push @coords,$val;
                                $biomarkersFound{$uid} = 1;
                                if ( scalar(@temp) == 2 ) {
                                    push @geneSet, $gene;
                                    push @coordSet,$val;
                                }
                            }
                        }
                    }
                }
            }
            if ( scalar(@geneSet) == 2 ) {
                $uid = "$geneSet[0]_$geneSet[1]-$biomarkerRef->{'variants'}->[0]->{'filepath'}";
                unless ( exists( $biomarkersFound{$uid} ) ) {
                      push (@genes,"$geneSet[0]_$geneSet[1]");
                      push @coords,"$coordSet[0]_$coordSet[1]";
                      push (@genes,"$geneSet[1]_$geneSet[0]");
                      push @coords,"$coordSet[1]|_coordSet[0]";
                      push (@biomarkers, "$geneSet[0]_$geneSet[1]");
                      push (@biomarkers, "$geneSet[1]_$geneSet[0]");
                    $biomarkersFound{$uid} = 1;
                    $biomarkersFound{'$geneSet[1]' . "|". $geneSet[0] . '-' . $ProjectRun} = 1;
                }
            }
            ##  Copy Number ###
        }elsif ($biomarkerRef->{'type'} eq 'FocalDeletion'
            || $biomarkerRef->{'type'} eq 'FocalGain'
            || $biomarkerRef->{'type'} eq 'LargeInsertion'
            || $biomarkerRef->{'type'} eq 'LargeDeletion'
            || $biomarkerRef->{'type'} eq 'LargeBlock'
            || $biomarkerRef->{'type'} eq 'NovelInsertion'
            || $biomarkerRef->{'type'} eq 'Inversion'
            || $biomarkerRef->{'type'} eq 'CnvChange'
            || $biomarkerRef->{'type'} eq 'TandemDuplication'
            || $biomarkerRef->{'type'} eq 'MobileElementionDeletion'
            || $biomarkerRef->{'type'} eq 'MobileElementionInsertion'
            || $biomarkerRef->{'type'} eq 'LargeBlockSubstitution'
            || $biomarkerRef->{'type'} eq 'variant'
         )
        {
            my $chr = $biomarkerRef->{'variants'}->[0]->{'chr'};
            $chr =~ s/chr//;
            my $pos = 100 * int( $biomarkerRef->{'variants'}->[0]->{'hg19pos'} / 100 );
            my $end = 100 * int( $biomarkerRef->{'variants'}->[0]->{'hg19pos'} / 100 );
            if ( exists( $biomarkerRef->{'variants'}->[0]->{'END'} ) ) {
                $end = 100 * int( $biomarkerRef->{'variants'}->[0]->{'END'} / 100 );
                if ( $biomarkerRef->{'variants'}->[0]->{'END'} =~ /(.*?):(.*)/ ) {
                    if ( length($2) > 0 ) { $end = 100 * int( $2 / 100 ); }
                }
            }
            if ( abs( $pos - $end ) < 500000000 ) {
                for ( my $i = $pos - 100 ; $i <= $end + 100 ; $i = $i + 100 ) {
                    if ( exists( $isGene->{'chr' . $chr . ':' . $i} ) ) {
                        my @temp = split( /;+/, $isGene->{'chr' . $chr . ':' . $i} );
                        foreach my $gene (@temp) {
                            if ( (length($gene) > 0 && exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{$gene} )) 
                             || exists($geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{'ALL'}) ) {
                                if ( $biomarkerRef->{'type'} eq 'FocalGain' || $biomarkerRef->{'type'} eq 'LargeBlockSubstitution' || $biomarkerRef->{'type'} eq 'LargeInsertion') {
                                    $uid = $gene . '_gain-' . $biomarkerRef->{'variants'}->[0]->{'filepath'}. "-" . $biomarkerRef->{'variants'}->[0]->{'varName'};
                                    unless ( exists( $biomarkersFound{$uid} ) ) {
                                        push( @genes,      $gene );
                                        push( @biomarkers, $gene . "_gain" );
                                        push (@coords,"chr$chr:$i");
                                        $biomarkersFound{$uid} = 1;
                                    }
                                } elsif ($biomarkerRef->{'type'} eq 'FocalDeletion' || $biomarkerRef->{'type'} eq 'LargeDeletion') {
                                    $uid = $gene . '_loss-' . $biomarkerRef->{'variants'}->[0]->{'filepath'} . "-" . $biomarkerRef->{'variants'}->[0]->{'varName'};
                                    unless ( exists( $biomarkersFound{$uid} ) ) {
                                        push( @genes,      $gene );
                                        push( @biomarkers, $gene . "_loss" );
                                        push (@coords,"chr$chr:$i");
                                        $biomarkersFound{$uid} = 1;
                                    }
                                } else {
                                    $uid = $gene . '_' . $biomarkerRef->{'type'} . '-' . $biomarkerRef->{'variants'}->[0]->{'filepath'} . "-" . $biomarkerRef->{'variants'}->[0]->{'varName'};
                                    unless ( exists( $biomarkersFound{$uid} ) ) {
                                        push( @genes,      $gene );
                                        push( @biomarkers, $gene . "_" . $biomarkerRef->{'type'} );
                                        push (@coords,"chr$chr:$i");
                                        $biomarkersFound{$uid} = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        elsif ( $biomarkerRef->{'type'} eq 'FusedGenes' ) {
            if ( exists( $biomarkerRef->{'variants'}->[0]->{'FUSED_GENE'} ) ) {
                my @fusedGenes = split( /_/, $biomarkerRef->{'variants'}->[0]->{'FUSED_GENE'} );
                if ( scalar(@fusedGenes) == 2 ) {
                    if (   exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ $fusedGenes[0] } )
                        || exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ $fusedGenes[1] } ) 
                        || exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ "ALL" }))
                    {
                        $uid = "$fusedGenes[0]_$fusedGenes[1]" . '-' . $biomarkerRef->{'variants'}->[0]->{'filepath'};
                        unless ( exists( $biomarkersFound{$uid} ) || $fusedGenes[0] =~ /ENS/ || $fusedGenes[1] =~ /ENS/ ) {
                            push( @genes,      "$fusedGenes[0]_$fusedGenes[1]" );
                            push( @biomarkers, "$fusedGenes[0]_$fusedGenes[1]" );
                            $biomarkersFound{$uid} = 1;
                        }
                    }
                }
            }
        }
        elsif ( $biomarkerRef->{'type'} eq 'UnderExpressed' || $biomarkerRef->{'type'} eq 'OverExpressed' ) {
            if ( exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ $biomarkerRef->{'variants'}->[0]->{'gene'} } ) ||
                exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ "ALL" } ) ) {
                my $gene = $biomarkerRef->{'variants'}->[0]->{'gene'};
                if ( $biomarkerRef->{'type'} eq 'UnderExpressed' ) {
                    $uid = $gene . '_UnderExpressed-' . $biomarkerRef->{'variants'}->[0]->{'filepath'};
                    unless ( exists( $biomarkersFound{$uid} ) ) {
                        push( @genes,      $biomarkerRef->{'variants'}->[0]->{'gene'} );
                        push( @biomarkers, $biomarkerRef->{'variants'}->[0]->{'gene'} . "_UnderExpressed" );
                        push (@coords,$biomarkerRef->{'variants'}->[0]->{'gene'});
                        $biomarkersFound{$uid} = 1;
                    }
                }
                else {
                    $uid = $gene . '_OverExpressed-' . $biomarkerRef->{'variants'}->[0]->{'filepath'};
                    unless ( exists( $biomarkersFound{$uid} ) ) {
                        push( @genes,      $biomarkerRef->{'variants'}->[0]->{'gene'} );
                        push( @biomarkers, $biomarkerRef->{'variants'}->[0]->{'gene'} . "_OverExpressed" );
                        push (@coords,$biomarkerRef->{'variants'}->[0]->{'gene'});
                        $biomarkersFound{$uid} = 1;
                    }
                }
            }
        }
        elsif ( $biomarkerRef->{'type'} eq 'TranscriptVariant' ) {
            if ( exists( $biomarkerRef->{'variants'}->[0]->{'TRANSCRIPT'} ) ) {
                unless ( exists( $biomarkersFound{$uid} ) ) {
                    $uid = $biomarkerRef->{variants}->[0]->{'TRANSCRIPT'} . '_Variant' . $biomarkerRef->{'variants'}->[0]->{'filepath'} ;
                    push( @genes,      $biomarkerRef->{'variants'}->[0]->{'gene'} );
                    push( @biomarkers, $biomarkerRef->{'variants'}->[0]->{'TRANSCRIPT'} );
                    push (@coords,$biomarkerRef->{'variants'}->[0]->{'gene'});
                    $biomarkersFound{$uid} = 1;
                }
            }
        } elsif ($biomarkerRef->{'type'} eq 'GeneEvent' ) {
           if ( exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ $biomarkerRef->{'variants'}->[0]->{'gene'} } ) 
                ||  exists( $geneListRef->{ $biomarkerRef->{'variants'}->[0]->{'assay'} }->{ "ALL" } ) ) {
                my $gene = $biomarkerRef->{'variants'}->[0]->{'gene'};
                $uid = $gene . '_GeneEvent-' . $biomarkerRef->{'variants'}->[0]->{'filepath'};
                unless ( exists( $biomarkersFound{$uid} ) ) {
                    push( @genes,      $biomarkerRef->{'variants'}->[0]->{'gene'} );
                    push( @biomarkers, $biomarkerRef->{'variants'}->[0]->{'gene'} . "_GeneEvent" );
                    push (@coords,$biomarkerRef->{'variants'}->[0]->{'gene'});
                    $biomarkersFound{$uid} = 1;
                }
            }
        } elsif ( $biomarkerRef->{'type'} eq 'Unknown' ) {
        }
        $biomarkerRef->{'toAnnotateByCoord'} = int(0);
        $biomarkerRef->{'toAnnotateByChrPos'} = int(0);
        $biomarkerRef->{'toAnnotateByGene'}  = int(0);
        delete( $biomarkerRef->{'_id'} );
        my $geneCounter = 0;
        while ( my $gene = shift @genes ) {
            my $biomarker = shift @biomarkers;
            my $lcoord=shift @coords;
            $biomarkerRef->{'gene'}      = $gene;
            $biomarkerRef->{'variants'}->[0]->{'gene'}      = $gene;
            $biomarkerRef->{'biomarker'} = $biomarker;
            $biomarkerRef->{'variants'}->[0]->{'biomarker'} = $biomarker;
            if (exists($RunParameters->{'AnnotateByGeneFields'})) {
                if ( $RunParameters->{'AnnotateByGeneFields'}->{'ALL'}==0) {
                     $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'},$RunParameters->{'AnnotateByGeneFields'}  );
                 } else {
                     $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'}  );
                 }
            } else {
                     $biomarkerRef->{'geneInfo'}          = Annotate->joinOn( 'gene', $biomarkerRef->{'gene'}, $RunParameters->{'annotateByGeneCollection'} );
            }
            $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRuns;
            $biomarkerRef = Annotate->assignAberrationDoc( $biomarkerRef, $RunParameters );
            if (defined($lcoord)) {
                $biomarkerRef->{'aberration'}->{'aberration_coord'}=$lcoord;
            }
            Maintain->assignGroup( $biomarkerRef, $RunParameters );
            if ($biomarkerRef->{'gene'} ne "none" && $biomarkerRef->{'gene'} ne "Unknown") {
                $biomarkerRef->{'report'}='standard';
                $biomarkerRef->{'impact'}='significant';
            } else {
                $biomarkerRef->{'report'}='comprehensive';
                $biomarkerRef->{'impact'}='low';
            }
            $connStage->insert( $biomarkerRef, { safe => 1 } );
            ++$geneCounter;
        }
        if ( $geneCounter == 0 && exists( $RunParameters->{'allBiomarkers'} ) ) {
            $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRuns;
            $biomarkerRef = Annotate->assignAberrationDoc( $biomarkerRef, $RunParameters );
            $biomarkerRef->{'biomarkerProjectRun'} = $biomarkerRef->{'biomarker'} . ":" . $biomarkerRef->{'variants'}->[0]->{ $RunParameters->{'uidName'} };    #
            Maintain->assignGroup( $biomarkerRef, $RunParameters );
            $connStage->insert( $biomarkerRef, { safe => 1 } );
        }
    }
    $conn->remove();
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t-Completed annotateGene: $diff seconds\n";
    return 1;
}

sub cleanVariants {
    my $variantRef=shift;
    if ($variantRef->{'varName'} eq ".") {
        undef $variantRef->{'varName'};
        delete $variantRef->{'varName'};
    }
    unless (exists($variantRef->{'genotype'}->{'s1'})) {
        undef $variantRef->{'genotype'};
        delete $variantRef->{'genotype'};
    }
    undef $variantRef->{'valid'};
    delete $variantRef->{'valid'};
    undef $variantRef->{'status'};
    delete $variantRef->{'status'};
}

sub finalize {
    my $function      = shift;
    my $RunParameters = shift;
    my $collectionDb = shift;
    print "\t+Finalizing Variants in Markers Database\n";
    my $start = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec();
    my $conn  = $RunParameters->{'vcfCollection'}->{ 'staged_' . $collectionDb };
    $conn->ensure_index( { 'info.group' => 1 } );
    my %ids    = ();
    my $cursor = $conn->find();
    $cursor->immortal(1);
    my $c = 0;
  LOOP: while ( my $doc = $cursor->next ) {
        ++$c;
        if ( $c % 10000 == 0 ) { print "\t\t+Count for finalizing is at $c\n"; }
        my @ProjectRun = ();
        my @variants   = ();
        unless ( exists( $doc->{'_id'} ) )         { next LOOP }
        if     ( exists( $ids{ $doc->{'_id'} } ) ) { next LOOP }
        $ids{ $doc->{'_id'} } = 1;        
#        my $biomarkerRef = dclone($doc);
        my $biomarkerRef = $doc;
        foreach my $variantRef ( @{ $biomarkerRef->{'variants'} } ) {
            my $includeGene = 1;
            if ( exists( $RunParameters->{'filterOutGenes'} ) ) {
                $includeGene = Annotate->filterOutGenes( $variantRef, $RunParameters, $biomarkerRef->{'gene'} );
            }
            if ( $includeGene eq "1" ) {
                push( @variants, $variantRef );
                push( @ProjectRun, $variantRef->{ $RunParameters->{'uidName'} });
            }
        }
        my $cursorGroup;
        Maintain->assignGroup( $biomarkerRef, $RunParameters );
        if ( exists($RunParameters->{'multipledVariant'})) {
            $cursorGroup = $conn->find( { 'info.group' => $biomarkerRef->{'info'}->{'group'} } );
            $cursorGroup->immortal(1);
          LOOP2: while ( my $newDoc = $cursorGroup->next ) {
                if ( exists( $ids{ $newDoc->{'_id'} } ) ) { next LOOP2 }
                $ids{ $newDoc->{'_id'} } = 1;
                my $biomarkerRef2 = $newDoc;                
                foreach my $variantRef ( @{ $biomarkerRef2->{'variants'} } ) {
                    my $includeGene = 1;
                    if ( exists( $RunParameters->{'filterOutGenes'} ) ) {
                        $includeGene = Annotate->filterOutGenes( $variantRef, $RunParameters, $biomarkerRef->{'gene'} );
                    }
                    if ( $includeGene eq "1" ) {
                        push( @variants, $variantRef );
                        push ( @ProjectRun,$variantRef->{ $RunParameters->{'uidName'} } );
                        $biomarkerRef->{'count'}->{'all'}                           = $biomarkerRef->{'count'}->{'all'} + $biomarkerRef2->{'count'}->{'all'};
                        $biomarkerRef->{'count'}->{'hets'}                          = $biomarkerRef->{'count'}->{'hets'} + $biomarkerRef2->{'count'}->{'hets'};
                        $biomarkerRef->{'count'}->{'homs'}                          = $biomarkerRef->{'count'}->{'homs'} + $biomarkerRef2->{'count'}->{'homs'};
                        $biomarkerRef->{'count'}->{'oths'}                          = $biomarkerRef->{'count'}->{'oths'} + $biomarkerRef2->{'count'}->{'oths'};
                        $biomarkerRef->{'count'}->{'refs'}                          = $biomarkerRef->{'count'}->{'refs'} + $biomarkerRef2->{'count'}->{'refs'};
                        ++$biomarkerRef->{'info'}->{'varArray'};
                    }
                }
            }
        }
        $biomarkerRef->{ $RunParameters->{'UidName'} } = \@ProjectRun;
        $biomarkerRef->{'variants'} = \@variants;
        delete( $biomarkerRef->{'_id'} );
        delete( $biomarkerRef->{'toAnnotateByGene'} );
        delete( $biomarkerRef->{'toAnnotateByGeneCodon'} );
        delete( $biomarkerRef->{'toAnnotateByCoord'} );
        delete( $biomarkerRef->{'toAnnotateByChrPos'} );
        $biomarkerRef->{'info'}->{'varMax'} = int(1);
        if ( $biomarkerRef->{'info'}->{'varArray'} >= $biomarkerRef->{'info'}->{'varMax'} && ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 < $RunParameters->{'varGlobalMax'} ) {
            $biomarkerRef->{'info'}->{'varMax'} = int( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 4 );
        } elsif ( int( ( $biomarkerRef->{'info'}->{'varArray'} + 1 ) * 2 ) >= $RunParameters->{'varGlobalMax'} ) {
            $biomarkerRef->{'info'}->{'varMax'} = int( $RunParameters->{'varGlobalMax'} );
        }
        if ( $biomarkerRef->{'info'}->{'varArray'} >= $biomarkerRef->{'info'}->{'varMax'} ) {
            $biomarkerRef->{'info'}->{'varFilled'} = int(1);
        } else {
            $biomarkerRef->{'info'}->{'varFilled'} = int(0);
        }
        unless (exists ($RunParameters->{'multipledVariant'})) {
            $biomarkerRef->{'info'}->{'varFilled'} = int(1);
        }
        forceBiomarkerType( 'geneInfo', $biomarkerRef, $RunParameters );
        forceBiomarkerType( 'annotate', $biomarkerRef, $RunParameters );
        $biomarkerRef->{'dbFreq'} = 0;
        MongoDB::force_double($biomarkerRef->{'dbFreq'});
        if (exists($biomarkerRef->{'geneInfo'}->{'geneId'})) {
            $biomarkerRef->{'geneId'}=$biomarkerRef->{'geneInfo'}->{'geneId'};
        }
        if (exists($biomarkerRef->{'snpeff'}->{'geneId'})) {
            $biomarkerRef->{'geneId'}=$biomarkerRef->{'snpeff'}->{'geneId'};
        }
        $biomarkerRef->{'dbFreqCat'} = 'rare';
        $biomarkerRef->{'biomarkerCount'} = $biomarkerRef->{'info'}->{'varArray'};
        if (exists($RunParameters->{'FilterAnnotatePipeline'})) {
             if (  $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "PASS" ||
                   $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "LOWQC"
                 ) {
                   Annotate->FilterAnnotatePipeline($biomarkerRef,$RunParameters);
             }
        }
        if (exists($RunParameters->{'FilterGeneInfoPipeline'})) {
             if ( $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "PASS" ||
                   $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "LOWQC"
                 ) {
                 Annotate->FilterGeneInfoPipeline($biomarkerRef,$RunParameters);
             }
        }
        if ( ( $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "PASS" ||
               $biomarkerRef->{'aberration'}->{'aberration_filter'} eq "LOWQC"
             ) && scalar(@variants) > 0 && exists( $biomarkerRef->{'collectionDb'} )
            ) {
#  my $num=  total_size($biomarkerRef); print "$biomarkerRef->{'biomarker'}:$num\n";
               $RunParameters->{$collectionDb}->insert( $biomarkerRef, { safe => 1 } );
        }
    }
    $conn->remove();
    my $diff = localtime->hour * 3600 + localtime->min() * 60 + localtime->sec() - $start;
    print "\t\t-Completed finalized: $diff seconds\n";
    return 1;
}

sub FilterGeneInfoPipeline {
    my $function      = shift;
    my $biomarkerRef    = shift;
    my $RunParameters = shift;
    my $collectionDb=$biomarkerRef->{'collectionDb'};
    my $newDatabaseFilter = $biomarkerRef->{'aberration'}->{'aberration_filter'};
    my $newVariantFilter = $biomarkerRef->{'aberration'}->{'aberration_filter'};
    if ( exists( $RunParameters->{'FilterGeneInfoPipeline'} ) ) {
        LOOP: for ( my $i = 0 ; $i < scalar( @{ $RunParameters->{'FilterGeneInfoPipeline'} } ) ; ++$i ) {
             my $state = 'none';
             foreach my $message ( sort keys %{ $RunParameters->{'FilterGeneInfoPipeline'}->[$i] } ) {
                 if ( $message =~ /SetBiomarkerTo(.*)/ ) {
                     $state = $1;
                     ( $newDatabaseFilter, $newVariantFilter ) = Annotate->SetBiomarkerTo( $state, $biomarkerRef, $RunParameters, $i, $newDatabaseFilter, $newVariantFilter,'GeneInfo','geneInfo');
                 } else {
                     print "\t\t+Warning, we shouldn't be here.  Did you set up FilterGeneInfoPipeline with SetBiomarkerToPASS?\n";
                 }
             }
        }
    }
    $biomarkerRef->{'aberration'}->{'aberration_filter'}=$newDatabaseFilter ;
}

sub FilterAnnotatePipeline {
    my $function      = shift;
    my $biomarkerRef    = shift;
    my $RunParameters = shift;
    my $collectionDb=$biomarkerRef->{'collectionDb'};
    unless ( exists( $biomarkerRef->{$collectionDb} ) ) { return }
    my $newDatabaseFilter = $biomarkerRef->{'aberration'}->{'aberration_filter'};
    my $newVariantFilter = $biomarkerRef->{'aberration'}->{'aberration_filter'};
    if ( exists( $RunParameters->{'FilterAnnotatePipeline'} ) ) {
        LOOP: for ( my $i = 0 ; $i < scalar( @{ $RunParameters->{'FilterAnnotatePipeline'} } ) ; ++$i ) {
              my $state = 'none';
              foreach my $message ( sort keys %{ $RunParameters->{'FilterAnnotatePipeline'}->[$i] } ) {
                  if ( $message =~ /SetBiomarkerTo(.*)/ ) {
                      $state = $1;
                      ( $newDatabaseFilter, $newVariantFilter ) = Annotate->SetBiomarkerTo( $state, $biomarkerRef, $RunParameters, $i, $newDatabaseFilter, $newVariantFilter,'Annotate','annotate' );
                  } else {
                      print "\t\t+Warning, we shouldn't be here.  Did you set up FilterAnnotatePipeline with SetBiomarkerToPASS?\n";
                  }
              }
         }
     }
     $biomarkerRef->{'aberration'}->{'aberration_filter'}=$newDatabaseFilter ;
}

sub SetBiomarkerTo {
    my $function          = shift;
    my $state             = shift;
    my $biomarkerRef        = shift;
    my $RunParameters     = shift;
    my $i                 = shift;
    my $newDatabaseFilter = shift;
    my $newVariantFilter  = shift;
    my $Type=shift;
    my $type=shift;
    if ( length($state) == 0 ) { return ( $newDatabaseFilter, $newVariantFilter ) }
    my $change = "SetBiomarkerTo" . $state;
    my $collectionDb = $biomarkerRef->{'collectionDb'};
    unless ( exists( $RunParameters->{'Filter' . $Type . 'Pipeline'}->[$i]->{$change}->{$collectionDb} ) ) {
        return ( $newDatabaseFilter, $newVariantFilter );
    }
    my $dbState = $RunParameters->{'Filter' . $Type . 'Pipeline'}->[$i]->{$change}->{'DatabaseIfFilter'};
    for ( my $j = 0 ; $j < scalar( @{ $RunParameters->{'Filter' . $Type . 'Pipeline'}->[$i]->{$change}->{$collectionDb} } ) ; ++$j ) {
         my $ref=$RunParameters->{'Filter' . $Type . 'Pipeline'}->[$i]->{$change}->{$collectionDb}->[$j];
         my $filter=$RunParameters->{'Filter' . $Type . 'Pipeline'}->[$i]->{$change}->{$collectionDb}->[$j]->{'field'};
         if ( exists( $biomarkerRef->{$type}->{$filter} )) {
              if ( $newVariantFilter ne $state ) { #success,fail
                  $newVariantFilter = Maintain->checkConditions( $biomarkerRef->{$type}, $filter, $ref, $state, $newVariantFilter );
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

sub assignAberrationDocForSnvs {
    my $function      = shift;
    my $reportRef     = shift;
    my $biomarkerRef  = shift;
    my $RunParameters = shift;
    $reportRef->{'aberration_description'}   = $biomarkerRef->{'snpeff'}->{'functionalEffect'};
    $reportRef->{'aberration_description'}=~s/_/ /g;
    $biomarkerRef->{'functionalEffect'}   = "Unknown";
    if (exists($biomarkerRef->{'snpeff'}->{'functionalEffect'})) {
        $biomarkerRef->{'functionalEffect'}=$biomarkerRef->{'snpeff'}->{'functionalEffect'};
    }
    if (exists($biomarkerRef->{'snpeff'}->{'functionalImpact'})) {
        if ($biomarkerRef->{'snpeff'}->{'functionalImpact'} eq "HIGH" || $biomarkerRef->{'snpeff'}->{'functionalImpact'} eq "MODERATE") {
            $reportRef->{'aberration_name'}   = "High Or Moderate Variant";
            $biomarkerRef->{'significance'}   = "Yes";
            $reportRef->{'aberration_effect'} = "GeneLost";
            $reportRef->{'aberration_type'}   = 'SNV';
            $reportRef->{'biomarker_type'}    = "gene";
            $reportRef->{'aberration_value'}  = $biomarkerRef->{'snpeff'}->{'functionalEffect'};
            $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
            $reportRef->{'effect'}            = $reportRef->{'gene'} . "(" . $biomarkerRef->{'snpeff'}->{'functionalEffect'} . ")";
            $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
        }
    }
    if (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/exon_loss_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/frameshift_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/stop_gained/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/stop_lost/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/rare_amino_acid_variant/
          ) {
        $reportRef->{'aberration_name'}   = "Loss Of Function";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $reportRef->{'aberration_type'}   = 'CNV_LOSS';
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'snpeff'}->{'AminoAcidChange'};
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(" . $biomarkerRef->{'snpeff'}->{'AminoAcidChange'} . ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/splice_acceptor_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/splice_donor_variant/
        ) {
        $reportRef->{'aberration_name'}   = "Splice-site Loss";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $reportRef->{'aberration_type'}   = 'CNV_LOSS';
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'snpeff'}->{'hgvs'};
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(Splicing Acceptor/Donor Loss)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'snpeff'}->{'functionalEffect'}=~/missense_variant/) {
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'snpeff'}->{'AminoAcidChange'};
        $reportRef->{'aberration_type'}   = 'SNV';
        $reportRef->{'aberration_name'}   = "Missense";
        $reportRef->{'aberration_effect'} = "AminoAcidChange";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(" . $biomarkerRef->{'snpeff'}->{'AminoAcidChange'} . ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_value'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/inframe_insertion/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/disruptive_inframe_insertion/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/inframe_deletion/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/disruptive_inframe_deletion/
          ) {
        $reportRef->{'aberration_effect'} = "GeneAltered";
        $reportRef->{'aberration_type'}   = 'INDEL';
        $reportRef->{'aberration_name'}   = "Inframe Indel";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'snpeff'}->{'AminoAcidChange'};
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(" . $biomarkerRef->{'snpeff'}->{'AminoAcidChange'} . ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    } elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/splice_branch_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/splice_region_variant/
          ) {
        $reportRef->{'aberration_effect'} = "SplicingAltered";
        $reportRef->{'aberration_name'}   = "Splicing Altered";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_type'}   = 'SPLICE';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(Splicing Altered)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'snpeff'}->{'functionalEffect'} eq "synonymous_variant" ) {
        $reportRef->{'aberration_value'} = $biomarkerRef->{'snpeff'}->{'AminoAcidChange'};
        $reportRef->{'aberration_type'}   = 'SNV';
        $reportRef->{'aberration_name'}   = "Synonymous";
        $reportRef->{'aberration_effect'} = "AminoAcidChange";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(" . $biomarkerRef->{'snpeff'}->{'AminoAcidChange'} . ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_value'} );
    } elsif ($biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/initiator_codon_variant/)       {

        $reportRef->{'aberration_effect'} = "LowImpact";
        $reportRef->{'aberration_name'}   = "Initiator Codon Variant";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_type'}   = 'CODON';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(initiator codon variant)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }  elsif ($biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/stop_retained_variant/) {
        $reportRef->{'aberration_effect'} = "StopRetainedVariant";
        $reportRef->{'aberration_name'}   = "Stop Retained Variant";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_type'}   = 'CODON';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(stop retained variant)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    } elsif ($biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/coding_sequence_variant/) {
        $reportRef->{'aberration_effect'} = "CodingSequenceVariant";
        $reportRef->{'aberration_name'}   = "Coding Sequence Variant";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'aberration_type'}   = 'CODON';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(coding sequence variant)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );

    } elsif ($biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/UTR.*truncation/) {
        $reportRef->{'aberration_effect'} = "UTR";
        $reportRef->{'aberration_name'}   = "UTR Truncation";
        $reportRef->{'biomarker_type'}    = "UTR";
        $reportRef->{'aberration_type'}   = 'UTR';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(UTR)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/5_prime_UTR_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/3_prime_UTR_variant/
           ) {
        $reportRef->{'aberration_effect'} = "UTR";
        $reportRef->{'aberration_name'}   = "UTR";
        $reportRef->{'biomarker_type'}    = "UTR";
        $reportRef->{'aberration_type'}   = 'UTR';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (UTR)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/upstream_gene_variant/ ||
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/downstream_gene_variant/
           ) {
        $reportRef->{'aberration_effect'} = "Upstream-Downstream";
        $reportRef->{'aberration_name'}   = "Upstream-Downstream";
        $reportRef->{'biomarker_type'}    = "Upstream-Downstream";
        $reportRef->{'aberration_type'}   = 'Upstream-Downstream';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Upstream/Downstream)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/TF_binding_site_variant/
           ) {
        $reportRef->{'aberration_effect'} = "TFBindingSite";
        $reportRef->{'aberration_name'}   = "TF Binding Site";
        $reportRef->{'biomarker_type'}    = "TFBindingSite";
        $reportRef->{'aberration_type'}   = 'TFBindingSite';
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(". $reportRef->{'aberration_name'}. ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/miRNA/
           ) {
        $reportRef->{'aberration_effect'} = "miRNA";
        $reportRef->{'aberration_name'}   = "miRNA";
        $reportRef->{'biomarker_type'}    = "miRNA";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(". $reportRef->{'aberration_name'}. ")";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/sequence_feature/
           ) {
        $reportRef->{'aberration_effect'} = "sequence_feature";
        $reportRef->{'aberration_name'}   = "Sequence Feature";
        $reportRef->{'biomarker_type'}    = "sequence_feature";
        $reportRef->{'biomarker_type'}    = "FEATURE";
        $reportRef->{'aberration_type'}    = "Unknown";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(". $reportRef->{'aberration_name'}. ")";
        unless (exists ($reportRef->{'gene'})) {
            $reportRef->{'gene'}="none";
        }
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_name'} );
    }    elsif (
            $biomarkerRef->{'snpeff'}->{'functionalEffect'} =~/intron/
           ) {
        $reportRef->{'aberration_effect'} = "intron";
        $reportRef->{'aberration_name'}   = "Intron Variant";
        $reportRef->{'biomarker_type'}    = "intron_variant";
        $reportRef->{'biomarker_type'}    = "INTRON";
        $reportRef->{'aberration_type'}    = "intron";
        $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(". $reportRef->{'aberration_name'}. ")";
        unless (exists ($reportRef->{'gene'})) {
            $reportRef->{'gene'}="none";
        }
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_name'} );
    }
      elsif ( !exists($reportRef->{'lookup'})) {
        $reportRef->{'aberration_effect'} = "Misc";
        $reportRef->{'aberration_name'}   = "Other Variant";
        $reportRef->{'biomarker_type'}    = "Unknown";
        $reportRef->{'aberration_type'}    = "Unknown";
        if (defined(         $biomarkerRef->{'snpeff'}->{'functionalEffect'})) {
            $biomarkerRef->{'functionalEffect'}    = $biomarkerRef->{'snpeff'}->{'functionalEffect'};
        }
        if (exists($biomarkerRef->{'snpeff'}->{'gene'})) {
            $reportRef->{'gene'}              = $biomarkerRef->{'snpeff'}->{'gene'};
        } elsif (!exists($reportRef->{'gene'})) {
             $reportRef->{'gene'} = 'none';
        }
        $reportRef->{'effect'}            = $reportRef->{'aberration_name'} ;
        $reportRef->{'lookup'}            = $reportRef->{'aberration_name'} ;
    }
    if (exists($biomarkerRef->{'snpeff'}->{'hgvs'}) && exists($biomarkerRef->{'snpeff'}->{'hgvs_c'}) ) {
        $reportRef->{'hgvs'}=$biomarkerRef->{'snpeff'}->{'hgvs_c'} . "(" . $biomarkerRef->{'snpeff'}->{'hgvs'} . ")";
    } elsif (exists($biomarkerRef->{'snpeff'}->{'hgvs_c'})) {
        $reportRef->{'hgvs'}=$biomarkerRef->{'snpeff'}->{'hgvs_c'} ;
    }
    return $reportRef;
}

sub assignAberrationDoc {
    my $function      = shift;
    my $biomarkerRef  = shift;
    my $RunParameters = shift;
    unless ( exists( $biomarkerRef->{'effect'} ) ) { $biomarkerRef->{'effect'} = "Unknown"; }
    unless ( exists( $biomarkerRef->{'effect'} ) ) { $biomarkerRef->{'effect'} = "Unknown"; }
    my $reportRef = { 'aberration_filter' => "FAIL", 'effect' => $biomarkerRef->{'effect'}, 'aberration_name' => "Unknown" };
    if ( $biomarkerRef->{'type'} eq "smallDeletion" || $biomarkerRef->{'type'} eq "smallInsertion" || $biomarkerRef->{'type'} eq "SNV" || $biomarkerRef->{'type'} eq 'blockSubstitution' ) {
        if ( exists( $biomarkerRef->{'snpeff'}->{'functionalEffect'} ) ) {
            $reportRef = Annotate->assignAberrationDocForSnvs( $reportRef, $biomarkerRef, $RunParameters );
        }
    }
    elsif ( $biomarkerRef->{'type'} eq "StructuralVariant" ) {
        $reportRef->{'aberration_type'}   = "CNV_LOSS";
        $reportRef->{'aberration_name'}   = "Structural Variant Breakpoint";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $biomarkerRef->{'significance'}    = "Yes";
        $biomarkerRef->{'functionalEffect'}    = "StructuralVariantBreakpoint";
        $reportRef->{'aberration_value'}  = "Structural Variant Breakpoint";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Breakpoint)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq 'FocalDeletion' ) {
        $reportRef->{'aberration_type'}   = "CNV_LOSS";
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_name'}   = "Focal Copy Number Loss";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $biomarkerRef->{'functionalEffect'}    = "FocalDeletion";
        if ( exists( $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} ) ) {
            $reportRef->{'aberration_value'} = sprintf( "%0.1f", $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} );
        }
        else { $reportRef->{'aberration_value'} = "-"; }
        $reportRef->{'biomarker_type'} = "gene";
        $reportRef->{'gene'}           = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}         = $reportRef->{'gene'} . " (Deleted $reportRef->{'aberration_value'})";
        $reportRef->{'lookup'}         = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq 'FocalGain' ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'} = "CNV_GAIN";
        $reportRef->{'aberration_name'} = "Focal Copy Number Gain";
        if ( exists( $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} ) ) {
            $reportRef->{'aberration_value'} = sprintf( "%0.1f", $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} );
        }
        else { $reportRef->{'aberration_value'} = "-"; }
        $biomarkerRef->{'functionalEffect'}    = "FocalCopyNumberGain";
        $reportRef->{'aberration_effect'} = "GeneGained";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Gained $reportRef->{'aberration_value'})";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq 'LargeInsertion' ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'}   = "LARGE_INSERTION";
        $reportRef->{'aberration_name'}   = "Large Insertion";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $biomarkerRef->{'functionalEffect'}    = "LargeInsertion";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(Insert)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq 'GeneEvent' ) {
        $reportRef->{'aberration_type'}   = "GENE_EVENT";
        $reportRef->{'aberration_name'}   = "Gene Event";
        $reportRef->{'aberration_effect'} = "GeneEvent";
        $biomarkerRef->{'functionalEffect'}    = "GeneEvent";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(event)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq 'LargeDeletion' ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'}   = "CNV_LOSS";
        $reportRef->{'aberration_name'}   = "Large Deletion";
        $biomarkerRef->{'functionalEffect'}    = "LargeDeletion";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Large Deletion)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq "OverExpressed" ) {
        $biomarkerRef->{'significance'}    = "Yes";
        if ( exists( $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} ) ) {
            $reportRef->{'aberration_value'} = sprintf( "%0.1f", $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} );
        }
        else { $reportRef->{'aberration_value'} = "-"; }
        $reportRef->{'aberration_name'}   = "Over Expressed Gene";
        $reportRef->{'aberration_type'}   = "EXP_OVER";
        $biomarkerRef->{'functionalEffect'}    = "OverExpressedGene";
        $reportRef->{'aberration_effect'} = "GeneGained";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Overexpressed $reportRef->{'aberration_value'})";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq "UnderExpressed" ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'} = "EXP_UNDER";
        if ( exists( $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} ) ) {
            $reportRef->{'aberration_value'} = sprintf( "%0.1f", $biomarkerRef->{'variants'}->[0]->{'LOG2FC'} );
        }
        else { $reportRef->{'aberration_value'} = "-"; }
        $reportRef->{'aberration_name'}   = "Under Expressed Gene";
        $reportRef->{'aberration_effect'} = "GeneLost";
        $biomarkerRef->{'functionalEffect'}    = "UnderExpressedGene";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . " (Underexpressed $reportRef->{'aberration_value'})";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    elsif ( $biomarkerRef->{'type'} eq "FusedGenes" ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'}   = "FUSED_GENE";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'variants'}->[0]->{'qual'};
        $reportRef->{'aberration_name'}   = "Fused Genes";
        $reportRef->{'aberration_effect'} = "FusedGenes";
        $biomarkerRef->{'functionalEffect'}    = "FusedGenes";
        $reportRef->{'biomarker_type'}    = "fusion";
        $reportRef->{'gene'}              = $biomarkerRef->{'biomarker'};
        $reportRef->{'effect'}            = $biomarkerRef->{'biomarker'} . "(Fusion)";
        $reportRef->{'lookup'}            = $biomarkerRef->{'biomarker'};
    }
    elsif ( $biomarkerRef->{'type'} eq "TranscriptVariant" ) {
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_type'}   = "TRANSCRIPT";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'variants'}->[0]->{'qual'};
        $reportRef->{'aberration_name'}   = "Transcript Variant";
        $reportRef->{'aberration_effect'} = "TranscriptVariant";
        $reportRef->{'biomarker_type'}    = "transcript";
        $biomarkerRef->{'functionalEffect'}    = "TranscriptVariant";
        $reportRef->{'gene'}              = $biomarkerRef->{'biomarker'};
        $reportRef->{'effect'}            = $biomarkerRef->{'biomarker'} . "(Altered Transcript)";
        $reportRef->{'lookup'}            = $biomarkerRef->{'biomarker'};
    } elsif (
            $biomarkerRef->{'type'} eq 'LargeInsertion'
            || $biomarkerRef->{'type'} eq 'LargeDeletion'
            || $biomarkerRef->{'type'} eq 'LargeBlock'
            || $biomarkerRef->{'type'} eq 'NovelInsertion'
            || $biomarkerRef->{'type'} eq 'Inversion'
            || $biomarkerRef->{'type'} eq 'CnvChange'
            || $biomarkerRef->{'type'} eq 'TandemDuplication'
            || $biomarkerRef->{'type'} eq 'MobileElementionDeletion'
            || $biomarkerRef->{'type'} eq 'MobileElementionInsertion'
            || $biomarkerRef->{'type'} eq 'LargeBlockSubstitution'
            || $biomarkerRef->{'type'} eq 'variant'){
        $reportRef->{'aberration_type'}   = "OTHERVARIANT";
        $biomarkerRef->{'significance'}    = "Yes";
        $reportRef->{'aberration_value'}  = $biomarkerRef->{'variants'}->[0]->{'qual'};
        $reportRef->{'aberration_name'}   = "Other Variant";
        $reportRef->{'aberration_effect'} = "OtherVariant";
        $reportRef->{'biomarker_type'}    = "gene";
        $reportRef->{'gene'}              = $biomarkerRef->{'gene'};
        $reportRef->{'effect'}            = $reportRef->{'gene'} . "(Other)";
        $reportRef->{'lookup'}            = join( "_", $reportRef->{'gene'}, $reportRef->{'aberration_type'} );
    }
    ####  CHECK TO SEE IF REPORTED #####
    if (exists($reportRef->{'gene'})) {
        if ($reportRef->{'gene'} eq "None" || $reportRef->{'gene'} eq "none") {
            $biomarkerRef->{'significance'}    = "None";
        }
    }
    if ( exists( $RunParameters->{'reported_aberration_names'}->{ $reportRef->{'aberration_name'} } ) ) {
        if ( $RunParameters->{'reported_aberration_names'}->{ $reportRef->{'aberration_name'} } eq "1" ) {
            $reportRef->{'aberration_filter'} = "PASS";
        }
        elsif ( $RunParameters->{'reported_aberration_names'}->{ $reportRef->{'aberration_name'} } eq "0" ) {
            $reportRef->{'aberration_filter'} = "FAIL";
        }
    }
    if ( exists( $RunParameters->{'reported_aberration_names'}->{'ALL'} ) ) {
        if ( $RunParameters->{'reported_aberration_names'}->{'ALL'} == 1 ) {
            $reportRef->{'aberration_filter'} = "PASS";
        }
    }
    
    if ( $biomarkerRef->{'variants'}->[0]->{'filter'} eq "LOWQC" && $reportRef->{'aberration_filter'} eq "PASS" ) {
        $reportRef->{'aberration_filter'} = "LOWQC";
    }
    if ( exists( $reportRef->{'aberration_name'} ) ) { $reportRef->{'aberration_type2'} = $reportRef->{'aberration_name'} }
    $biomarkerRef->{'aberration'} = $reportRef;
    $biomarkerRef->{'effect'}     = $biomarkerRef->{'aberration'}->{'effect'};
    return $biomarkerRef;
}

sub isCommon {
    my $function     = shift;
    my $biomarkerRef = shift;
    my $maxMaf       = -1;
    my $freq         = $maxMaf;
    my $common       = 'private';
    my $database="Unknown";
    foreach my $filter ('1000Gp3_AFR_AF','1000Gp3_AF','1000Gp3_AFR_AF','1000Gp3_EUR_AF','1000Gp3_AMR_AF','1000Gp3_EAS_AF','1000Gp3_SAS_AF','TWINSUK_AF', 'ALSPAC_AF','ESP6500_AA_AF','ESP6500_EA_AF','ExAC_AF','ExAC_Adj_AF','ExAC_AFR_AF','ExAC_AMR_AF','ExAC_EAS_AF','ExAC_FIN_AF','ExAC_NFE_AF', 'ExAC_SAS_AF','TWINSUK_AF','SAS_maf','FIN_maf','NFE_maf','AMR_maf','EAS_maf','AFR_maf') {
        if ( exists( $biomarkerRef->{'annotate'}->{$filter} ) ) {
            if ( $biomarkerRef->{'annotate'}->{$filter} > $maxMaf ) { $maxMaf = $biomarkerRef->{'annotate'}->{$filter} ; $database=$filter}
        }
    }
    if ( $maxMaf >= 0 && $maxMaf < 0.025 ) {
        $common = 'rare';
        $freq   = $maxMaf;
    }
    elsif ( $maxMaf > 0.025 || $biomarkerRef->{'count'}->{'all'} > 50 ) {
        $common = 'common';
        $freq   = $maxMaf;
    }
    elsif ( exists( $biomarkerRef->{'annotate'}->{'dbsnp'} ) ) {
        $common = 'documented';
    }
    else {
        $common = 'private';
    }
    $biomarkerRef->{'annotate'}->{'frequency'}         = $freq;
    $biomarkerRef->{'annotate'}->{'maxPopFreq'}         = $freq;
    $biomarkerRef->{'annotate'}->{'maxPopFreqDatabase'}         = $database;
    MongoDB::force_double($biomarkerRef->{'annotate'}->{'frequency'});
    MongoDB::force_double($biomarkerRef->{'annotate'}->{'maxPopFreq'});
    $biomarkerRef->{'annotate'}->{'frequencyCategory'} = $common;
    return $biomarkerRef;
}

sub addSnpeffInfo {
    my ( $function, $val, $snpeffVariantRef,$snpeffMultiplesVariantRef, $RunParameters, $geneListRef ) = @_;
    my $deg = 0;
    my @snpeffs = split( /\,/, $val );
    my $max=0;
    if (exists($RunParameters->{'allSnpeffTranscripts'})) {
        if ($RunParameters->{'allSnpeffTranscripts'}==1) {
            $max=$#snpeffs;
        }
    } else {
        $RunParameters->{'allSnpeffTranscripts'}=0;
    }
    LOOP1: for (my $l=0;$l<=$#snpeffs;++$l) {
        if ($snpeffs[$l] =~ /sequence_feature\|(.*?)\|.*?\|.*?\|(.*?)\|.*/) {
            my $match1=$1;
            my $match2=$2;
            if (exists($snpeffVariantRef->{'featureStrength'})) {
                $snpeffVariantRef->{'featureStrength'}.=", " . $match1;
                $snpeffVariantRef->{'featureDescription'}.=", " . $match2;
            } else {
                $snpeffVariantRef->{'featureStrength'}=$1;
                $snpeffVariantRef->{'featureDescription'}=$2;
            }
        }
    }

    LOOP: for (my $k=0;$k<=$#snpeffs;++$k) {
        my @effAr = ( $snpeffs[$k] =~ /(.*?)\|/g );
        unless (exists($RunParameters->{'snpeffValidAnnotations'}->{'ALL'})) {$RunParameters->{'snpeffValidAnnotations'}->{'ALL'}=0}
        if (defined($effAr[7]) || $RunParameters->{'allSnpeffTranscripts'}==1 || $RunParameters->{'snpeffValidAnnotations'}->{'ALL'}==1) {
            $snpeffVariantRef->{'functionalEffect'}="Unknown";
            $snpeffVariantRef->{'functionalImpact'}="NONE";
            if ($effAr[7] eq "protein_coding") {
               (
                   $snpeffVariantRef->{'allele'}, $snpeffVariantRef->{'functionalEffect'}, $snpeffVariantRef->{'functionalImpact'},$snpeffVariantRef->{'gene'},
                   $snpeffVariantRef->{'geneId'},  $snpeffVariantRef->{'featureType'},$snpeffVariantRef->{'featureId'}, $snpeffVariantRef->{'transcriptBiotype'},
                   $snpeffVariantRef->{'rank'},     $snpeffVariantRef->{'hgvs_c'},$snpeffVariantRef->{'hgvs'},$snpeffVariantRef->{'cDna/length'},$snpeffVariantRef->{'aa/length'},
                   $snpeffVariantRef->{'distance'},$snpeffVariantRef->{'errors'}
               ) = @effAr;
            }
            if ($RunParameters->{'snpeffValidAnnotations'}->{'ALL'}==1) {
               (
                   $snpeffVariantRef->{'allele'}, $snpeffVariantRef->{'functionalEffect'}, $snpeffVariantRef->{'functionalImpact'},$snpeffVariantRef->{'gene'},
                   $snpeffVariantRef->{'geneId'},  $snpeffVariantRef->{'featureType'},$snpeffVariantRef->{'featureId'}, $snpeffVariantRef->{'transcriptBiotype'},
                   $snpeffVariantRef->{'rank'},     $snpeffVariantRef->{'hgvs_c'},$snpeffVariantRef->{'hgvs'},$snpeffVariantRef->{'cDna/length'},$snpeffVariantRef->{'aa/length'},
                   $snpeffVariantRef->{'distance'},$snpeffVariantRef->{'errors'}
               ) = @effAr;
            }
            $snpeffVariantRef->{'report'}='comprehensive';
            $snpeffVariantRef->{'impact'}='low';
            if ( $snpeffVariantRef->{'functionalImpact'} eq "HIGH" || $snpeffVariantRef->{'functionalImpact'} eq "MODERATE") {
                $snpeffVariantRef->{'report'}='standard';
                $snpeffVariantRef->{'impact'}='significant';
            }
            unless (defined($snpeffVariantRef->{'gene'})) {
                  $snpeffVariantRef->{'gene'}="none";
            }
            if (!(exists($geneListRef->{ $RunParameters->{'defaultAssay'} }->{$snpeffVariantRef->{'gene'}} ) ) && 
            !(exists($geneListRef->{ $RunParameters->{'defaultAssay'} }->{'ALL'} ) ) ) {
                            next LOOP
                        }
            unless (defined($snpeffVariantRef->{'functionalEffect'})) {
                  $snpeffVariantRef->{'functionalEffect'}="Unknown";
            }
            $snpeffVariantRef->{'functionalEffects'}=$snpeffVariantRef->{'functionalEffect'};
            $snpeffVariantRef->{'functionalEffect'}=~s/\&.*//g;
            foreach my $key (keys %{$snpeffVariantRef}) {
                unless (defined($snpeffVariantRef->{$key})) {delete $snpeffVariantRef->{$key}}
            }
            if (exists($RunParameters->{'snpeffValidAnnotations'}->{'ALL'})){
                    if (exists($snpeffVariantRef->{'hgvs'})) {
                        $snpeffVariantRef->{'AminoAcidChange'} = $snpeffVariantRef->{'hgvs'};
                    }
                    $snpeffVariantRef->{'toJoinSnpeff'}     = int(1);
                    $snpeffVariantRef->{'toAnnotateByGene'} = int(1);
            } elsif(exists ($snpeffVariantRef->{'functionalEffect'})) {
                 if ( exists($RunParameters->{'snpeffValidAnnotations'}->{$snpeffVariantRef->{'functionalEffect'}})) {
                     if ($RunParameters->{'snpeffValidAnnotations'}->{$snpeffVariantRef->{'functionalEffect'}} >0 ) {
                          if (exists($snpeffVariantRef->{'hgvs'})) {
                              $snpeffVariantRef->{'AminoAcidChange'} = $snpeffVariantRef->{'hgvs'};
                          }
                          $snpeffVariantRef->{'toJoinSnpeff'}     = int(1);
                          $snpeffVariantRef->{'toAnnotateByGene'} = int(1);
                     } else {
                          $snpeffVariantRef->{'toJoinSnpeff'}     = int(0);
                          $snpeffVariantRef->{'toAnnotateByGene'} = int(0);
                     }
                 }
            }
            if (exists($snpeffVariantRef->{'distance'})) {
                if ($snpeffVariantRef->{'distance'}=~/(.*?)\/(.*)/) {
                     $snpeffVariantRef->{'codon'}=$1;
                     $snpeffVariantRef->{'geneCodon'}=$snpeffVariantRef->{'gene'} . "_" .$snpeffVariantRef->{'codon'};
                }
            }
            if (exists($snpeffVariantRef->{'rank'})) {
                 if ($snpeffVariantRef->{'rank'}=~/(.*?)\/(.*)/) {
                      $snpeffVariantRef->{'exon'}=$1;
                      $snpeffVariantRef->{'gene_exon'}=$snpeffVariantRef->{'gene'} . "_" .$snpeffVariantRef->{'exon'};
                      $snpeffVariantRef->{'geneId_exon'}=$snpeffVariantRef->{'geneId'} . "_" .$snpeffVariantRef->{'exon'};
                 }
            }
            if ($RunParameters->{'allSnpeffTranscripts'}==1 ) {
                 for ( my $i = 0 ; $i <= $max ; ++$i ) {
                     my @effAr = ( $snpeffs[$i] =~ /(.*?)\|/g );
                     my $temp={};
                     (
                         $temp->{'allele'}, $temp->{'functionalEffect'}, $temp->{'functionalImpact'},$temp->{'gene'},
                         $temp->{'geneId'},  $temp->{'featureType'},$temp->{'featureId'}, $temp->{'transcriptBiotype'},
                         $temp->{'rank'},     $temp->{'hgvs_c'},$temp->{'hgvs'},$temp->{'cDna/length'},$temp->{'aa/length'},
                         $temp->{'distance'},$temp->{'errors'}
                     ) = @effAr;
                     foreach my $key (keys %{$temp}) {
                         unless (defined($temp->{$key})) {delete $temp->{$key}}
                     }
                     push (@{$snpeffMultiplesVariantRef},$temp);
                }
            }
            if ( exists( $snpeffVariantRef->{'AminoAcidChange'} ) && exists($snpeffVariantRef->{'toJoinSnpeff'} )) {
                if ( length( $snpeffVariantRef->{'AminoAcidChange'} ) > 0 ) {
                    if ( $snpeffVariantRef->{'AminoAcidChange'} =~ /p\.(.*)/ ) {
                        $snpeffVariantRef->{'AminoAcidChange'} = aaConv($1);
                    }
                } else {
                    $snpeffVariantRef->{'AminoAcidChange'}="None";
                }
            } else {
                $snpeffVariantRef->{'AminoAcidChange'}="None";
            }
             last LOOP;
        }
    }
    return $snpeffVariantRef, $snpeffMultiplesVariantRef
}

sub aaConv {
    my $aa = shift;
    $aa =~ s/Ala/A/g;
    $aa =~ s/Arg/R/g;
    $aa =~ s/Asn/N/g;
    $aa =~ s/Asp/D/g;
    $aa =~ s/Asx/B/g;
    $aa =~ s/Cys/C/g;
    $aa =~ s/Glu/E/g;
    $aa =~ s/Gln/Q/g;
    $aa =~ s/Glx/Z/g;
    $aa =~ s/Gly/G/g;
    $aa =~ s/His/H/g;
    $aa =~ s/Ile/I/g;
    $aa =~ s/Lys/K/g;
    $aa =~ s/Leu/L/g;
    $aa =~ s/Met/M/g;
    $aa =~ s/Phe/F/g;
    $aa =~ s/Pro/P/g;
    $aa =~ s/Ser/S/g;
    $aa =~ s/Thr/T/g;
    $aa =~ s/Trp/W/g;
    $aa =~ s/Tyr/Y/g;
    $aa =~ s/Val/V/g;
    return $aa;
}

sub joinOn {
#    my ( $function, $key, $val, $variantType ) = @_;
    my $function=shift;
    my $key=shift;
    my $val=shift;
    my $variantType=shift;
    my $findFieldsRef;
    my $hashref = {};
    my $hashref2 = {};

    if ($findFieldsRef=shift) {
        foreach my $FieldsKey (keys %{$findFieldsRef}) {
            if ($findFieldsRef->{$FieldsKey}==1) {
                 $hashref2->{$FieldsKey}=1;
            }
        }
        my $cur = $variantType->find( { $key => $val },{limit=>1} )->fields($hashref2);
        if ( my $col = $cur->next ) {
            my %hashcopy = %$col;
            delete $hashcopy{'_id'};
            $hashref = \%hashcopy;
            $hashref->{'join'} = 1;
            foreach my $hashEle ( keys %{$hashref} ) {
                    if ( $hashref->{$hashEle} eq "." || $hashref->{$hashEle} eq "" || $hashref->{$hashEle} eq ".\r" ) {
                        undef $hashref->{$hashEle};
                        delete $hashref->{$hashEle};
                    }
            }
        }else {
            $hashref->{'join'} = 0;
        }
    } else {
        my $cur = $variantType->find( { $key => $val } ,{limit=>1});
        if ( my $col = $cur->next ) {
            my %hashcopy = %$col;
            delete $hashcopy{'_id'};
            $hashref = \%hashcopy;
            $hashref->{'join'} = 1;
            foreach my $hashEle ( keys %{$hashref} ) {
                if (!defined($hashref->{$hashEle})) { print "\t\tWarn $hashEle not defined in $hashref->{'_id'} for $key=$val\n";}
                if ( $hashref->{$hashEle} eq "." || $hashref->{$hashEle} eq "" || $hashref->{$hashEle} eq ".\r" ) {
                    if (exists($hashref->{$hashEle})) {
                      undef $hashref->{$hashEle};
                      delete $hashref->{$hashEle};
                    }
                }
            }
        }else {
            $hashref->{'join'} = 0;
        }
    }
    return $hashref;
}    # joinOn ($key,$val,$variantType)=$hashRef

sub loadLibrary {    #Annotate->loadLibrary($filename,$geneListRef,$libraryToPrune,$RunParameters);

    my $function      = shift;
    my $filename      = shift;
    my $geneListRef   = shift;
    my $lib           = shift;
    my $RunParameters = shift;
    open( DB, "$filename" ) or die "$filename\n";
    while (<DB>) {
        my $line = $_;
        chomp($line);
        $geneListRef->{$lib}->{$line} = $lib;
    }
    close(DB);
    return $geneListRef;
}

sub filterOutGenes {
##            'filterOutGenes'         => { 'seuratPointVariantFile' => {'IsOneOf'=>'pathToFile'}, 
##                                          'seuratPointMutationFile' => {'IsNotOneOf'=>'pathToFile'}} 
##                                         },      
    my $function=shift;
    my $variantRef=shift;
    my $RunParameters       = shift; 
    my $gene=shift;
    my $isGene=$RunParameters->{'isGene'};    
    unless (defined($gene)) {return 1}
    if (length($gene)==0) {return 1}
    Maintain->reloadStartupRunParameters($RunParameters); 
    my $geneListRef={};
    foreach my $libraryToPrune (keys %{$RunParameters->{'reportableGenes'}}) {
         my $filename=$RunParameters->{'reportableGenes'}->{$libraryToPrune};
         Annotate->loadLibrary($filename,$geneListRef,$libraryToPrune,$RunParameters);
    }
#    print "\t+Filtering out genes\n";       
    CYCLE: foreach my $fileFilterType (keys %{$RunParameters->{'filterOutGenes'}}) {
        unless ($variantRef->{'filetype'} eq $fileFilterType) {
            next CYCLE;
        }
        my $keyCount=0;
        foreach my $operator (keys %{$RunParameters->{'filterOutGenes'}->{$fileFilterType}}) {
            if ($operator ne "IsOneOf" && $operator ne "IsNotOneOf") {
                die "\t\t-DIE $operator is not valid operator!\n"; 
            }
            if ($keyCount>0) {
                die "\t\t DIE  We haven't added support for multiple operators, best to not do this even if it could work. Lets just die because Ahmet is right\n";
            }
            ++$keyCount;
            ## Load genelist if it isn't already
            unless (exists($RunParameters->{'filterOutGenes_' . $variantRef->{filetype}})) {
                my $filename=$RunParameters->{'filterOutGenes'}->{$fileFilterType}->{$operator};
                print "\t\t\t-Loading file to filterOutGenes: $filename\n";
                open (FILENAME, "$filename") or die "\t\tDIE! Can't open filename: $filename\n";
                while (<FILENAME>) {
                   chomp;
                   my $geneToFilterOut=$_;
                   if (exists($geneListRef->{$variantRef->{'assay'}}->{$geneToFilterOut})) {
                       $RunParameters->{'filterOutGenes_' . $variantRef->{filetype}}->{$geneToFilterOut}=1;
                       print "\t\t+Gene being filtered from filetype $variantRef->{filetype} if $operator($geneToFilterOut)\n";
                   } else {
                      die "\t\tDIE: $filename should only contain valid genes for the appropriate library '$geneToFilterOut' was not found!\n"
                   }
                }
                close (FILENAME);
            }
            #  Checking operators
            if ($operator eq "IsOneOf") {
                my $found=0;
                foreach my $geneToFilterOut (keys %{$RunParameters->{'filterOutGenes_' . $variantRef->{filetype}}}) {
                   if ($gene eq $geneToFilterOut) {
                       ++$found;
                   }
                }                          
                if ($found>0) {
                    return 0
                }                            
            } elsif ($operator eq "IsNotOneOf") {
                my $found=0;
                foreach my $geneToFilterOut (keys %{$RunParameters->{'filterOutGenes_' . $variantRef->{filetype}}}) {
                    if ($gene eq $geneToFilterOut) {
                       ++$found;
                    }                        
                }
                if ($found<1) {
                    return 0
                }
            } else {
               die "\t\tDIE  I don't think this $operator is valid\n"; 
            }            
        }
    }
    return 1;
}

sub forceBiomarkerType {
    my $function      = shift;
    my $forceTypeOf   = shift;
    my $biomarkerRef  = shift;
    my $RunParameters = shift;
    if ( exists( $RunParameters->{ $forceTypeOf . 'ForceType' } ) ) {
        foreach my $key ( keys %{ $RunParameters->{ $forceTypeOf . 'ForceType' } } ) {
            if ( exists( $RunParameters->{ $forceTypeOf . 'ForceType' }->{$key} ) ) {
                if ( $RunParameters->{ $forceTypeOf . 'ForceType' }->{$key} == 2 ) {
                    if ( exists( $biomarkerRef->{$forceTypeOf}->{$key} ) ) {
                        unless ( looks_like_number( $biomarkerRef->{$forceTypeOf}->{$key} ) ) {
                            $biomarkerRef->{$forceTypeOf}->{$key} = -1;
                        }
                        MongoDB::force_double( 1.00 * $biomarkerRef->{$forceTypeOf}->{$key} );
                    }
                }
                elsif ( $RunParameters->{ $forceTypeOf . 'ForceType' }->{$key} == 3 ) {
                    if ( exists( $biomarkerRef->{$forceTypeOf}->{$key} ) ) {
                        unless ( looks_like_number( $biomarkerRef->{$forceTypeOf}->{$key} ) ) {
                            $biomarkerRef->{$forceTypeOf}->{$key} = -1;
                        }
                        $biomarkerRef->{$forceTypeOf}->{$key} = int( sprintf( "%.0f", $biomarkerRef->{$forceTypeOf}->{$key} ) );
                    }
                }
                elsif ( $RunParameters->{ $forceTypeOf . 'ForceType' }->{$key} == 1 ) {
                    $biomarkerRef->{$forceTypeOf}->{$key} = "$biomarkerRef->{$forceTypeOf}->{$key}";
                }
                elsif ( $RunParameters->{ $forceTypeOf . 'ForceType' }->{$key} == 0 ) {
                    undef( $biomarkerRef->{$forceTypeOf}->{$key} );
                    delete( $biomarkerRef->{$forceTypeOf}->{$key} );
                }
            }
        }
    }
}
1;
