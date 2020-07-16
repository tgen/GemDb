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
# * Major Contributor(s):
#    David Craig
#/
######################################################################################

$|=1;
package insertRules;
use POSIX;
sub loadTo {
      my $filesInsert=addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $dbName=shift;
      my $filename=$filePath;
      $filename=~s/.*\///g;
      print "-Module insertRules function:$function filepath:$filePath collection:$collection hostname:$hostName  adding rules $filename\n";      
      my ($CONN,$DB,$RULES)=initMongo->new('AnnotateBy',$collection,$hostName);
      $RULES -> remove({filename=>$filename});
      ##CCDS###
      open (CCDS,"CCDS.txt") or die "Can't find CCDS.txt\n";
      while (<CCDS>) {
         chomp;
         $ccds{$_}=1;
      }
      close (CCDS);
      open (RULEFILE,"$filePath") or die "Can't open file\n";
      $headerLine=(<RULEFILE>);
      chomp($headerLine);
      $headerLine=~s/\./\+/g;
	@headArray=split(/\t/,$headerLine); 
      $RULES->ensure_index({'filename'=>1});      
      $RULES->remove({'RulesVersion'=>$filename});      
      LOOP: while (<RULEFILE>) {
          $line=$_;
	  chomp($line);
	  @fields=split(/\t/,$line);
          my %ref=();
          my $toInsert=\%ref;
          foreach ($i=0;$i<=$#headArray;++$i) {
		if (!defined($fields[$i])){ $fields[$i]="NULL"}      
               if (length($fields[$i])==0) { $fields[$i]="NULL"}
               $toInsert->{$headArray[$i]}=$fields[$i];
          }          
          $toInsert->{'filepath'}=$filePath;
          $toInsert->{'filename'}=$filename;
          $toInsert->{'RulesVersion'}=$filename;          
          $toInsert->{'drug_rules_version'}=$filename;
	if ($toInsert->{'biomarker_type'}=~/fused_gene/i) {
                     $toInsert->{'class'}="biomarker";
                     $toInsert->{'type'}=uc($toInsert->{"biomarker_type"});
                     $toInsert->{'entrez'}=$toInsert->{"biomarker_entrez"};
                     $toInsert->{"biomarker_aberration_type"}=uc($toInsert->{"biomarker_type"});
                     $toInsert->{'aberration_type'}=uc($toInsert->{"biomarker_aberration_type"});
                     if ($toInsert->{'biomarker_symbol'}=~/\[CCDS\]/) {
                           $toInsert->{"original_symbol"} = $toInsert->{"biomarker_symbol"};
                           foreach my $gene (keys %ccds) {
                              my $new =$toInsert->{"original_symbol"};
                              $new=~s/_\[CCDS\]/_$gene/;
                              $new=~s/\[CCDS\]_/$gene\_/;                              
                              $toInsert->{'biomarker_symbol'}=$new;    
                              $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};                                                                   
                              $toInsert->{'drug_rule_fusion_wild_card_match'}="CCDS";                              
                              $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'}; 
                              $toInsert->{'Drug Name'}=$toInsert->{'drug'};              
                              $RULES->insert($toInsert,{'safe'=>1});        
                           }    
                           next LOOP;                                                                    
                     } else {
                           $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};                 
                         $toInsert->{'aberration_value'}=$toInsert->{'biomarker_symbol'} ;
                         $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'};  
                         $toInsert->{'Drug Name'}=$toInsert->{'drug'};              
                         $RULES->insert($toInsert,{'safe'=>1});        
                         next LOOP;                       
                     }
          } elsif ($toInsert->{'biomarker_type'} eq "gene" || $toInsert->{'biomarker_aberration_type'}=~/snv/i) {
                 if ($toInsert->{'biomarker_transcript_id'} ne "EMPTY") {
                     $toInsert->{'class'}="biomarker";
                     $toInsert->{'type'}="TRANSCRIPT";
                     $toInsert->{'entrez'}=$toInsert->{"biomarker_entrez"};
                     $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};
                     $toInsert->{'aberration_type'}="TRANSCRIPT";
                     $toInsert->{'aberration_value'}=$toInsert->{'biomarker_transcript_id'} ;   
                     $toInsert->{'aberration_lookup'}=$toInsert->{'biomarker_transcript_id'};     
                     $toInsert->{'Drug Name'}=$toInsert->{'drug'};           
                     $RULES->insert($toInsert,{'safe'=>1});        next LOOP;             
                 } elsif ($toInsert->{'biomarker_aberration_type'}=~/snv/i) {
                     $toInsert->{'class'}="modifier";
                     $toInsert->{'type'}=uc($toInsert->{"biomarker_type"});
                     $toInsert->{'entrez'}=$toInsert->{"biomarker_entrez"};
                     $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};
                     $toInsert->{'aberration_type'}=uc($toInsert->{"biomarker_aberration_type"});                     
                     $toInsert->{'aberration_value'}=$toInsert->{"biomarker_aberration_value"};
                     if ($toInsert->{'aberration_value'}=~/#/) {
                            my $origVal=$toInsert->{'aberration_value'};
                            $toInsert->{'originalValue'}=$origVal;
                            foreach my $aa ('A','R','N','D','B','C','E','Q','Z','G','H','I','L','M','F','P','S','T','W','Y','V') {
                               $toInsert->{'aberration_value'}=~s/#/$aa/; 
                               $toInsert->{'drug_rule_snv_wild_card_match'}="#";
                               $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . uc($toInsert->{'aberration_value'});   
                               $RULES->insert($toInsert,{'safe'=>1});                     
                               $toInsert->{'aberration_value'}=$origVal;
                            }
                            next LOOP;
                     } elsif ($toInsert->{'aberration_value'} eq "change") {
                               $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                               $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'} . "_" . uc($toInsert->{'aberration_value'});       
                               $toInsert->{'drug_rule_snv_wild_card_match'}="#";                              
                               $RULES->insert($toInsert,{'safe'=>1});    next LOOP;                 
                     } else {
                               $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                               $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . uc($toInsert->{'aberration_value'});   
                               $RULES->insert($toInsert,{'safe'=>1});    next LOOP;             
                     }    
                } elsif ($toInsert->{'biomarker_aberration_type'}=~/cnv/i || $toInsert->{'biomarker_aberration_type'}=~/exp/i || $toInsert->{'biomarker_aberration_type'}=~/prot_exp/i ) {           
                   $toInsert->{'class'}="biomarker";
                   $toInsert->{'type'}=uc($toInsert->{"biomarker_type"});
                   $toInsert->{'entrez'}=$toInsert->{"biomarker_entrez"};
                   $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};
                   $toInsert->{'aberration_type'}=uc($toInsert->{"biomarker_aberration_type"});
                   $toInsert->{'aberration_value'}=uc($toInsert->{"biomarker_aberration_value"});
                   $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                   $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'} . "_" . uc($toInsert->{'aberration_value'});       
                   $RULES->insert($toInsert,{'safe'=>1});    next LOOP;                         
                } elsif ($toInsert->{'biomarker_aberration_type'}=~/small_insertion/ || $toInsert->{'biomarker_aberration_type'} eq "small_deletion" ) {     
                   $toInsert->{'class'}="biomarker";
                   $toInsert->{'type'}=uc($toInsert->{"biomarker_type"});
                   $toInsert->{'entrez'}=$toInsert->{"biomarker_entrez"};
                   $toInsert->{'symbol'}=$toInsert->{"biomarker_symbol"};
                   $toInsert->{'aberration_type'}=uc($toInsert->{"biomarker_aberration_type"});
                   $toInsert->{'aberration_value'}=uc($toInsert->{"biomarker_aberration_value"});
                   $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'};  
                   $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                   $RULES->insert($toInsert,{'safe'=>1});      next LOOP;           
                }    
           } elsif ($toInsert->{'modifier_type'}=~/fused_gene/i) {
                     $toInsert->{'class'}="modifier";
                     $toInsert->{'type'}=uc($toInsert->{"modifier_type"});
                     $toInsert->{'entrez'}=$toInsert->{"modifier_entrez"};
                     $toInsert->{'symbol'}=$toInsert->{"modifier_symbol"};
                     $toInsert->{'aberration_type'}=uc($toInsert->{"modifier_aberration_type"});
                     $toInsert->{'aberration_value'}=$toInsert->{'modifier_symbol'} ;
                     $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'};              
                     $toInsert->{'Drug Name'}=$toInsert->{'drug'};  
                     $RULES->insert($toInsert,{'safe'=>1});        next LOOP;              
           } elsif ($toInsert->{'modifier_type'} eq "gene" || $toInsert->{'modifier_aberration_type'}=~/snv/) {
                 if ($toInsert->{'modifier_transcript_id'} ne "EMPTY") {
                     $toInsert->{'class'}="modifier";
                     $toInsert->{'type'}="TRANSCRIPT";
                     $toInsert->{'entrez'}=$toInsert->{"modifier_entrez"};
                     $toInsert->{'symbol'}=$toInsert->{"modifier_symbol"};
                     $toInsert->{'aberration_type'}="TRANSCRIPT";
                     $toInsert->{'aberration_value'}=$toInsert->{'modifier_transcript_id'} ;   
                     $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'};   
                     $toInsert->{'Drug Name'}=$toInsert->{'drug'};             
                     $RULES->insert($toInsert,{'safe'=>1});    next LOOP;                                                                                    
                 } elsif ($toInsert->{'modifier_aberration_type'}=~/snv/) {
                         $toInsert->{'class'}="modifier";
                         $toInsert->{'type'}=uc($toInsert->{"modifier_type"});
                         $toInsert->{'entrez'}=$toInsert->{"modifier_entrez"};
                         $toInsert->{'symbol'}=$toInsert->{"modifier_symbol"};
                         $toInsert->{'aberration_type'}=uc($toInsert->{"modifier_aberration_type"});                                          
                         $toInsert->{'aberration_value'}=$toInsert->{"modifier_aberration_value"};                     
                         if ($toInsert->{'aberration_value'}=~/#/) {
                                my $origVal=$toInsert->{'aberration_value'};
                                $toInsert->{'originalValue'}=$origVal;
                                foreach my $aa ('A','R','N','D','B','C','E','Q','Z','G','H','I','L','M','F','P','S','T','W','Y','V') {
                                   $toInsert->{'aberration_value'}=~s/#/$aa/; 
                                   $toInsert->{'drug_rule_snv_wild_card_match'}="#";
                                   $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . uc($toInsert->{'aberration_value'});   
                                   $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                                   $RULES->insert($toInsert,{'safe'=>1});                     
                                   $toInsert->{'aberration_value'}=$origVal;
                                }
                                next LOOP;
                         } elsif ($toInsert->{'aberration_value'} eq "change") {
                                  $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'} . "_" . uc($toInsert->{'aberration_value'});       
                                   $toInsert->{'drug_rule_snv_wild_card_match'}="CHANGE";    
                                   $toInsert->{'Drug Name'}=$toInsert->{'drug'};                          
                                   $RULES->insert($toInsert,{'safe'=>1});next LOOP;                         
                         } else {
                                   $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . uc($toInsert->{'aberration_value'});   
                                                          $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                                   $RULES->insert($toInsert,{'safe'=>1});
                                   next LOOP;                 
                         } 
                } elsif ($toInsert->{'modifier_aberration_type'}=~/cnv/ || $toInsert->{'modifier_aberration_type'}=~/exp/ ) {           
                       $toInsert->{'class'}="modifier";
                       $toInsert->{'type'}=uc($toInsert->{"modifier_aberration_type"});
                       $toInsert->{'entrez'}=$toInsert->{"modifier_entrez"};
                       $toInsert->{'symbol'}=$toInsert->{"modifier_symbol"};
                       $toInsert->{'aberration_type'}=uc($toInsert->{"modifier_aberration_type"});
                       $toInsert->{'aberration_value'}=uc($toInsert->{"modifier_aberration_value"});
                       $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'} . "_" . uc($toInsert->{'aberration_value'});       
                       $toInsert->{'Drug Name'}=$toInsert->{'drug'};                       
                       $RULES->insert($toInsert,{'safe'=>1});    next LOOP;                 
                } elsif ($toInsert->{'modifier_aberration_type'} eq "small_insertion" || $toInsert->{'modifier_aberration_type'} eq "small_deletion" ) {     
                       $toInsert->{'class'}="modifier";
                       $toInsert->{'type'}=uc($toInsert->{"modifier_aberration_type"});
                       $toInsert->{'entrez'}=$toInsert->{"modifier_entrez"};
                       $toInsert->{'symbol'}=$toInsert->{"modifier_symbol"};
                       $toInsert->{'aberration_type'}=uc($toInsert->{"modifier_aberration_type"});
                       $toInsert->{'aberration_value'}=uc($toInsert->{"modifier_aberration_value"});
                       $toInsert->{'aberration_lookup'}=$toInsert->{'symbol'} . "_" . $toInsert->{'aberration_type'};  
                       $toInsert->{'Drug Name'}=$toInsert->{'drug'};
                       $RULES->insert($toInsert,{'safe'=>1});    next LOOP;                                 
                }
            }
      }
      $RULES->ensure_index({'aberration_lookup'=>1});
      $RULES->ensure_index({'drug'=>1});      
      return 1;
}
sub joinTo {
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;
      my $dbName=shift;
      my $joinTo=shift;
      my $filename=$filePath;
      $filename=~s/.*\///g;
      print "-Module insertRules function:$function filepath:$filePath collection:$collection hostname:$hostName  adding rules $filename\n";      
      my ($CONN,$DB,$RULES)=initMongo->new('AnnotateBy',$collection,$hostName);
      open (RULEFILE,"$filePath") or die "Can't open file\n";
      $headerLine=(<RULEFILE>);
      chomp($headerLine);
      $headerLine=~s/\./\+/g;
      @headArray=split(/\t/,$headerLine); 
      $RULES->ensure_index({'filename'=>1});      
      LOOP: while (<RULEFILE>) { 
          $line=$_;
          chomp($line);
          @fields=split(/\t/,$line);
           print "$fields[0]\tcount:$#fields @fields\n\n\n\n";
          $toInsert={};
          foreach ($i=0;$i<=$#headArray;++$i) {
               unless (length($fields[$i])>=0) { print "going to die\n";die "empty unhappy field\n"; }
               $toInsert->{$headArray[$i]}=$fields[$i];
          }
          $test=$RULES->update({$joinTo=>$toInsert->{$joinTo}},{'$set'=>$toInsert},{'safe'=>1,multiple=>1});            
      }           
}
sub addToCollections {
    my $function=shift;
    my $filePath=shift;
    my  $collection=shift;
    my  $hostName=shift;
    my $type='rules';
    print "filePath:$filePath\n";
    my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
    my $md5sum=`md5sum $filePath`;
    chomp($md5sum);
    $md5sum=~s/ .*//g;
    my $wc=`wc -l $filePath`;
    $wc=~s/ .*//g; 
    print "mdsub: $md5sum wc: $wc\n";
    my $filesInsert={'date'=>$date,'type'=>$type,'filename'=>$filePath,'md5sum'=>$md5sum,'collection'=>$collection,'lines'=>$wc,'md5sum'=>$md5sum};
    return $filesInsert;
}
1;
