#!/usr/bin/perl -w
#dbnsfp2.4
package insertDbnsfp;
use initMongo;
use POSIX;
use Time::HiRes qw( usleep );
sub loadTo {
      print "-Module: @_\n";
#      my $md5=addToCollections(@_);
      my $function=shift;
      my $filePath=shift;
      my $collection=shift;
      my $hostName=shift;	
      my ($conn,$db,$dbnsfp)=initMongo->new('AnnotateBy',$collection,$hostName);
      my @FILES=`ls -1 db/dbnsfp3.4/dbNSFP3*_variant.chr*`;
	  foreach $file (@FILES) {
		  my $done={};
		  chomp($file);
		  print "-dbnsfp opening $file\n";
		  open(FILE,"<$file") or die "Can't open $file\n";
		  $header=(<FILE>);
		  chomp($header);
		  &addHead;
		  LOOP: while(<FILE>) {
			  next if $_ =~ m/^#/g;
			  $line=$_;
			  chomp($line);
			  @temp = split("\t",$line);                   
			  %hash=();
			  %hashAnnotate=();
			  $hash{'coord'}=join(':',"chr" . $temp[7],$temp[8],$temp[2],$temp[3]);                  
			  if (exists($done->{$hash{'coord'}})){next LOOP}               
			  foreach $i (11,14,19,20,21,22,25,31,34,41,48,51,54,56,61,64,73,176,177,178,179,180,181,182) {
				if ($temp[$i] ne ".") {
					$hash{$head[$i]}=$temp[$i];
				}
			  }
			  for $i (68,69,76,80,82,85,88,94,95,97,99,101,106) {
				  if ($temp[$i] ne "." ) {
					 $temp[$i]=~s/\;.*//g;
					 $temp[$i]=$temp[$i]*1.00;
					 MongoDB::force_double($temp[$i]);
					 $hash{$head[$i]}=$temp[$i];
				   }  		  
			  }
			  for $i (109,111,113,115,117,119,121,123,125,127,129,131,133,135,137,139,141,143,145,147,149,153,155,157,159,161,163,165,167,169,171,173,175){
				  if ($temp[$i] ne "." || $temp[$i] ne "0" ) {
					 $temp[$i]=~s/\;.*//g;
					 $temp[$i]=$temp[$i]*1.00;
					 MongoDB::force_double($temp[$i]);
					 $hash{$head[$i]}=$temp[$i];
				   }  
			  }		  
			  for $i (0,1,2,6,7,8,9,10) {
				 if ($temp[$i] ne ".") {
					 $hash{$head[$i]}=int($temp[$i]);
				 }
			  }			  
			  $done->{$hash{'coord'}}=1;
			  $dbnsfp->update({'coord'=>$hash{'coord'}},{'$set'=>\%hash},{'upsert'=>1,'safe'=>1});
     		}
		close (FILE);
	 }
  $idx = Tie::IxHash->new("coord"=> 1); 
  $dbnsfp->ensure_index($idx);
}

sub addToCollections {
      my $function=shift;
      my $filePath=shift;
     my  $collection=shift;
     my  $hostName=shift;
      my $type='coord';
     # $filePath='db/dbnsfp/dbNSFP2.8_variant.chr22';
     $filepath=
      my $date = POSIX::strftime("%d%m%y", localtime( ( stat $filePath )[9]));
#      my $md5sum=`md5sum $filePath`;
#      $md5sum=~s/ .*//g;
#      chomp($md5sum);
      my $wc=`wc -l $filePath`;
      $wc=~s/ .*//g;
      my ($conn,$db,$collections)=initMongo->new('denote','collections',$hostName);
      $collections->update({'md5sum'=>$md5sum},{'date'=>$date,'type'=>$type,'collection'=>$collection},{'upsert'=>1,'safe'=>1});
      return;
}
sub addHead{
	$head[0]='GCh38chr';
	$head[1]='GCh38pos';
	$head[2]='ref';
	$head[3]='alt';
	$head[4]='aaref';
	$head[5]='aaalt';
	$head[6]='dbsnp144';
	$head[7]='chr';
	$head[8]='hg19pos';
	$head[9]='hg18chr';
	$head[10]='hg18pos';
	$head[11]='genename';
	$head[12]='cdsStrand';
	$head[13]='refcodon';
	$head[14]='codonpos';
	$head[15]='codon_degeneracy';
	$head[16]='Ancestral_allele';
	$head[17]='AltaiNeandertal';
	$head[18]='Denisova';
	$head[19]='Ensembl_geneid';
	$head[20]='Ensembl_transcriptid';
	$head[21]='Ensembl_proteinid';
	$head[22]='aapos';
	$head[23]='sift_score';
	$head[24]='sift_converted_rankscore';
	$head[25]='sift';
	$head[26]='Uniprot_acc_Polyphen2';
	$head[27]='Uniprot_id_Polyphen2';
	$head[28]='Uniprot_aapos_Polyphen2';
	$head[29]='Polyphen2_HDIV_score';
	$head[30]='Polyphen2_HDIV_rankscore';
	$head[31]='Polyphen2_HDIV';
	$head[32]='Polyphen2_HVAR_score';
	$head[33]='Polyphen2_HVAR_rankscore';
	$head[34]='Polyphen2_HVAR';
	$head[35]='LRT_score';
	$head[36]='LRT_converted_rankscore';
	$head[37]='LRT';
	$head[38]='LRT_Omega';
	$head[39]='MutationTaster_score';
	$head[40]='MutationTaster_converted_rankscore';
	$head[41]='MutationTaster';
	$head[42]='MutationTaster_model';
	$head[43]='MutationTaster_AAE';
	$head[44]='Uniprot_id_MutationAssessor';
	$head[45]='Uniprot_variant_MutationAssessor';
	$head[46]='MutationAssessor_score';
	$head[47]='MutationAssessor_rankscore';
	$head[48]='MutationAssessor';
	$head[49]='FATHMM_score';
	$head[50]='FATHMM_converted_rankscore';
	$head[51]='FATHMM';
	$head[52]='PROVEAN_score';
	$head[53]='PROVEAN_converted_rankscore';
	$head[54]='PROVEAN';
	$head[55]='Transcript_id_VEST3';
	$head[56]='Transcript_var_VEST3';
	$head[57]='VEST3';
	$head[58]='VEST3_rankscore';
	$head[59]='MetaSVM_score';
	$head[60]='MetaSVM_rankscore';
	$head[61]='MetaSVM_pred';
	$head[62]='MetaLR_score';
	$head[63]='MetaLR_rankscore';
	$head[64]='MetaLR_pred';
	$head[65]='Reliability_index';
	$head[66]='CADD_raw';
	$head[67]='CADD_raw_rankscore';
	$head[68]='CADD_phred';
	$head[69]='DANN_score';
	$head[70]='DANN_rankscore';
	$head[71]='fathmm-MKL_coding_score';
	$head[72]='fathmm-MKL_coding_rankscore';
	$head[73]='fathmm-MKL_coding_pred';
	$head[74]='fathmm-MKL_coding_group';
	$head[75]='Eigen-raw';
	$head[76]='Eigen-phred';
	$head[77]='Eigen-raw_rankscore';
	$head[78]='Eigen-PC-raw';
	$head[79]='Eigen-PC-raw_rankscore';
	$head[80]='GenoCanyon_score';
	$head[81]='GenoCanyon_score_rankscore';
	$head[82]='integrated_fitCons_score';
	$head[83]='integrated_fitCons_score_rankscore';
	$head[84]='integrated_confidence_value';
	$head[85]='GM12878_fitCons_score';
	$head[86]='GM12878_fitCons_score_rankscore';
	$head[87]='GM12878_confidence_value';
	$head[88]='H1-hESC_fitCons_score';
	$head[89]='H1-hESC_fitCons_score_rankscore';
	$head[90]='H1-hESC_confidence_value';
	$head[91]='HUVEC_fitCons_score';
	$head[92]='HUVEC_fitCons_score_rankscore';
	$head[93]='HUVEC_confidence_value';
	$head[94]='GERP++_NR';
	$head[95]='GERP++_RS';
	$head[96]='GERP++_RS_rankscore';
	$head[97]='phyloP100way_vertebrate';
	$head[98]='phyloP100way_vertebrate_rankscore';
	$head[99]='phyloP20way_mammalian';
	$head[100]='phyloP20way_mammalian_rankscore';
	$head[101]='phastCons100way_vertebrate';
	$head[102]='phastCons100way_vertebrate_rankscore';
	$head[103]='phastCons20way_mammalian';
	$head[104]='phastCons20way_mammalian_rankscore';
	$head[105]='SiPhy_29way_pi';
	$head[106]='SiPhy_29way_logOdds';
	$head[107]='SiPhy_29way_logOdds_rankscore';
	$head[108]='1000Gp3_AC';
	$head[109]='1000Gp3_AF';
	$head[110]='1000Gp3_AFR_AC';
	$head[111]='1000Gp3_AFR_AF';
	$head[112]='1000Gp3_EUR_AC';
	$head[113]='1000Gp3_EUR_AF';
	$head[114]='1000Gp3_AMR_AC';
	$head[115]='1000Gp3_AMR_AF';
	$head[116]='1000Gp3_EAS_AC';
	$head[117]='1000Gp3_EAS_AF';
	$head[118]='1000Gp3_SAS_AC';
	$head[119]='1000Gp3_SAS_AF';
	$head[120]='TWINSUK_AC';
	$head[121]='TWINSUK_AF';
	$head[122]='ALSPAC_AC';
	$head[123]='ALSPAC_AF';
	$head[124]='ESP6500_AA_AC';
	$head[125]='ESP6500_AA_AF';
	$head[126]='ESP6500_EA_AC';
	$head[127]='ESP6500_EA_AF';
	$head[128]='ExAC_AC';
	$head[129]='ExAC_AF';
	$head[130]='ExAC_Adj_AC';
	$head[131]='ExAC_Adj_AF';
	$head[132]='ExAC_AFR_AC';
	$head[133]='ExAC_AFR_AF';
	$head[134]='ExAC_AMR_AC';
	$head[135]='ExAC_AMR_AF';
	$head[136]='ExAC_EAS_AC';
	$head[137]='ExAC_EAS_AF';
	$head[138]='ExAC_FIN_AC';
	$head[139]='ExAC_FIN_AF';
	$head[140]='ExAC_NFE_AC';
	$head[141]='ExAC_NFE_AF';
	$head[142]='ExAC_SAS_AC';
	$head[143]='ExAC_SAS_AF';
	$head[144]='ExAC_nonTCGA_AC';
	$head[145]='ExAC_nonTCGA_AF';
	$head[146]='ExAC_nonTCGA_Adj_AC';
	$head[147]='ExAC_nonTCGA_Adj_AF';
	$head[148]='ExAC_nonTCGA_AFR_AC';
	$head[149]='ExAC_nonTCGA_AFR_AF';
	$head[150]='ExAC_nonTCGA_AMR_AC';
	$head[151]='ExAC_nonTCGA_AMR_AF';
	$head[152]='ExAC_nonTCGA_EAS_AC';
	$head[153]='ExAC_nonTCGA_EAS_AF';
	$head[154]='ExAC_nonTCGA_FIN_AC';
	$head[155]='ExAC_nonTCGA_FIN_AF';
	$head[156]='ExAC_nonTCGA_NFE_AC';
	$head[157]='ExAC_nonTCGA_NFE_AF';
	$head[158]='ExAC_nonTCGA_SAS_AC';
	$head[159]='ExAC_nonTCGA_SAS_AF';
	$head[160]='ExAC_nonpsych_AC';
	$head[161]='ExAC_nonpsych_AF';
	$head[162]='ExAC_nonpsych_Adj_AC';
	$head[163]='ExAC_nonpsych_Adj_AF';
	$head[164]='ExAC_nonpsych_AFR_AC';
	$head[165]='ExAC_nonpsych_AFR_AF';
	$head[166]='ExAC_nonpsych_AMR_AC';
	$head[167]='ExAC_nonpsych_AMR_AF';
	$head[168]='ExAC_nonpsych_EAS_AC';
	$head[169]='ExAC_nonpsych_EAS_AF';
	$head[170]='ExAC_nonpsych_FIN_AC';
	$head[171]='ExAC_nonpsych_FIN_AF';
	$head[172]='ExAC_nonpsych_NFE_AC';
	$head[173]='ExAC_nonpsych_NFE_AF';
	$head[174]='ExAC_nonpsych_SAS_AC';
	$head[175]='ExAC_nonpsych_SAS_AF';
	$head[176]='clinvar_rs';
	$head[177]='clinvar_clnsig';
	$head[178]='clinvar_trait';
	$head[179]='clinvar_golden_stars';
	$head[180]='Interpro_domain';
	$head[181]='GTEx_V6_gene';
	$head[182]='GTEx_V6_tissue';
}
1;

