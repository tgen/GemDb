#!/usr/bin/perl -w
#dbnsfp2.4
package insertDbnsfp35;
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
      my @FILES=`ls -1 db/dbnsfp3.5/dbNSFP3*_variant.chr*`;
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
			  foreach $i (11,14,19,20,21,22,25,31,34,41,48,51,54,56,61,64,83,238,239,240,241,242,243,244) {
				if ($temp[$i] ne ".") {
					$hash{$head[$i]}=$temp[$i];
				}
			  }
			  for $i (78,79,87,91,93,96,99,102,105,106,108,110,112,117) {
				  if ($temp[$i] ne "." ) {
					 $temp[$i]=~s/\;.*//g;
					 $temp[$i]=$temp[$i]*1.00;
					 MongoDB::force_double($temp[$i]);
					 $hash{$head[$i]}=$temp[$i];
				   }  		  
			  }
			  for $i (120,121,122,124,126,128,130,132,134,136,138,140,142,143,144,146,148,150,152,154,156,158,159,160,162,164,166,168,170,172,174,175,176,178,180,182,184,186,189,190,191,192,195,198,201,204,207,210,213,216,217,218,219,222,225,228,231,234,237){
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
$head[6]='rs_dbSNP150';
$head[7]='chr';
$head[8]='hg19pos';
$head[9]='hg18chr';
$head[10]='hg18pos';
$head[11]='genename';
$head[12]='cds_strand';
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
$head[23]='SIFT_score';
$head[24]='SIFT_converted_rankscore';
$head[25]='SIFT_pred';
$head[26]='Uniprot_acc_Polyphen2';
$head[27]='Uniprot_id_Polyphen2';
$head[28]='Uniprot_aapos_Polyphen2';
$head[29]='Polyphen2_HDIV_score';
$head[30]='Polyphen2_HDIV_rankscore';
$head[31]='Polyphen2_HDIV_pred';
$head[32]='Polyphen2_HVAR_score';
$head[33]='Polyphen2_HVAR_rankscore';
$head[34]='Polyphen2_HVAR_pred';
$head[35]='LRT_score';
$head[36]='LRT_converted_rankscore';
$head[37]='LRT_pred';
$head[38]='LRT_Omega';
$head[39]='MutationTaster_score';
$head[40]='MutationTaster_converted_rankscore';
$head[41]='MutationTaster_pred';
$head[42]='MutationTaster_model';
$head[43]='MutationTaster_AAE';
$head[44]='MutationAssessor_UniprotID';
$head[45]='MutationAssessor_variant';
$head[46]='MutationAssessor_score';
$head[47]='MutationAssessor_score_rankscore';
$head[48]='MutationAssessor_pred';
$head[49]='FATHMM_score';
$head[50]='FATHMM_converted_rankscore';
$head[51]='FATHMM_pred';
$head[52]='PROVEAN_score';
$head[53]='PROVEAN_converted_rankscore';
$head[54]='PROVEAN_pred';
$head[55]='Transcript_id_VEST3';
$head[56]='Transcript_var_VEST3';
$head[57]='VEST3_score';
$head[58]='VEST3_rankscore';
$head[59]='MetaSVM_score';
$head[60]='MetaSVM_rankscore';
$head[61]='MetaSVM_pred';
$head[62]='MetaLR_score';
$head[63]='MetaLR_rankscore';
$head[64]='MetaLR_pred';
$head[65]='Reliability_index';
$head[66]='M-CAP_score';
$head[67]='M-CAP_rankscore';
$head[68]='M-CAP_pred';
$head[69]='REVEL_score';
$head[70]='REVEL_rankscore';
$head[71]='MutPred_score';
$head[72]='MutPred_rankscore';
$head[73]='MutPred_protID';
$head[74]='MutPred_AAchange';
$head[75]='MutPred_Top5features';
$head[76]='CADD_raw';
$head[77]='CADD_raw_rankscore';
$head[78]='CADD_phred';
$head[79]='DANN_score';
$head[80]='DANN_rankscore';
$head[81]='fathmm-MKL_coding_score';
$head[82]='fathmm-MKL_coding_rankscore';
$head[83]='fathmm-MKL_coding_pred';
$head[84]='fathmm-MKL_coding_group';
$head[85]='Eigen_coding_or_noncoding';
$head[86]='Eigen-raw';
$head[87]='Eigen-phred';
$head[88]='Eigen-PC-raw';
$head[89]='Eigen-PC-phred';
$head[90]='Eigen-PC-raw_rankscore';
$head[91]='GenoCanyon_score';
$head[92]='GenoCanyon_score_rankscore';
$head[93]='integrated_fitCons_score';
$head[94]='integrated_fitCons_score_rankscore';
$head[95]='integrated_confidence_value';
$head[96]='GM12878_fitCons_score';
$head[97]='GM12878_fitCons_score_rankscore';
$head[98]='GM12878_confidence_value';
$head[99]='H1-hESC_fitCons_score';
$head[100]='H1-hESC_fitCons_score_rankscore';
$head[101]='H1-hESC_confidence_value';
$head[102]='HUVEC_fitCons_score';
$head[103]='HUVEC_fitCons_score_rankscore';
$head[104]='HUVEC_confidence_value';
$head[105]='GERP++_NR';
$head[106]='GERP++_RS';
$head[107]='GERP++_RS_rankscore';
$head[108]='phyloP100way_vertebrate';
$head[109]='phyloP100way_vertebrate_rankscore';
$head[110]='phyloP20way_mammalian';
$head[111]='phyloP20way_mammalian_rankscore';
$head[112]='phastCons100way_vertebrate';
$head[113]='phastCons100way_vertebrate_rankscore';
$head[114]='phastCons20way_mammalian';
$head[115]='phastCons20way_mammalian_rankscore';
$head[116]='SiPhy_29way_pi';
$head[117]='SiPhy_29way_logOdds';
$head[118]='SiPhy_29way_logOdds_rankscore';
$head[119]='1000Gp3_AC';
$head[120]='1000Gp3_AF';
$head[121]='1000Gp3_AFR_AC';
$head[122]='1000Gp3_AFR_AF';
$head[123]='1000Gp3_EUR_AC';
$head[124]='1000Gp3_EUR_AF';
$head[125]='1000Gp3_AMR_AC';
$head[126]='1000Gp3_AMR_AF';
$head[127]='1000Gp3_EAS_AC';
$head[128]='1000Gp3_EAS_AF';
$head[129]='1000Gp3_SAS_AC';
$head[130]='1000Gp3_SAS_AF';
$head[131]='TWINSUK_AC';
$head[132]='TWINSUK_AF';
$head[133]='ALSPAC_AC';
$head[134]='ALSPAC_AF';
$head[135]='ESP6500_AA_AC';
$head[136]='ESP6500_AA_AF';
$head[137]='ESP6500_EA_AC';
$head[138]='ESP6500_EA_AF';
$head[139]='ExAC_AC';
$head[140]='ExAC_AF';
$head[141]='ExAC_Adj_AC';
$head[142]='ExAC_Adj_AF';
$head[143]='ExAC_AFR_AC';
$head[144]='ExAC_AFR_AF';
$head[145]='ExAC_AMR_AC';
$head[146]='ExAC_AMR_AF';
$head[147]='ExAC_EAS_AC';
$head[148]='ExAC_EAS_AF';
$head[149]='ExAC_FIN_AC';
$head[150]='ExAC_FIN_AF';
$head[151]='ExAC_NFE_AC';
$head[152]='ExAC_NFE_AF';
$head[153]='ExAC_SAS_AC';
$head[154]='ExAC_SAS_AF';
$head[155]='ExAC_nonTCGA_AC';
$head[156]='ExAC_nonTCGA_AF';
$head[157]='ExAC_nonTCGA_Adj_AC';
$head[158]='ExAC_nonTCGA_Adj_AF';
$head[159]='ExAC_nonTCGA_AFR_AC';
$head[160]='ExAC_nonTCGA_AFR_AF';
$head[161]='ExAC_nonTCGA_AMR_AC';
$head[162]='ExAC_nonTCGA_AMR_AF';
$head[163]='ExAC_nonTCGA_EAS_AC';
$head[164]='ExAC_nonTCGA_EAS_AF';
$head[165]='ExAC_nonTCGA_FIN_AC';
$head[166]='ExAC_nonTCGA_FIN_AF';
$head[167]='ExAC_nonTCGA_NFE_AC';
$head[168]='ExAC_nonTCGA_NFE_AF';
$head[169]='ExAC_nonTCGA_SAS_AC';
$head[170]='ExAC_nonTCGA_SAS_AF';
$head[171]='ExAC_nonpsych_AC';
$head[172]='ExAC_nonpsych_AF';
$head[173]='ExAC_nonpsych_Adj_AC';
$head[174]='ExAC_nonpsych_Adj_AF';
$head[175]='ExAC_nonpsych_AFR_AC';
$head[176]='ExAC_nonpsych_AFR_AF';
$head[177]='ExAC_nonpsych_AMR_AC';
$head[178]='ExAC_nonpsych_AMR_AF';
$head[179]='ExAC_nonpsych_EAS_AC';
$head[180]='ExAC_nonpsych_EAS_AF';
$head[181]='ExAC_nonpsych_FIN_AC';
$head[182]='ExAC_nonpsych_FIN_AF';
$head[183]='ExAC_nonpsych_NFE_AC';
$head[184]='ExAC_nonpsych_NFE_AF';
$head[185]='ExAC_nonpsych_SAS_AC';
$head[186]='ExAC_nonpsych_SAS_AF';
$head[187]='gnomAD_exomes_AC';
$head[188]='gnomAD_exomes_AN';
$head[189]='gnomAD_exomes_AF';
$head[190]='gnomAD_exomes_AFR_AC';
$head[191]='gnomAD_exomes_AFR_AN';
$head[192]='gnomAD_exomes_AFR_AF';
$head[193]='gnomAD_exomes_AMR_AC';
$head[194]='gnomAD_exomes_AMR_AN';
$head[195]='gnomAD_exomes_AMR_AF';
$head[196]='gnomAD_exomes_ASJ_AC';
$head[197]='gnomAD_exomes_ASJ_AN';
$head[198]='gnomAD_exomes_ASJ_AF';
$head[199]='gnomAD_exomes_EAS_AC';
$head[200]='gnomAD_exomes_EAS_AN';
$head[201]='gnomAD_exomes_EAS_AF';
$head[202]='gnomAD_exomes_FIN_AC';
$head[203]='gnomAD_exomes_FIN_AN';
$head[204]='gnomAD_exomes_FIN_AF';
$head[205]='gnomAD_exomes_NFE_AC';
$head[206]='gnomAD_exomes_NFE_AN';
$head[207]='gnomAD_exomes_NFE_AF';
$head[208]='gnomAD_exomes_SAS_AC';
$head[209]='gnomAD_exomes_SAS_AN';
$head[210]='gnomAD_exomes_SAS_AF';
$head[211]='gnomAD_exomes_OTH_AC';
$head[212]='gnomAD_exomes_OTH_AN';
$head[213]='gnomAD_exomes_OTH_AF';
$head[214]='gnomAD_genomes_AC';
$head[215]='gnomAD_genomes_AN';
$head[216]='gnomAD_genomes_AF';
$head[217]='gnomAD_genomes_AFR_AC';
$head[218]='gnomAD_genomes_AFR_AN';
$head[219]='gnomAD_genomes_AFR_AF';
$head[220]='gnomAD_genomes_AMR_AC';
$head[221]='gnomAD_genomes_AMR_AN';
$head[222]='gnomAD_genomes_AMR_AF';
$head[223]='gnomAD_genomes_ASJ_AC';
$head[224]='gnomAD_genomes_ASJ_AN';
$head[225]='gnomAD_genomes_ASJ_AF';
$head[226]='gnomAD_genomes_EAS_AC';
$head[227]='gnomAD_genomes_EAS_AN';
$head[228]='gnomAD_genomes_EAS_AF';
$head[229]='gnomAD_genomes_FIN_AC';
$head[230]='gnomAD_genomes_FIN_AN';
$head[231]='gnomAD_genomes_FIN_AF';
$head[232]='gnomAD_genomes_NFE_AC';
$head[233]='gnomAD_genomes_NFE_AN';
$head[234]='gnomAD_genomes_NFE_AF';
$head[235]='gnomAD_genomes_OTH_AC';
$head[236]='gnomAD_genomes_OTH_AN';
$head[237]='gnomAD_genomes_OTH_AF';
$head[238]='clinvar_rs';
$head[239]='clinvar_clnsig';
$head[240]='clinvar_trait';
$head[241]='clinvar_golden_stars';
$head[242]='Interpro_domain';
$head[243]='GTEx_V6p_gene';
$head[244]='GTEx_V6p_tissue';
}

1;

