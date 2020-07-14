#!/usr/bin/perl -w
package insertDbnsfp40_b38;

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
      my @FILES=`ls -1 db/dbnsfp/dbnsfp4.0_academic/dbNSFP*chr*`;
	  foreach $file (@FILES) {
		  my $done={};
		  chomp($file);
		  print "-dbnsfp opening $file\n";
		  open(FILE,"gunzip -c $file |") or die "Can't open $file\n";
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
			  $hash{'coord'}=join(':',"chr" . $temp[0],$temp[1],$temp[2],$temp[3]);                  
			  if (exists($done->{$hash{'coord'}})){next LOOP}               
			  

			  # String 
			  for $i (6){
			  #for $i ( 6, 13, 14, 15, 16, 17, 41, 44, 47, 54, 62, 65, 70, 363, 364, 365, 366, 369, 370){
				if ($temp[$i] !~ m/^[.].*/) {
					$temp[$i]=~s/\;.*//g;
					$hash{$head[$i]}=$temp[$i];
				}
			  }
			  
			  # Double variable
			  #for $i (500){
			  #for $i (39, 40, 42, 43, 45, 46, 101, 102, 103, 104, 115, 135, 136, 137,155,156,157,148,159,160,161,162,163,164,165,166,226, 227, 228, 229, 299, 300, 301, 302 ){
	  
			#	if ($temp[$i] !~ m/^[.].*/ || $temp[$i] ne "0") {
#if ($temp[$i] ne "." || $temp[$i] ne "./." || $temp[$i] ne "././." || $temp[$i] ne "0" ) {
			#		 $temp[$i]=~s/\;.*//g;
			#		 $temp[$i]=$temp[$i]*1.00;
			#		 MongoDB::force_double($temp[$i]);
			#		 if ($temp[$i] != 0.0){
			#		 	$hash{$head[$i]}=$temp[$i];
			#		 }
			#	   }  		  
			#  }
			  # Integer values
			 # for $i (500){
			  #for $i (7,8) {
			#	 if ($temp[$i] ne ".") {
			#		 $hash{$head[$i]}=int($temp[$i]);
			#	 }
			 # }			  
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
$head[0]='chr';  #int
$head[1]='pos(1-based)';#int
$head[2]='ref';#string
$head[3]='alt';#string
$head[4]='aaref';#string
$head[5]='aaalt';
$head[6]='rs_dbSNP151';
$head[7]='hg19_chr';
$head[8]='hg19_pos(1-based)';
$head[9]='hg18_chr';
$head[10]='hg18_pos(1-based)';
$head[11]='aapos';
$head[12]='genename';
$head[13]='Ensembl_geneid';
$head[14]='Ensembl_transcriptid';
$head[15]='Ensembl_proteinid';
$head[16]='Uniprot_acc';
$head[17]='Uniprot_entry';
$head[18]='HGVSc_ANNOVAR';
$head[19]='HGVSp_ANNOVAR';
$head[20]='HGVSc_snpEff';
$head[21]='HGVSp_snpEff';
$head[22]='HGVSc_VEP';
$head[23]='HGVSp_VEP';
$head[24]='APPRIS';
$head[25]='GENCODE_basic';
$head[26]='TSL';
$head[27]='VEP_canonical';
$head[28]='cds_strand';
$head[29]='refcodon';
$head[30]='codonpos';
$head[31]='codon_degeneracy';
$head[32]='Ancestral_allele';
$head[33]='AltaiNeandertal';
$head[34]='Denisova';
$head[35]='VindijiaNeandertal';
$head[36]='SIFT_score';
$head[37]='SIFT_converted_rankscore';
$head[38]='SIFT_pred';
$head[39]='SIFT4G_score';
$head[40]='SIFT4G_converted_rankscore';
$head[41]='SIFT4G_pred';
$head[42]='Polyphen2_HDIV_score';
$head[43]='Polyphen2_HDIV_rankscore';
$head[44]='Polyphen2_HDIV_pred';
$head[45]='Polyphen2_HVAR_score';
$head[46]='Polyphen2_HVAR_rankscore';
$head[47]='Polyphen2_HVAR_pred';
$head[48]='LRT_score';
$head[49]='LRT_converted_rankscore';
$head[50]='LRT_pred';
$head[51]='LRT_Omega';
$head[52]='MutationTaster_score';
$head[53]='MutationTaster_converted_rankscore';
$head[54]='MutationTaster_pred';
$head[55]='MutationTaster_model';
$head[56]='MutationTaster_AAE';
$head[57]='MutationAssessor_score';
$head[58]='MutationAssessor_rankscore';
$head[59]='MutationAssessor_pred';
$head[60]='FATHMM_score';
$head[61]='FATHMM_converted_rankscore';
$head[62]='FATHMM_pred';
$head[63]='PROVEAN_score';
$head[64]='PROVEAN_converted_rankscore';
$head[65]='PROVEAN_pred';
$head[66]='VEST4_score';
$head[67]='VEST4_rankscore';
$head[68]='MetaSVM_score';
$head[69]='MetaSVM_rankscore';
$head[70]='MetaSVM_pred';
$head[71]='MetaLR_score';
$head[72]='MetaLR_rankscore';
$head[73]='MetaLR_pred';
$head[74]='Reliability_index';
$head[75]='M-CAP_score';
$head[76]='M-CAP_rankscore';
$head[77]='M-CAP_pred';
$head[78]='REVEL_score';
$head[79]='REVEL_rankscore';
$head[80]='MutPred_score';
$head[81]='MutPred_rankscore';
$head[82]='MutPred_protID';
$head[83]='MutPred_AAchange';
$head[84]='MutPred_Top5features';
$head[85]='MVP_score';
$head[86]='MVP_rankscore';
$head[87]='MPC_score';
$head[88]='MPC_rankscore';
$head[89]='PrimateAI_score';
$head[90]='PrimateAI_rankscore';
$head[91]='PrimateAI_pred';
$head[92]='DEOGEN2_score';
$head[93]='DEOGEN2_rankscore';
$head[94]='DEOGEN2_pred';
$head[95]='Aloft_Fraction_transcripts_affected';
$head[96]='Aloft_prob_Tolerant';
$head[97]='Aloft_prob_Recessive';
$head[98]='Aloft_prob_Dominant';
$head[99]='Aloft_pred';
$head[100]='Aloft_Confidence';
$head[101]='CADD_raw';
$head[102]='CADD_raw_rankscore';
$head[103]='CADD_phred';
$head[104]='DANN_score';
$head[105]='DANN_rankscore';
$head[106]='fathmm-MKL_coding_score';
$head[107]='fathmm-MKL_coding_rankscore';
$head[108]='fathmm-MKL_coding_pred';
$head[109]='fathmm-MKL_coding_group';
$head[110]='fathmm-XF_coding_score';
$head[111]='fathmm-XF_coding_rankscore';
$head[112]='fathmm-XF_coding_pred';
$head[113]='Eigen-raw_coding';
$head[114]='Eigen-raw_coding_rankscore';
$head[115]='Eigen-phred_coding';
$head[116]='Eigen-PC-raw_coding';
$head[117]='Eigen-PC-raw_coding_rankscore';
$head[118]='Eigen-PC-phred_coding';
$head[119]='GenoCanyon_score';
$head[120]='GenoCanyon_score_rankscore';
$head[121]='integrated_fitCons_score';
$head[122]='integrated_fitCons_rankscore';
$head[123]='integrated_confidence_value';
$head[124]='GM12878_fitCons_score';
$head[125]='GM12878_fitCons_rankscore';
$head[126]='GM12878_confidence_value';
$head[127]='H1-hESC_fitCons_score';
$head[128]='H1-hESC_fitCons_rankscore';
$head[129]='H1-hESC_confidence_value';
$head[130]='HUVEC_fitCons_score';
$head[131]='HUVEC_fitCons_rankscore';
$head[132]='HUVEC_confidence_value';
$head[133]='LINSIGHT';
$head[134]='LINSIGHT_rankscore';
$head[135]='GERP++_NR';
$head[136]='GERP++_RS';
$head[137]='GERP++_RS_rankscore';
$head[138]='phyloP100way_vertebrate';
$head[139]='phyloP100way_vertebrate_rankscore';
$head[140]='phyloP30way_mammalian';
$head[141]='phyloP30way_mammalian_rankscore';
$head[142]='phyloP17way_primate';
$head[143]='phyloP17way_primate_rankscore';
$head[144]='phastCons100way_vertebrate';
$head[145]='phastCons100way_vertebrate_rankscore';
$head[146]='phastCons30way_mammalian';
$head[147]='phastCons30way_mammalian_rankscore';
$head[148]='phastCons17way_primate';
$head[149]='phastCons17way_primate_rankscore';
$head[150]='SiPhy_29way_pi';
$head[151]='SiPhy_29way_logOdds';
$head[152]='SiPhy_29way_logOdds_rankscore';
$head[153]='bStatistic';
$head[154]='bStatistic_rankscore';
$head[155]='1000Gp3_AC';
$head[156]='1000Gp3_AF';
$head[157]='1000Gp3_AFR_AC';
$head[158]='1000Gp3_AFR_AF';
$head[159]='1000Gp3_EUR_AC';
$head[160]='1000Gp3_EUR_AF';
$head[161]='1000Gp3_AMR_AC';
$head[162]='1000Gp3_AMR_AF';
$head[163]='1000Gp3_EAS_AC';
$head[164]='1000Gp3_EAS_AF';
$head[165]='1000Gp3_SAS_AC';
$head[166]='1000Gp3_SAS_AF';
$head[167]='TWINSUK_AC';
$head[168]='TWINSUK_AF';
$head[169]='ALSPAC_AC';
$head[170]='ALSPAC_AF';
$head[171]='UK10K_AC';
$head[172]='UK10K_AF';
$head[173]='ESP6500_AA_AC';
$head[174]='ESP6500_AA_AF';
$head[175]='ESP6500_EA_AC';
$head[176]='ESP6500_EA_AF';
$head[177]='ExAC_AC';
$head[178]='ExAC_AF';
$head[179]='ExAC_Adj_AC';
$head[180]='ExAC_Adj_AF';
$head[181]='ExAC_AFR_AC';
$head[182]='ExAC_AFR_AF';
$head[183]='ExAC_AMR_AC';
$head[184]='ExAC_AMR_AF';
$head[185]='ExAC_EAS_AC';
$head[186]='ExAC_EAS_AF';
$head[187]='ExAC_FIN_AC';
$head[188]='ExAC_FIN_AF';
$head[189]='ExAC_NFE_AC';
$head[190]='ExAC_NFE_AF';
$head[191]='ExAC_SAS_AC';
$head[192]='ExAC_SAS_AF';
$head[193]='ExAC_nonTCGA_AC';
$head[194]='ExAC_nonTCGA_AF';
$head[195]='ExAC_nonTCGA_Adj_AC';
$head[196]='ExAC_nonTCGA_Adj_AF';
$head[197]='ExAC_nonTCGA_AFR_AC';
$head[198]='ExAC_nonTCGA_AFR_AF';
$head[199]='ExAC_nonTCGA_AMR_AC';
$head[200]='ExAC_nonTCGA_AMR_AF';
$head[201]='ExAC_nonTCGA_EAS_AC';
$head[202]='ExAC_nonTCGA_EAS_AF';
$head[203]='ExAC_nonTCGA_FIN_AC';
$head[204]='ExAC_nonTCGA_FIN_AF';
$head[205]='ExAC_nonTCGA_NFE_AC';
$head[206]='ExAC_nonTCGA_NFE_AF';
$head[207]='ExAC_nonTCGA_SAS_AC';
$head[208]='ExAC_nonTCGA_SAS_AF';
$head[209]='ExAC_nonpsych_AC';
$head[210]='ExAC_nonpsych_AF';
$head[211]='ExAC_nonpsych_Adj_AC';
$head[212]='ExAC_nonpsych_Adj_AF';
$head[213]='ExAC_nonpsych_AFR_AC';
$head[214]='ExAC_nonpsych_AFR_AF';
$head[215]='ExAC_nonpsych_AMR_AC';
$head[216]='ExAC_nonpsych_AMR_AF';
$head[217]='ExAC_nonpsych_EAS_AC';
$head[218]='ExAC_nonpsych_EAS_AF';
$head[219]='ExAC_nonpsych_FIN_AC';
$head[220]='ExAC_nonpsych_FIN_AF';
$head[221]='ExAC_nonpsych_NFE_AC';
$head[222]='ExAC_nonpsych_NFE_AF';
$head[223]='ExAC_nonpsych_SAS_AC';
$head[224]='ExAC_nonpsych_SAS_AF';
$head[225]='gnomAD_exomes_flag';
$head[226]='gnomAD_exomes_AC';
$head[227]='gnomAD_exomes_AN';
$head[228]='gnomAD_exomes_AF';
$head[229]='gnomAD_exomes_nhomalt';
$head[230]='gnomAD_exomes_AFR_AC';
$head[231]='gnomAD_exomes_AFR_AN';
$head[232]='gnomAD_exomes_AFR_AF';
$head[233]='gnomAD_exomes_AFR_nhomalt';
$head[234]='gnomAD_exomes_AMR_AC';
$head[235]='gnomAD_exomes_AMR_AN';
$head[236]='gnomAD_exomes_AMR_AF';
$head[237]='gnomAD_exomes_AMR_nhomalt';
$head[238]='gnomAD_exomes_ASJ_AC';
$head[239]='gnomAD_exomes_ASJ_AN';
$head[240]='gnomAD_exomes_ASJ_AF';
$head[241]='gnomAD_exomes_ASJ_nhomalt';
$head[242]='gnomAD_exomes_EAS_AC';
$head[243]='gnomAD_exomes_EAS_AN';
$head[244]='gnomAD_exomes_EAS_AF';
$head[245]='gnomAD_exomes_EAS_nhomalt';
$head[246]='gnomAD_exomes_FIN_AC';
$head[247]='gnomAD_exomes_FIN_AN';
$head[248]='gnomAD_exomes_FIN_AF';
$head[249]='gnomAD_exomes_FIN_nhomalt';
$head[250]='gnomAD_exomes_NFE_AC';
$head[251]='gnomAD_exomes_NFE_AN';
$head[252]='gnomAD_exomes_NFE_AF';
$head[253]='gnomAD_exomes_NFE_nhomalt';
$head[254]='gnomAD_exomes_SAS_AC';
$head[255]='gnomAD_exomes_SAS_AN';
$head[256]='gnomAD_exomes_SAS_AF';
$head[257]='gnomAD_exomes_SAS_nhomalt';
$head[258]='gnomAD_exomes_POPMAX_AC';
$head[259]='gnomAD_exomes_POPMAX_AN';
$head[260]='gnomAD_exomes_POPMAX_AF';
$head[261]='gnomAD_exomes_POPMAX_nhomalt';
$head[262]='gnomAD_exomes_controls_AC';
$head[263]='gnomAD_exomes_controls_AN';
$head[264]='gnomAD_exomes_controls_AF';
$head[265]='gnomAD_exomes_controls_nhomalt';
$head[266]='gnomAD_exomes_controls_AFR_AC';
$head[267]='gnomAD_exomes_controls_AFR_AN';
$head[268]='gnomAD_exomes_controls_AFR_AF';
$head[269]='gnomAD_exomes_controls_AFR_nhomalt';
$head[270]='gnomAD_exomes_controls_AMR_AC';
$head[271]='gnomAD_exomes_controls_AMR_AN';
$head[272]='gnomAD_exomes_controls_AMR_AF';
$head[273]='gnomAD_exomes_controls_AMR_nhomalt';
$head[274]='gnomAD_exomes_controls_ASJ_AC';
$head[275]='gnomAD_exomes_controls_ASJ_AN';
$head[276]='gnomAD_exomes_controls_ASJ_AF';
$head[277]='gnomAD_exomes_controls_ASJ_nhomalt';
$head[278]='gnomAD_exomes_controls_EAS_AC';
$head[279]='gnomAD_exomes_controls_EAS_AN';
$head[280]='gnomAD_exomes_controls_EAS_AF';
$head[281]='gnomAD_exomes_controls_EAS_nhomalt';
$head[282]='gnomAD_exomes_controls_FIN_AC';
$head[283]='gnomAD_exomes_controls_FIN_AN';
$head[284]='gnomAD_exomes_controls_FIN_AF';
$head[285]='gnomAD_exomes_controls_FIN_nhomalt';
$head[286]='gnomAD_exomes_controls_NFE_AC';
$head[287]='gnomAD_exomes_controls_NFE_AN';
$head[288]='gnomAD_exomes_controls_NFE_AF';
$head[289]='gnomAD_exomes_controls_NFE_nhomalt';
$head[290]='gnomAD_exomes_controls_SAS_AC';
$head[291]='gnomAD_exomes_controls_SAS_AN';
$head[292]='gnomAD_exomes_controls_SAS_AF';
$head[293]='gnomAD_exomes_controls_SAS_nhomalt';
$head[294]='gnomAD_exomes_controls_POPMAX_AC';
$head[295]='gnomAD_exomes_controls_POPMAX_AN';
$head[296]='gnomAD_exomes_controls_POPMAX_AF';
$head[297]='gnomAD_exomes_controls_POPMAX_nhomalt';
$head[298]='gnomAD_genomes_flag';
$head[299]='gnomAD_genomes_AC';
$head[300]='gnomAD_genomes_AN';
$head[301]='gnomAD_genomes_AF';
$head[302]='gnomAD_genomes_nhomalt';
$head[303]='gnomAD_genomes_AFR_AC';
$head[304]='gnomAD_genomes_AFR_AN';
$head[305]='gnomAD_genomes_AFR_AF';
$head[306]='gnomAD_genomes_AFR_nhomalt';
$head[307]='gnomAD_genomes_AMR_AC';
$head[308]='gnomAD_genomes_AMR_AN';
$head[309]='gnomAD_genomes_AMR_AF';
$head[310]='gnomAD_genomes_AMR_nhomalt';
$head[311]='gnomAD_genomes_ASJ_AC';
$head[312]='gnomAD_genomes_ASJ_AN';
$head[313]='gnomAD_genomes_ASJ_AF';
$head[314]='gnomAD_genomes_ASJ_nhomalt';
$head[315]='gnomAD_genomes_EAS_AC';
$head[316]='gnomAD_genomes_EAS_AN';
$head[317]='gnomAD_genomes_EAS_AF';
$head[318]='gnomAD_genomes_EAS_nhomalt';
$head[319]='gnomAD_genomes_FIN_AC';
$head[320]='gnomAD_genomes_FIN_AN';
$head[321]='gnomAD_genomes_FIN_AF';
$head[322]='gnomAD_genomes_FIN_nhomalt';
$head[323]='gnomAD_genomes_NFE_AC';
$head[324]='gnomAD_genomes_NFE_AN';
$head[325]='gnomAD_genomes_NFE_AF';
$head[326]='gnomAD_genomes_NFE_nhomalt';
$head[327]='gnomAD_genomes_POPMAX_AC';
$head[328]='gnomAD_genomes_POPMAX_AN';
$head[329]='gnomAD_genomes_POPMAX_AF';
$head[330]='gnomAD_genomes_POPMAX_nhomalt';
$head[331]='gnomAD_genomes_controls_AC';
$head[332]='gnomAD_genomes_controls_AN';
$head[333]='gnomAD_genomes_controls_AF';
$head[334]='gnomAD_genomes_controls_nhomalt';
$head[335]='gnomAD_genomes_controls_AFR_AC';
$head[336]='gnomAD_genomes_controls_AFR_AN';
$head[337]='gnomAD_genomes_controls_AFR_AF';
$head[338]='gnomAD_genomes_controls_AFR_nhomalt';
$head[339]='gnomAD_genomes_controls_AMR_AC';
$head[340]='gnomAD_genomes_controls_AMR_AN';
$head[341]='gnomAD_genomes_controls_AMR_AF';
$head[342]='gnomAD_genomes_controls_AMR_nhomalt';
$head[343]='gnomAD_genomes_controls_ASJ_AC';
$head[344]='gnomAD_genomes_controls_ASJ_AN';
$head[345]='gnomAD_genomes_controls_ASJ_AF';
$head[346]='gnomAD_genomes_controls_ASJ_nhomalt';
$head[347]='gnomAD_genomes_controls_EAS_AC';
$head[348]='gnomAD_genomes_controls_EAS_AN';
$head[349]='gnomAD_genomes_controls_EAS_AF';
$head[350]='gnomAD_genomes_controls_EAS_nhomalt';
$head[351]='gnomAD_genomes_controls_FIN_AC';
$head[352]='gnomAD_genomes_controls_FIN_AN';
$head[353]='gnomAD_genomes_controls_FIN_AF';
$head[354]='gnomAD_genomes_controls_FIN_nhomalt';
$head[355]='gnomAD_genomes_controls_NFE_AC';
$head[356]='gnomAD_genomes_controls_NFE_AN';
$head[357]='gnomAD_genomes_controls_NFE_AF';
$head[358]='gnomAD_genomes_controls_NFE_nhomalt';
$head[359]='gnomAD_genomes_controls_POPMAX_AC';
$head[360]='gnomAD_genomes_controls_POPMAX_AN';
$head[361]='gnomAD_genomes_controls_POPMAX_AF';
$head[362]='gnomAD_genomes_controls_POPMAX_nhomalt';
$head[363]='clinvar_id';
$head[364]='clinvar_clnsig';
$head[365]='clinvar_trait';
$head[366]='clinvar_review';
$head[367]='clinvar_hgvs';
$head[368]='clinvar_var_source';
$head[369]='clinvar_MedGen_id';
$head[370]='clinvar_OMIM_id';
$head[371]='clinvar_Orphanet_id';
$head[372]='Interpro_domain';
$head[373]='GTEx_V7_gene';
$head[374]='GTEx_V7_tissue';
$head[375]='Geuvadis_eQTL_target_gene';
}

1;

