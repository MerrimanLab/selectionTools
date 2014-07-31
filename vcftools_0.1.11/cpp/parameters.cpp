/*
 * parameters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */

// Class for reading in, checking and storing user parameters
#include "parameters.h"

parameters::parameters(int argc, char *argv[])
{
	string tmp;
	for (int i=0; i<argc; i++)
	{
		tmp = argv[i];
		this->argv.push_back(tmp);
	}

	bcf_format = false;
	BED_exclude = false;
	BED_file = "";
	derived = false;
	diff_discordance_matrix = false;
	diff_file = "";
	diff_file_bcf = false;
	diff_file_compressed = false;
	diff_indv_discordance = false;
	diff_indv_map_file = "";
	diff_site_discordance = false;
	diff_switch_error = false;
	end_pos = numeric_limits<int>::max();
	exclude_positions_file = "";
	force_write_index = false;
	fst_window_size = -1;
	fst_window_step = -1;
	gatk = false;
	indv_exclude_file = "";
	indv_keep_file = "";
	invert_mask = false;
	keep_only_indels = false;
	recode_all_INFO = false;
	ld_bp_window_size = numeric_limits<int>::max();
	ld_snp_window_size = numeric_limits<int>::max();
	ld_bp_window_min = -1;
	ld_snp_window_min = -1;
	min_mac = -1;
	min_maf = -1.0;
	mask_file = "";
	max_alleles = numeric_limits<int>::max();
	max_genotype_depth = numeric_limits<int>::max();
	max_indv_mean_depth = numeric_limits<double>::max();
	max_mac = numeric_limits<int>::max();
	max_maf = numeric_limits<double>::max();
	max_mean_depth = numeric_limits<double>::max();
	max_missing_call_count = numeric_limits<int>::max();
	max_non_ref_ac = numeric_limits<int>::max();
	max_non_ref_af = numeric_limits<double>::max();
	max_N_indv = -1;
	min_alleles = -1;
	min_genotype_depth = -1;
	min_genotype_quality = -1.0;
	min_HWE_pvalue = -1.0;
	min_indv_call_rate = 0;
	min_indv_mean_depth = -1.0;
	min_interSNP_distance = -1;
	min_kept_mask_value = 0;
	min_mean_depth = -1.0;
	min_quality = -1.0;
	min_r2 = -1.0;
	min_site_call_rate = 0;
	min_non_ref_ac = -1;
	min_non_ref_af = -1.0;
	output_012_matrix = false;
	output_as_IMPUTE = false;
	output_as_ldhat_phased = false;
	output_as_ldhat_unphased = false;
	output_BEAGLE_genotype_likelihoods_GL = false;
	output_BEAGLE_genotype_likelihoods_PL = false;
	output_counts = false;
	output_filter_summary = false;
	output_filtered_sites = false;
	output_freq = false;
	output_geno_depth = false;
	output_geno_chisq = false;
	output_geno_rsq = false;
	output_hap_rsq = false;
	output_het = false;
	output_HWE = false;
	output_indel_hist = false;
	output_indv_depth = false;
	output_interchromosomal_hap_rsq = false;
	output_interchromosomal_geno_rsq = false;
	output_LROH = false;
	output_missingness = false;
	output_N_PCA_SNP_loadings = -1;
	output_PCA = false;
	output_prefix="out";
	output_relatedness = false;
	output_singletons = false;
	output_site_depth = false;
	output_site_mean_depth = false;
	output_site_pi=false;
	output_site_quality = false;
	output_SNP_density_bin_size = 0;
	output_Tajima_D_bin_size = 0;
	output_TsTv_bin_size = 0;
	output_TsTv_by_count = false;
	output_TsTv_by_qual = false;
	phased_only = false;
	PCA_no_normalisation = false;
	pi_window_size = 0;
	pi_window_step = 0;
	plink_output = false;
	plink_tped_output = false;
	positions_file = "";
	recode = false;
	recode_bcf = false;
	recode_to_stream = false;
	recode_bcf_to_stream = false;
	remove_all_filtered_genotypes = false;
	remove_all_filtered_sites = false;
	remove_indels = false;
	snps_to_exclude_file = "";
	snps_to_keep_file = "";
	start_pos = -1;
	suppress_allele_output = false;
	vcf_filename="";
	vcf_compressed = false;
}

void parameters::read_parameters()
{
	unsigned int i=1;
	string in_str;
	while (i<argv.size())
	{
		in_str = argv[i];
		if (in_str == "--vcf") { vcf_filename = get_arg(i+1); vcf_compressed = false; i++; } 				// VCF file to process
		else if (in_str == "--bcf") { vcf_filename = get_arg(i+1); bcf_format = true; vcf_compressed = false; i++; } 				// BCF file to process
		else if (in_str == "--012") output_012_matrix = true;							// Output as 0/1/2 matrix
		else if (in_str == "--BEAGLE-GL") { output_BEAGLE_genotype_likelihoods_GL = true; min_alleles=2; max_alleles=2; }	// Output as BEAGLE Genotype Likelihood format
		else if (in_str == "--BEAGLE-PL") { output_BEAGLE_genotype_likelihoods_PL = true; min_alleles=2; max_alleles=2; }	// Output as BEAGLE Genotype Likelihood format
		else if (in_str == "--bed") { BED_file = get_arg(i+1); i++; BED_exclude=false; }
		else if (in_str == "--chr") { chrs_to_keep.insert(get_arg(i+1)); i++; }					// Chromosome to process
		else if (in_str == "--counts") output_counts = true;								// Output per-site allele count statistics
		else if (in_str == "--counts2") {output_counts = true; suppress_allele_output = true; }								// Output per-site allele count statistics
		else if (in_str == "--depth") output_indv_depth = true;							// Output per-individual coverage statistics
		else if (in_str == "--derived") derived = true;								// Output frequencies so that AA is first.
		else if (in_str == "--diff-discordance-matrix") { diff_discordance_matrix = true; }	// Calculate some concensus statistics
		else if (in_str == "--diff-indv-discordance") { diff_indv_discordance = true; }	// Calculate some concensus statistics
		else if (in_str == "--diff-indv-map") { diff_indv_map_file = get_arg(i+1); i++; }	// Map individual ids in the second file
		else if (in_str == "--diff-site-discordance") { diff_site_discordance = true; }	// Calculate some concensus statistics
		else if (in_str == "--diff-switch-error") { diff_switch_error = true; }	// Calculate some concensus statistics
		else if (in_str == "--diff") { diff_file = get_arg(i+1); diff_file_compressed = false; i++; }	// Calculate some concensus statistics
		else if (in_str == "--diff-bcf") { diff_file = get_arg(i+1); diff_file_bcf = true; i++; }	// Calculate some concensus statistics
		else if (in_str == "--exclude-bed") { BED_file = get_arg(i+1); i++; BED_exclude=true; }
		else if (in_str == "--exclude") { snps_to_exclude_file = get_arg(i+1); i++; }				// List of SNPs to exclude
		else if (in_str == "--exclude-positions") { exclude_positions_file = get_arg(i+1); i++; }
		else if (in_str == "--extract-FORMAT-info") { FORMAT_id_to_extract = get_arg(i+1); i++; }
		else if (in_str == "--FILTER-summary") {output_filter_summary = true;}
		else if (in_str == "--filtered-sites") {output_filtered_sites = true;}
		else if (in_str == "--force-index-write") {force_write_index = true;}
		else if (in_str == "--freq") output_freq = true;									// Output per-site frequency statistics
		else if (in_str == "--freq2") {output_freq = true;	suppress_allele_output = true; }							// Output per-site frequency statistics
		else if (in_str == "--from-bp") { start_pos = atoi(get_arg(i+1).c_str()); i++; }					// Start position
		else if (in_str == "--fst-window-size") { fst_window_size = atoi(get_arg(i+1).c_str()); i++; }                  // Window size for Fst calculation
		else if (in_str == "--fst-window-step") { fst_window_step = atoi(get_arg(i+1).c_str()); i++; }                  // Window step for Fst calculation
		else if (in_str == "--gatk") gatk = true;								//Denote BCF file came from GATK
		else if (in_str == "--geno-depth") output_geno_depth = true;						// Output Depth for each genoptype
		else if (in_str == "--geno-r2") { output_geno_rsq = true; min_alleles = 2; max_alleles = 2; } // Output pairwise LD (r^2)
		else if (in_str == "--geno-chisq") { output_geno_chisq = true; } // Output pairwise LD (r^2)
		else if (in_str == "--geno-r2-positions") { geno_rsq_position_list = get_arg(i+1); i++; min_alleles = 2; max_alleles = 2; }	// Output pairwise LD (r^2)
		else if (in_str == "--geno") { min_site_call_rate = atof(get_arg(i+1).c_str()); i++; }				// Minimum per-site call rate
		else if (in_str == "--get-INFO") { INFO_to_extract.push_back(get_arg(i+1)); i++; }	// Add to list of INFO fields to extract
		else if (in_str == "--gzdiff") { diff_file = get_arg(i+1); diff_file_compressed = true; i++; }	// Calculate some concensus statistics
		else if (in_str == "--gzvcf") { vcf_filename = get_arg(i+1); vcf_compressed = true; i++; } 				// VCF file to process
		else if (in_str == "--hap-r2") { output_hap_rsq = true; phased_only = true; min_alleles = 2; max_alleles = 2; }	// Output pairwise LD (r^2)
		else if (in_str == "--hap-r2-positions") { hap_rsq_position_list = get_arg(i+1); i++; phased_only = true; min_alleles = 2; max_alleles = 2; }	// Output pairwise LD (r^2)
		else if (in_str == "--hapmap-fst-pop") { hapmap_fst_populations.push_back(get_arg(i+1)); i++; }
		else if (in_str == "--hardy") output_HWE = true;									// Output HWE statistics
		else if (in_str == "--het") output_het = true;									// Output heterozygosity statistics
		else if (in_str == "--hist-indel-len") output_indel_hist = true;
		else if (in_str == "--hwe") { max_alleles = 2; min_HWE_pvalue = atof(get_arg(i+1).c_str()); i++; }					// Minimum per-site HWE p-value
		else if (in_str == "--IMPUTE") { output_as_IMPUTE = true; phased_only=true; min_site_call_rate=1.0; min_alleles=2; max_alleles=2; }// Output as IMPUTE format
		else if (in_str == "--indv") { indv_to_keep.insert(get_arg(i+1)); i++; }						// List of individuals to keep
		else if (in_str == "--interchrom-hap-r2") { output_interchromosomal_hap_rsq = true; phased_only = true;  min_alleles = 2; max_alleles = 2;	}	// Output pairwise LD (r^2)
		else if (in_str == "--interchrom-geno-r2") { output_interchromosomal_geno_rsq = true; min_alleles = 2; max_alleles = 2;	}	// Output pairwise LD (r^2)
		else if (in_str == "--invert-mask") { mask_file = get_arg(i+1); i++; invert_mask = true; }
		else if (in_str == "--keep-filtered") { site_filter_flags_to_keep.insert(get_arg(i+1)); i++; }	// Remove a specific filter flag
		else if (in_str == "--keep") { indv_keep_file = get_arg(i+1); i++; }						// List of individuals to keep
		else if (in_str == "--keep-only-indels") { keep_only_indels = true; }
		else if (in_str == "--keep-INFO") { site_INFO_flags_to_keep.insert(get_arg(i+1)); i++; }	// Filter sites by INFO flags
		else if (in_str == "--keep-INFO-all") { recode_all_INFO=true; }	// Old command (soon to be depreciated)
		else if (in_str == "--ld-window-bp") { ld_bp_window_size = atoi(get_arg(i+1).c_str()); i++; }	// Max bp distance for LD output
		else if (in_str == "--ld-window") { ld_snp_window_size = atoi(get_arg(i+1).c_str()); i++; }		// Max SNP distance for LD output
		else if (in_str == "--ld-window-bp-min") { ld_bp_window_min = atoi(get_arg(i+1).c_str()); i++; }	// Max bp distance for LD output
		else if (in_str == "--ld-window-min") { ld_snp_window_min = atoi(get_arg(i+1).c_str()); i++; }		// Max SNP distance for LD output
		else if (in_str == "--ldhat-geno") { output_as_ldhat_unphased = true; }
		else if (in_str == "--ldhat") { output_as_ldhat_phased = true; phased_only = true; } // Output as LDhat format
		else if (in_str == "--LROH") {output_LROH = true;}
		else if (in_str == "--mac") { min_mac = atoi(get_arg(i+1).c_str()); i++; }								// Minimum Site MAC
		else if (in_str == "--maf") { min_maf = atof(get_arg(i+1).c_str()); i++; }								// Minimum Site MAF
		else if (in_str == "--mask-min") { min_kept_mask_value = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--mask") { mask_file = get_arg(i+1); i++; invert_mask = false; }
		else if (in_str == "--max-alleles") { max_alleles = atoi(get_arg(i+1).c_str()); i++; }				// Maximum number of alleles per-site
		else if (in_str == "--max-indv-meanDP") { max_indv_mean_depth = atof(get_arg(i+1).c_str()); i++; }	// Maximum mean depth for an individual
		else if (in_str == "--max-mac") { max_mac = atoi(get_arg(i+1).c_str()); i++; }						// Maximum site MAC
		else if (in_str == "--max-maf") { max_maf = atof(get_arg(i+1).c_str()); i++; }						// Maximum Site MAF
		else if (in_str == "--max-meanDP") { max_mean_depth = atof(get_arg(i+1).c_str()); i++; }			// Site Maximum mean depth across individuals
		else if (in_str == "--max-missing-count") { max_missing_call_count = atoi(get_arg(i+1).c_str()); i++; } // Site maximum missing genotypes
		else if (in_str == "--max-non-ref-ac") { max_non_ref_ac = atoi(get_arg(i+1).c_str()); i++; }		// Minimum Site non-ref AC
		else if (in_str == "--max-non-ref-af") { max_non_ref_af = atof(get_arg(i+1).c_str()); i++; }		// Minimum Site non-ref AF
		else if (in_str == "--maxDP") { max_genotype_depth = atoi(get_arg(i+1).c_str()); i++; }				// Maximum genotype depth
		else if (in_str == "--max-indv") {max_N_indv = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--min-alleles") { min_alleles = atoi(get_arg(i+1).c_str()); i++; }				// Minimum number of alleles per-site
		else if (in_str == "--min-indv-meanDP") { min_indv_mean_depth = atof(get_arg(i+1).c_str()); i++; }	// Minimum mean depth for an individual
		else if (in_str == "--min-meanDP") { min_mean_depth = atof(get_arg(i+1).c_str()); i++; }			// Site Minimum mean depth
		else if (in_str == "--min-r2") { min_r2 = atof(get_arg(i+1).c_str()); i++; }					// Min r^2 for LD output
		else if (in_str == "--mind") { min_indv_call_rate = atof(get_arg(i+1).c_str()); i++; }				// Minimum per-individual call rate
		else if (in_str == "--minDP") { min_genotype_depth = atoi(get_arg(i+1).c_str()); i++; }				// Minimum genotype depth
		else if (in_str == "--minGQ") { min_genotype_quality = atof(get_arg(i+1).c_str()); i++; }			// Minimum genotype quality
		else if (in_str == "--minQ") { min_quality = atof(get_arg(i+1).c_str()); i++; }						// Minimum per-site quality
		else if (in_str == "--missing") output_missingness = true;						// Output Individual and Site missingness summaries
		else if (in_str == "--non-ref-ac") { min_non_ref_ac = atoi(get_arg(i+1).c_str()); i++; }				// Minimum Site non-ref AC
		else if (in_str == "--non-ref-af") { min_non_ref_af = atof(get_arg(i+1).c_str()); i++; }				// Minimum Site non-ref AF
		else if (in_str == "--not-chr") { chrs_to_exclude.insert(get_arg(i+1)); i++; }					// Chromosome to process
		else if (in_str == "--out") { output_prefix = get_arg(i+1); i++; }							// Output file prefix
		else if (in_str == "--pca") { output_PCA = true; min_alleles=2; max_alleles=2; }
		else if (in_str == "--pca-no-norm") { output_PCA = true; PCA_no_normalisation = true; min_alleles=2; max_alleles=2; }
		else if (in_str == "--pca-snp-loadings") { output_N_PCA_SNP_loadings = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--phased") phased_only = true;								// Keep only phased individuals / sites
		else if (in_str == "--plink") plink_output = true;								// Output as PLINK file
		else if (in_str == "--plink-tped") plink_tped_output = true;	// Output as PLINK tped file
		else if (in_str == "--positions") { positions_file = get_arg(i+1); i++; }
		else if (in_str == "--recode") recode = true;										// Output VCF file
		else if (in_str == "--recode-bcf") recode_bcf = true;										// Output VCF file
		else if (in_str == "--recode-bcf-to-stream") {  recode_bcf_to_stream = true; LOG.set_screen_output(false); }
		else if (in_str == "--recode-to-stream") {  recode_to_stream = true; LOG.set_screen_output(false); }
		else if (in_str == "--recode-INFO-all") { recode_all_INFO=true; }		// Specify INFO to keep when recoding
		else if (in_str == "--recode-INFO") { recode_INFO_to_keep.insert(get_arg(i+1)); i++; }		// Specify INFO to keep when recoding
		else if (in_str == "--relatedness") { output_relatedness = true; }	// Estimate relatedness between individuals
		else if (in_str == "--remove-filtered-all") remove_all_filtered_sites = true;							// Remove sites flagged as filtered
		else if (in_str == "--remove-filtered-geno-all") remove_all_filtered_genotypes = true;			// Remove genotypes flagged as filtered
		else if (in_str == "--remove-filtered-geno") { geno_filter_flags_to_exclude.insert(get_arg(i+1)); i++; }		// Remove genotypes flagged as filtered
		else if (in_str == "--remove-filtered") { site_filter_flags_to_exclude.insert(get_arg(i+1)); i++; }	// Remove a specific filter flag
		else if (in_str == "--remove-indels") { remove_indels = true; }
		else if (in_str == "--remove-indv") { indv_to_exclude.insert(get_arg(i+1)); i++; }						// List of individuals to keep
		else if (in_str == "--remove-INFO") { site_INFO_flags_to_remove.insert(get_arg(i+1)); i++; }	// Filter sites by INFO flags
		else if (in_str == "--remove") { indv_exclude_file = get_arg(i+1); i++; }					// List of individuals to exclude
		else if (in_str == "--singletons") output_singletons = true;							// Output as 0/1/2 mat
		else if (in_str == "--site-depth") output_site_depth = true;						// Output Depth for each site
		else if (in_str == "--site-mean-depth") output_site_mean_depth = true;			// Output Mean Depth for each site
		else if (in_str == "--site-pi") { output_site_pi = true; }
		else if (in_str == "--site-quality") output_site_quality = true;								// Output per-site qualities
		else if (in_str == "--snp") { snps_to_keep.insert(get_arg(i+1)); i++; }						// SNP to keep
		else if (in_str == "--SNPdensity") { output_SNP_density_bin_size = atoi(get_arg(i+1).c_str()); i++; }	// Output SNP density using Bin Size
		else if (in_str == "--snps") { snps_to_keep_file = get_arg(i+1); i++; }						// List of SNPs to keep
		else if (in_str == "--TajimaD") { output_Tajima_D_bin_size = atoi(get_arg(i+1).c_str()); i++;}
		else if (in_str == "--to-bp") { end_pos = atoi(get_arg(i+1).c_str()); i++; }						// End position
		else if (in_str == "--thin") { min_interSNP_distance = atoi(get_arg(i+1).c_str()); i++; }	// Set minimum distance between SNPs
		else if (in_str == "--TsTv") {output_TsTv_bin_size = atoi(get_arg(i+1).c_str()); i++;}						// Output Ts/Tv stats
		else if (in_str == "--TsTv-by-count") {output_TsTv_by_count = true; }						// Output Ts/Tv stats
		else if (in_str == "--TsTv-by-qual") {output_TsTv_by_qual = true; }
		else if (in_str == "--weir-fst-pop") { weir_fst_populations.push_back(get_arg(i+1)); i++; }
		else if (in_str == "--window-pi") { pi_window_size = atoi(get_arg(i+1).c_str()); i++; }
		else if (in_str == "--window-pi-step") { pi_window_step = atoi(get_arg(i+1).c_str()); i++; }
		else
			error("Unknown option: " + string(in_str), 0);
		i++;
	}
	check_parameters();
}

string parameters::get_arg(unsigned int i)
{
	if (i>=argv.size())
		error("Requested Missing Argument",76);
	return argv[i];
}

void parameters::print_params()
{
	parameters defaults(0, 0);

	LOG.printLOG("Parameters as interpreted:\n");
	if (vcf_filename != defaults.vcf_filename)
	{
		if (bcf_format == true)
			LOG.printLOG("\t--bcf " + vcf_filename + "\n");
		else if (vcf_compressed == false)
			LOG.printLOG("\t--vcf " + vcf_filename + "\n");
		else
			LOG.printLOG("\t--gzvcf " + vcf_filename + "\n");
	}

	if (chrs_to_keep.size() > 0)
	{
		for (set<string>::iterator it=chrs_to_keep.begin(); it != chrs_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--chr " + tmp + "\n");
		}
	}

	if (chrs_to_exclude.size() > 0)
	{
		for (set<string>::iterator it=chrs_to_exclude.begin(); it != chrs_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--not-chr " + tmp + "\n");
		}
	}

	if (derived != defaults.derived) LOG.printLOG("\t--derived\n");
	if (end_pos != defaults.end_pos) LOG.printLOG("\t--to-bp " + output_log::int2str(end_pos) + "\n");
	if (exclude_positions_file != defaults.exclude_positions_file) LOG.printLOG("\t--exclude-positions " + exclude_positions_file + "\n");
	if (force_write_index != defaults.force_write_index) LOG.printLOG("\t--force-index-write\n");
	if (FORMAT_id_to_extract != defaults.FORMAT_id_to_extract) LOG.printLOG("\t--extract-FORMAT-info " + FORMAT_id_to_extract + "\n");
	if (gatk != defaults.gatk) LOG.printLOG("\t--gatk\n");
	if (geno_rsq_position_list != defaults.geno_rsq_position_list) LOG.printLOG("\t--geno-r2-positions " + geno_rsq_position_list + "\n");
	if (hap_rsq_position_list != defaults.hap_rsq_position_list) LOG.printLOG("\t--hap-r2-positions " + hap_rsq_position_list + "\n");
	if (hapmap_fst_populations.size() != 0)
	{
		for (unsigned int ui=0; ui<hapmap_fst_populations.size(); ui++)
			LOG.printLOG("\t--hapmap-fst-pop " + hapmap_fst_populations[ui] + "\n");
	}
	if (weir_fst_populations.size() != 0)
	{
		for (unsigned int ui=0; ui<weir_fst_populations.size(); ui++)
			LOG.printLOG("\t--weir-fst-pop " + weir_fst_populations[ui] + "\n");
	}
	if (indv_exclude_file != defaults.indv_exclude_file) LOG.printLOG("\t--exclude " + indv_exclude_file + "\n");
	if (indv_keep_file != defaults.indv_keep_file) LOG.printLOG("\t--keep " + indv_keep_file + "\n");
	if (recode_all_INFO == true) LOG.printLOG("\t--recode-INFO-all\n");
	if (keep_only_indels != defaults.keep_only_indels) LOG.printLOG("\t--keep-only-indels\n");
	if (ld_bp_window_size != defaults.ld_bp_window_size) LOG.printLOG("\t--ld-window-bp " + output_log::int2str(ld_bp_window_size) + "\n");
	if (ld_snp_window_size != defaults.ld_snp_window_size) LOG.printLOG("\t--ld-window " + output_log::int2str(ld_snp_window_size) + "\n");
	if (ld_bp_window_min != defaults.ld_bp_window_min) LOG.printLOG("\t--ld-window-bp-min " + output_log::int2str(ld_bp_window_min) + "\n");
	if (ld_snp_window_min != defaults.ld_snp_window_min) LOG.printLOG("\t--ld-window-min " + output_log::int2str(ld_snp_window_min) + "\n");
	if (min_mac != defaults.min_mac) LOG.printLOG("\t--mac " + output_log::dbl2str(min_mac, 3) + "\n");
	if (min_maf != defaults.min_maf) LOG.printLOG("\t--maf " + output_log::dbl2str(min_maf, 3) + "\n");
	if (max_alleles != defaults.max_alleles) LOG.printLOG("\t--max-alleles " + output_log::int2str(max_alleles) + "\n");
	if (max_genotype_depth != defaults.max_genotype_depth) LOG.printLOG("\t--maxDP " + output_log::dbl2str(max_genotype_depth, 3) + "\n");
	if (max_indv_mean_depth != defaults.max_indv_mean_depth) LOG.printLOG("\t--max-indv-meanDP " + output_log::dbl2str(max_indv_mean_depth, 3) + "\n");
	if (max_mac != defaults.max_mac) LOG.printLOG("\t--max-mac " + output_log::dbl2str(max_mac, 3) + "\n");
	if (max_maf != defaults.max_maf) LOG.printLOG("\t--max-maf " + output_log::dbl2str(max_maf, 3) + "\n");
	if (max_missing_call_count != defaults.max_missing_call_count) LOG.printLOG("\t--max-missing-count " + output_log::dbl2str(max_missing_call_count, 3) + "\n");
	if (max_mean_depth != defaults.max_mean_depth) LOG.printLOG("\t--max-meanDP " + output_log::dbl2str(max_mean_depth, 3) + "\n");
	if (max_non_ref_ac != defaults.max_non_ref_ac) LOG.printLOG("\t--max-non-ref-ac " + output_log::dbl2str(max_non_ref_ac, 3) + "\n");
	if (max_non_ref_af != defaults.max_non_ref_af) LOG.printLOG("\t--max-non-ref-af " + output_log::dbl2str(max_non_ref_af, 3) + "\n");
	if (max_N_indv != defaults.max_N_indv) LOG.printLOG("\t--max-indv " + output_log::int2str(max_N_indv) + "\n");
	if (min_alleles != defaults.min_alleles) LOG.printLOG("\t--min-alleles " + output_log::int2str(min_alleles) + "\n");
	if (min_genotype_depth != defaults.min_genotype_depth) LOG.printLOG("\t--minDP " + output_log::dbl2str(min_genotype_depth, 3) + "\n");
	if (min_genotype_quality != defaults.min_genotype_quality) LOG.printLOG("\t--minGQ " + output_log::dbl2str(min_genotype_quality, 3) + "\n");
	if (min_HWE_pvalue != defaults.min_HWE_pvalue) LOG.printLOG("\t--hwe " + output_log::dbl2str(min_HWE_pvalue, 3) + "\n");
	if (min_indv_call_rate != defaults.min_indv_call_rate) LOG.printLOG("\t--mind " + output_log::dbl2str(min_indv_call_rate, 3) + "\n");
	if (min_indv_mean_depth != defaults.min_indv_mean_depth) LOG.printLOG("\t--min-indv-meanDP " + output_log::dbl2str(min_indv_mean_depth, 3) + "\n");
	if (min_interSNP_distance != defaults.min_interSNP_distance) LOG.printLOG("\t--thin " + output_log::int2str(min_interSNP_distance) + "\n");
	if (min_kept_mask_value != defaults.min_kept_mask_value) LOG.printLOG("\t--mask-min " + output_log::int2str(min_kept_mask_value) + "\n");
	if (min_mean_depth != defaults.min_mean_depth) LOG.printLOG("\t--min-meanDP " + output_log::dbl2str(min_mean_depth, 3) + "\n");
	if (min_quality != defaults.min_quality) LOG.printLOG("\t--minQ " + output_log::dbl2str(min_quality, 3) + "\n");
	if (min_r2 != defaults.min_r2) LOG.printLOG("\t--min-r2 " + output_log::dbl2str(min_r2, 3) + "\n");
	if (min_site_call_rate != defaults.min_site_call_rate) LOG.printLOG("\t--geno " + output_log::dbl2str(min_site_call_rate, 3) + "\n");
	if (min_non_ref_ac != defaults.min_non_ref_ac) LOG.printLOG("\t--non-ref-ac " + output_log::dbl2str(min_non_ref_ac, 3) + "\n");
	if (min_non_ref_af != defaults.min_non_ref_af) LOG.printLOG("\t--non-ref-af " + output_log::dbl2str(min_non_ref_af, 3) + "\n");
	if (output_012_matrix) LOG.printLOG("\t--012\n");
	if (output_as_IMPUTE) LOG.printLOG("\t--IMPUTE\n");
	if (output_BEAGLE_genotype_likelihoods_GL) LOG.printLOG("\t--BEAGLE-GL\n");
	if (output_BEAGLE_genotype_likelihoods_PL) LOG.printLOG("\t--BEAGLE-PL\n");
	if (output_counts && (suppress_allele_output==false)) LOG.printLOG("\t--counts\n");
	if (output_counts && (suppress_allele_output==true)) LOG.printLOG("\t--counts2\n");
	if (output_filter_summary) LOG.printLOG("\t--FILTER-summary\n");
	if (output_filtered_sites != defaults.output_filtered_sites) LOG.printLOG("\t--filtered-sites\n");
	if (output_freq && (suppress_allele_output==false)) LOG.printLOG("\t--freq\n");
	if (output_freq && (suppress_allele_output==true)) LOG.printLOG("\t--freq2\n");
	if (output_geno_depth) LOG.printLOG("\t--geno-depth\n");
	if (output_geno_rsq) LOG.printLOG("\t--geno-r2\n");
	if (output_geno_chisq) LOG.printLOG("\t--geno-chisq\n");
	if (output_hap_rsq) LOG.printLOG("\t--hap-r2\n");
	if (output_het) LOG.printLOG("\t--het\n");
	if (output_HWE) LOG.printLOG("\t--hardy\n");
	if (output_indel_hist != defaults.output_indel_hist) LOG.printLOG("\t--hist-indel-len \n");
	if (output_indv_depth) LOG.printLOG("\t--depth\n");
	if (output_interchromosomal_hap_rsq) LOG.printLOG("\t--interchrom-hap-r2\n");
	if (output_interchromosomal_geno_rsq) LOG.printLOG("\t--interchrom-geno-r2\n");
	if (output_LROH != defaults.output_LROH) LOG.printLOG("\t--LROH\n");
	if (output_missingness) LOG.printLOG("\t--missing\n");
	if (output_PCA != defaults.output_PCA)
	{
		if (PCA_no_normalisation == false)
			LOG.printLOG("\t--pca\n");
		else
			LOG.printLOG("\t--pca-no-norm\n");

		if (output_N_PCA_SNP_loadings != defaults.output_N_PCA_SNP_loadings)
			LOG.printLOG("\t--pca-snp-loadings " + output_log::int2str(output_N_PCA_SNP_loadings) + "\n");
	}
	if (output_prefix != defaults.output_prefix) LOG.printLOG("\t--out " + output_prefix + "\n");
	if (output_relatedness) LOG.printLOG("\t--relatedness\n");
	if (output_singletons) LOG.printLOG("\t--singletons\n");
	if (output_site_depth) LOG.printLOG("\t--site-depth\n");
	if (output_site_mean_depth) LOG.printLOG("\t--site-mean-depth\n");
	if (output_site_pi != defaults.output_site_pi) LOG.printLOG("\t--site-pi\n");
	if (output_site_quality) LOG.printLOG("\t--site-quality\n");
	if (output_SNP_density_bin_size != defaults.output_SNP_density_bin_size) LOG.printLOG("\t--SNPdensity " + output_log::int2str(output_SNP_density_bin_size) + "\n");
	if (output_TsTv_bin_size != defaults.output_TsTv_bin_size) LOG.printLOG("\t--TsTv " + output_log::int2str(output_TsTv_bin_size) + "\n");
	if (output_TsTv_by_count) LOG.printLOG("\t--TsTv-by-count\n");
	if (output_TsTv_by_qual) LOG.printLOG("\t--TsTv-by-qual\n");
	if (phased_only) LOG.printLOG("\t--phased\n");
	if (pi_window_size != defaults.pi_window_size) LOG.printLOG("\t--window-pi " + output_log::int2str(pi_window_size) + "\n");
	if (pi_window_step != defaults.pi_window_step) LOG.printLOG("\t--window-pi-step " + output_log::int2str(pi_window_step) + "\n");
	if (plink_output) LOG.printLOG("\t--plink\n");
	if (plink_tped_output) LOG.printLOG("\t--plink-tped\n");
	if (positions_file != defaults.positions_file) LOG.printLOG("\t--positions " + positions_file + "\n");
	if (recode) LOG.printLOG("\t--recode\n");
	if (recode_bcf) LOG.printLOG("\t--recode-bcf\n");
	if (recode_to_stream) LOG.printLOG("\t--recode-to-stream\n");
	if (recode_bcf_to_stream) LOG.printLOG("\t--recode-bcf-to-stream\n");
	if (remove_all_filtered_genotypes) LOG.printLOG("\t--remove-filtered-geno-all\n");
	if (remove_all_filtered_sites) LOG.printLOG("\t--remove-filtered-all\n");
	if (remove_indels != defaults.remove_indels) LOG.printLOG("\t--remove-indels\n");
	if (snps_to_exclude_file != defaults.snps_to_exclude_file) LOG.printLOG("\t--exclude " + snps_to_exclude_file + "\n");
	if (snps_to_keep_file != defaults.snps_to_keep_file) LOG.printLOG("\t--snps " + snps_to_keep_file + "\n");
	if (start_pos != defaults.start_pos) LOG.printLOG("\t--from-bp " + output_log::int2str(start_pos) + "\n");
	if (output_Tajima_D_bin_size != defaults.output_Tajima_D_bin_size) LOG.printLOG("\t--TajimaD " + output_log::int2str(output_Tajima_D_bin_size) + "\n");

	if (output_as_ldhat_phased) LOG.printLOG("\t--ldhat\n");
	if (output_as_ldhat_unphased) LOG.printLOG("\t--ldhat-geno\n");

	if (site_filter_flags_to_exclude.size() > 0)
		for (set<string>::iterator it=site_filter_flags_to_exclude.begin(); it != site_filter_flags_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-filtered " + tmp + "\n");
		}

	if (site_filter_flags_to_keep.size() > 0)
		for (set<string>::iterator it=site_filter_flags_to_keep.begin(); it != site_filter_flags_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--keep-filtered " + tmp + "\n");
		}

	if (geno_filter_flags_to_exclude.size() > 0)
		for (set<string>::iterator it=geno_filter_flags_to_exclude.begin(); it != geno_filter_flags_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-filtered-geno " + tmp + "\n");
		}

	if (INFO_to_extract.size() > 0)
		for (unsigned int ui=0; ui<INFO_to_extract.size(); ui++)
			LOG.printLOG("\t--get-INFO " + INFO_to_extract[ui] + "\n");

	if (diff_file != defaults.diff_file)
	{
		if ( diff_file_bcf == true )
			LOG.printLOG("\t--diff-bcf " + diff_file + "\n");
		else if (diff_file_compressed == false)
			LOG.printLOG("\t--diff " + diff_file + "\n");
		else
			LOG.printLOG("\t--gzdiff " + diff_file + "\n");
		if (diff_site_discordance == true) LOG.printLOG("\t--diff-site-discordance\n");
		if (diff_indv_discordance == true) LOG.printLOG("\t--diff-indv-discordance\n");
		if (diff_discordance_matrix == true) LOG.printLOG("\t--diff-discordance-matrix\n");
		if (diff_switch_error == true) LOG.printLOG("\t--diff-switch-error\n");
		if (diff_indv_map_file != defaults.diff_indv_map_file) LOG.printLOG("\t--diff-indv-map " + diff_indv_map_file + "\n");
	}

	if (recode_INFO_to_keep.size() > 0)
		for (set<string>::iterator it=recode_INFO_to_keep.begin(); it != recode_INFO_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--recode-INFO " + tmp + "\n");
		}

	if (site_INFO_flags_to_remove.size() > 0)
		for (set<string>::iterator it=site_INFO_flags_to_remove.begin(); it != site_INFO_flags_to_remove.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-INFO " + tmp + "\n");
		}

	if (site_INFO_flags_to_keep.size() > 0)
		for (set<string>::iterator it=site_INFO_flags_to_keep.begin(); it != site_INFO_flags_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--keep-INFO " + tmp + "\n*** Note: --keep-INFO has changed. Are you sure you don't want --recode-INFO? ***\n");
		}

	if (BED_file != defaults.BED_file)
	{
		if (BED_exclude == false)
			LOG.printLOG("\t--bed " + BED_file + "\n");
		else
			LOG.printLOG("\t--exclude-bed " + BED_file + "\n");
	}

	if (mask_file != defaults.mask_file)
	{
		if (invert_mask == false)
			LOG.printLOG("\t--mask " + mask_file + "\n");
		else
			LOG.printLOG("\t--invert-mask " + mask_file + "\n");
	}

	if (snps_to_keep.size() > 0)
		for (set<string>::iterator it=snps_to_keep.begin(); it != snps_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--snp " + tmp + "\n");
		}

	if (indv_to_keep.size() > 0)
		for (set<string>::iterator it=indv_to_keep.begin(); it != indv_to_keep.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--indv " + tmp + "\n");
		}

	if (indv_to_exclude.size() > 0)
		for (set<string>::iterator it=indv_to_exclude.begin(); it != indv_to_exclude.end(); ++it)
		{
			string tmp = *it;
			LOG.printLOG("\t--remove-indv " + tmp + "\n");
		}

	LOG.printLOG("\n");
}

void parameters::print_help()
{
	unsigned int i;
	string in_str;

	if (argv.size() <= 1)
	{	// If there are no user parameters, display help.
		argv.push_back("--?");
		print_help();
	}

	for(i = 0; i < argv.size(); i++)
	{
		in_str = argv[i];
		if ((in_str == "-h") || (in_str == "-?") || (in_str == "-help") || (in_str == "--?") || (in_str == "--help") || (in_str == "--h"))
		{
			cout << endl << "VCFtools (" << VCFTOOLS_VERSION << ")" << endl;
			cout << "\u00A9 Adam Auton 2009" << endl << endl;
			cout << "Process Variant Call Format files" << endl;
			cout << endl;
			cout << "For a list of options, please go to:" << endl;
			cout << "\thttp://vcftools.sourceforge.net/options.html" << endl;
			cout << endl;
			cout << "Questions, comments, and suggestions should be emailed to:" << endl;
			cout << "\tvcftools-help@lists.sourceforge.net" << endl;
			cout << endl;

			exit(0);
		}
	}
}

void parameters::check_parameters()
{
	parameters defaults(0, 0);
	if (vcf_filename == "") error("VCF required.", 0);
	if (end_pos < start_pos) error("End position must be greater than Start position.", 1);
	if (((end_pos != numeric_limits<int>::max()) || (start_pos != -1)) && (chrs_to_keep.size() != 1)) error("Require a single chromosome when specifying a range.", 2);
	if (max_maf < min_maf) error("Maximum MAF must be not be less than Minimum MAF.", 4);
	if (max_mac < min_mac) error("Maximum MAC must be not be less than Minimum MAC.", 4);
	if (min_maf != defaults.min_maf)
	{
		if ((min_maf < 0.0) || (min_maf > 1.0)) error("MAF must be between 0 and 1.", 4);
	}
	if (max_maf != defaults.max_maf)
	{
		if ((max_maf < 0.0) || (max_maf > 1.0)) error("Maximum MAF must be between 0 and 1.", 4);
	}
	if (min_non_ref_af != defaults.min_non_ref_af)
	{
		if ((min_non_ref_af < 0.0) || (min_non_ref_af > 1.0)) error("Non-Ref Allele Frequency must be between 0 and 1.", 4);
	}
	if (max_non_ref_af < min_non_ref_af) error("Maximum Non-Ref Allele Frequency must not be less that Minimum Non-Ref AF.", 4);
	if (max_non_ref_ac < min_non_ref_ac) error("Maximum Non-Ref Allele Count must not be less that Minimum Non-Ref AC.", 4);
	if ((min_site_call_rate > 1) || (min_indv_call_rate > 1)) error("Minimum Call rates cannot be greater than 1.", 5);
	if (max_alleles < min_alleles) error("Max Number of Alleles must be greater than Min Number of Alleles.", 6);
	if (max_mean_depth < min_mean_depth) error("Max Mean Depth must be greater the Min Mean Depth.", 7);
	if (max_indv_mean_depth < min_indv_mean_depth) error("Max Indv Mean Depth must be greater the Min Indv Mean Depth.", 8);
	if (max_genotype_depth < min_genotype_depth) error("Max Genotype Depth must be greater than Min Genotype Depth.", 9);
	if (((output_as_ldhat_phased == true) || (output_as_ldhat_unphased)) && (chrs_to_keep.size() != 1)) error("Require a chromosome (--chr) when outputting LDhat format.", 11);
	if ((output_BEAGLE_genotype_likelihoods_GL == true) && (chrs_to_keep.size() != 1)) error("Require a chromosome (--chr) when outputting Beagle likelihoods.", 11);
	if ((output_BEAGLE_genotype_likelihoods_PL == true) && (chrs_to_keep.size() != 1)) error("Require a chromosome (--chr) when outputting Beagle likelihoods.", 11);
	if (min_kept_mask_value > 9) error("Min Mask value must be between 0 and 9.", 14);
	if ((output_LROH == true) && (chrs_to_keep.size() != 1)) error("Require a chromosome (--chr) when outputting LROH.", 11);
	if (output_TsTv_bin_size < 0) error("TsTv bin size must be > 0",16);
	if (output_Tajima_D_bin_size < 0) error("Tajima D bin size must be > 0", 17);
	if (pi_window_size < 0) error("Pi Window size must be > 0", 18);
	if (output_SNP_density_bin_size < 0) error("SNP density bin size must be > 0", 18);
}

void parameters::error(string err_msg, int code)
{
	LOG.printLOG("\n\nError: " + err_msg + "\n\n");
	exit(code);
}
