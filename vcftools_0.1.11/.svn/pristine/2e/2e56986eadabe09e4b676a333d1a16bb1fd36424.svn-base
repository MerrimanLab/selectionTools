/*
 * vcftools.cpp
 */

#include "vcftools.h"

output_log LOG;

int main(int argc, char *argv[])
{
	time_t start,end;
	time(&start);

	// The following turns off sync between C and C++ streams.
	// Apparently it's faster to turn sync off, and as I don't use C streams, it's okay to turn off.
	ios_base::sync_with_stdio(false);

	parameters params(argc, argv);
	params.print_help();
	params.read_parameters();

	LOG.open(params.output_prefix);

	LOG.printLOG("\nVCFtools - " + VCFTOOLS_VERSION + "\n");
	LOG.printLOG("(C) Adam Auton 2009\n\n");

	params.print_params();

	variant_file *vf;
	if (!params.bcf_format)
		vf = new vcf_file(params.vcf_filename, params.vcf_compressed, params.chrs_to_keep, params.chrs_to_exclude, params.force_write_index);
	else
		vf = new bcf_file(params.vcf_filename, params.chrs_to_keep, params.chrs_to_exclude, params.force_write_index, params.gatk);

	vf->apply_filters(params);

	unsigned int N_indv = vf->N_kept_individuals();
	unsigned int N_sites = vf->N_kept_sites();
	unsigned int N_total_indv = vf->N_total_indv();
	unsigned int N_total_sites = vf->N_total_sites();
	LOG.printLOG("After filtering, kept " + output_log::int2str(N_indv) + " out of " + output_log::int2str(N_total_indv) + " Individuals\n");
	LOG.printLOG("After filtering, kept " + output_log::int2str(N_sites) + " out of a possible " + output_log::int2str(N_total_sites) + " Sites\n");
	if (N_sites == 0)
		LOG.error("No data left for analysis!");

	if (params.diff_file != "")
	{	// Merge files - cannot be run with other output options.
		variant_file *variant_diff;
		if (params.diff_file_bcf)
			variant_diff = new bcf_file(params.diff_file, params.chrs_to_keep, params.chrs_to_exclude, params.force_write_index, params.gatk);
		else
			variant_diff = new vcf_file(params.diff_file, params.diff_file_compressed, params.chrs_to_keep, params.chrs_to_exclude, params.force_write_index);

		variant_diff->apply_filters(params);	// Apply various filters as required.
		vf->output_indv_in_files(params.output_prefix, *variant_diff, params.diff_indv_map_file);
		vf->output_sites_in_files(params.output_prefix, *variant_diff);

		if (params.diff_site_discordance == true) vf->output_discordance_by_site(params.output_prefix, *variant_diff, params.diff_indv_map_file);
		if (params.diff_discordance_matrix == true) vf->output_discordance_matrix(params.output_prefix, *variant_diff, params.diff_indv_map_file);
		if (params.diff_indv_discordance == true) vf->output_discordance_by_indv(params.output_prefix, *variant_diff, params.diff_indv_map_file);
		if (params.diff_switch_error == true) vf->output_switch_error(params.output_prefix, *variant_diff, params.diff_indv_map_file);
		delete variant_diff;
	}

	vf->output_INFO_for_each_site(params.output_prefix, params.INFO_to_extract);
	vf->output_FORMAT_information(params.output_prefix, params.FORMAT_id_to_extract);
	if (params.output_indv_depth == true) vf->output_individuals_by_mean_depth(params.output_prefix);
	if (params.output_geno_depth == true) vf->output_genotype_depth(params.output_prefix);
	if (params.output_site_depth == true) vf->output_site_depth(params.output_prefix, false);
	if (params.output_site_mean_depth == true) vf->output_site_depth(params.output_prefix, true);
	if (params.output_freq == true) vf->output_frequency(params.output_prefix, false, params.suppress_allele_output, params.derived);
	if (params.output_counts == true) vf->output_frequency(params.output_prefix, true, params.suppress_allele_output, params.derived);
	if (params.plink_output == true) vf->output_as_plink(params.output_prefix);
	if (params.plink_tped_output == true) vf->output_as_plink_tped(params.output_prefix);
	if (params.output_HWE == true) vf->output_hwe(params.output_prefix);
	if (params.output_SNP_density_bin_size > 0) vf->output_SNP_density(params.output_prefix, params.output_SNP_density_bin_size);
	if (params.output_missingness == true) vf->output_missingness(params.output_prefix);
	if (params.output_geno_chisq == true) vf->output_genotype_chisq(params.output_prefix, params.ld_snp_window_size, params.ld_snp_window_min, params.ld_bp_window_size, params.ld_bp_window_min, -1.0);
	if (params.output_geno_rsq == true) vf->output_genotype_r2(params.output_prefix, params.ld_snp_window_size, params.ld_snp_window_min, params.ld_bp_window_size, params.ld_bp_window_min, params.min_r2);
	if (params.output_interchromosomal_hap_rsq == true) vf->output_interchromosomal_haplotype_r2(params.output_prefix, params.min_r2);
	if (params.output_interchromosomal_geno_rsq == true) vf->output_interchromosomal_genotype_r2(params.output_prefix, params.min_r2);
	if (params.output_hap_rsq == true) vf->output_haplotype_r2(params.output_prefix, params.ld_snp_window_size, params.ld_snp_window_min, params.ld_bp_window_size, params.ld_bp_window_min, params.min_r2);
	if (params.hap_rsq_position_list != "") vf->output_haplotype_r2_of_SNP_list_vs_all_others(params.output_prefix, params.hap_rsq_position_list, params.min_r2);
	if (params.geno_rsq_position_list != "") vf->output_genotype_r2_of_SNP_list_vs_all_others(params.output_prefix, params.geno_rsq_position_list, params.min_r2);
	if (params.output_het == true) vf->output_het(params.output_prefix);
	if (params.output_site_quality == true) vf->output_site_quality(params.output_prefix);
	if (params.output_012_matrix == true) vf->output_as_012_matrix(params.output_prefix);
	if (params.output_as_IMPUTE == true) vf->output_as_IMPUTE(params.output_prefix);
	if (params.output_BEAGLE_genotype_likelihoods_GL == true) vf->output_BEAGLE_genotype_likelihoods(params.output_prefix, 0);
	if (params.output_BEAGLE_genotype_likelihoods_PL == true) vf->output_BEAGLE_genotype_likelihoods(params.output_prefix, 1);
	if (params.output_as_ldhat_unphased == true) vf->output_as_LDhat_unphased(params.output_prefix);
	if (params.output_as_ldhat_phased == true) vf->output_as_LDhat_phased(params.output_prefix);
	if (params.output_singletons == true) vf->output_singletons(params.output_prefix);
	if (params.output_site_pi == true) vf->output_per_site_nucleotide_diversity(params.output_prefix);
	if (params.pi_window_size > 0) vf->output_windowed_nucleotide_diversity(params.output_prefix, params.pi_window_size, params.pi_window_step);
	if (params.output_Tajima_D_bin_size > 0) vf->output_Tajima_D(params.output_prefix, params.output_Tajima_D_bin_size);
	if (params.output_TsTv_bin_size > 0) vf->output_TsTv(params.output_prefix, params.output_TsTv_bin_size);
	if (params.output_TsTv_by_count) vf->output_TsTv_by_count(params.output_prefix);
	if (params.output_TsTv_by_qual) vf->output_TsTv_by_quality(params.output_prefix);
	if (params.recode == true) vf->print(params.output_prefix, params.recode_INFO_to_keep, params.recode_all_INFO);
	if (params.recode_bcf == true) vf->print_bcf(params.output_prefix, params.recode_INFO_to_keep, params.recode_all_INFO);
	if (params.recode_to_stream == true) vf->print(std::cout, params.recode_INFO_to_keep, params.recode_all_INFO);
	if (params.recode_bcf_to_stream == true) vf->print_bcf("", params.recode_INFO_to_keep, params.recode_all_INFO, true);
	if (params.output_filter_summary == true) vf->output_FILTER_summary(params.output_prefix);
	if (params.output_filtered_sites == true) vf->output_kept_and_removed_sites(params.output_prefix);
	if (params.output_LROH == true) vf->output_LROH(params.output_prefix);
	if (params.output_relatedness == true) vf->output_indv_relatedness(params.output_prefix);
	if (params.output_PCA == true) vf->output_PCA(params.output_prefix, !params.PCA_no_normalisation, params.output_N_PCA_SNP_loadings);
	if (params.fst_window_size <= 0)
	{
		if (params.hapmap_fst_populations.size() > 0) vf->output_hapmap_fst(params.output_prefix, params.hapmap_fst_populations);
		if (params.weir_fst_populations.size() > 0) vf->output_weir_and_cockerham_fst(params.output_prefix, params.weir_fst_populations);
	}
	else
	{
	      if (params.hapmap_fst_populations.size() > 0) vf->output_windowed_hapmap_fst(params.output_prefix, params.hapmap_fst_populations, params.fst_window_size, params.fst_window_step);
	      if (params.weir_fst_populations.size() > 0) vf->output_windowed_weir_and_cockerham_fst(params.output_prefix, params.weir_fst_populations, params.fst_window_size, params.fst_window_step);
	}
	if (params.output_indel_hist == true) vf->output_indel_hist(params.output_prefix);

	time(&end);
	double running_time = difftime(end,start);
	LOG.printLOG("Run Time = " + output_log::dbl2str_fixed(running_time, 2) + " seconds\n");
	LOG.close();
	delete vf;
	return 0;
}



