/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

void variant_file::apply_filters(const parameters &params)
{
	LOG.printLOG("Applying Required Filters.\n");
	// Apply all filters in turn.
	filter_individuals(params.indv_to_keep, params.indv_to_exclude, params.indv_keep_file, params.indv_exclude_file);
	filter_sites_by_allele_type(params.keep_only_indels, params.remove_indels);
	filter_sites(params.snps_to_keep, params.snps_to_keep_file, params.snps_to_exclude_file);
	filter_sites_by_filter_status(params.site_filter_flags_to_exclude, params.site_filter_flags_to_keep, params.remove_all_filtered_sites);
	string chr_to_keep = "";
	if (params.chrs_to_keep.size() == 1)
		chr_to_keep = *(params.chrs_to_keep.begin()); // Get first chromosome in list (there should only be one).
	filter_sites_by_position(chr_to_keep, params.start_pos, params.end_pos);
	filter_sites_by_positions(params.positions_file, params.exclude_positions_file);
	filter_sites_by_BED_file(params.BED_file, params.BED_exclude);
	filter_sites_by_number_of_alleles(params.min_alleles, params.max_alleles);
	filter_sites_by_INFO_flags(params.site_INFO_flags_to_remove, params.site_INFO_flags_to_keep);
	filter_sites_by_quality(params.min_quality);
	filter_sites_by_mean_depth(params.min_mean_depth, params.max_mean_depth);
	filter_sites_by_mask(params.mask_file, params.invert_mask, params.min_kept_mask_value);
	filter_individuals_by_mean_depth(params.min_indv_mean_depth, params.max_indv_mean_depth);
	if (params.phased_only == true)
	{
		filter_individuals_by_phase();
		filter_sites_by_phase();
	}
	filter_genotypes_by_quality(params.min_genotype_quality);
	filter_genotypes_by_depth(params.min_genotype_depth, params.max_genotype_depth);
	filter_genotypes_by_filter_flag(params.geno_filter_flags_to_exclude, params.remove_all_filtered_genotypes);
	filter_individuals_by_call_rate(params.min_indv_call_rate);
	filter_individuals_randomly(params.max_N_indv);
	filter_sites_by_frequency_and_call_rate(params.min_maf, params.max_maf, params.min_non_ref_af, params.max_non_ref_af, params.min_site_call_rate);
	filter_sites_by_allele_count(params.min_mac, params.max_mac, params.min_non_ref_ac, params.max_non_ref_ac, params.max_missing_call_count);
	filter_sites_by_HWE_pvalue(params.min_HWE_pvalue);
	filter_sites_by_thinning(params.min_interSNP_distance);
}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_entry.size(); ui++)
		if (include_entry[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_indv() const
{
	return N_indv;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_default_contigs(vector<string> &contig_vector)
{
	contig_vector.resize(0);
	contig_vector.push_back("##contig=<ID=1,length=249250621,assembly=b37>");
	contig_vector.push_back("##contig=<ID=2,length=243199373,assembly=b37>");
	contig_vector.push_back("##contig=<ID=3,length=198022430,assembly=b37>");
	contig_vector.push_back("##contig=<ID=4,length=191154276,assembly=b37>");
	contig_vector.push_back("##contig=<ID=5,length=180915260,assembly=b37>");
	contig_vector.push_back("##contig=<ID=6,length=171115067,assembly=b37>");
	contig_vector.push_back("##contig=<ID=7,length=159138663,assembly=b37>");
	contig_vector.push_back("##contig=<ID=8,length=146364022,assembly=b37>");
	contig_vector.push_back("##contig=<ID=9,length=141213431,assembly=b37>");
	contig_vector.push_back("##contig=<ID=10,length=135534747,assembly=b37>");
	contig_vector.push_back("##contig=<ID=11,length=135006516,assembly=b37>");
	contig_vector.push_back("##contig=<ID=12,length=133851895,assembly=b37>");
	contig_vector.push_back("##contig=<ID=13,length=115169878,assembly=b37>");
	contig_vector.push_back("##contig=<ID=14,length=107349540,assembly=b37>");
	contig_vector.push_back("##contig=<ID=15,length=102531392,assembly=b37>");
	contig_vector.push_back("##contig=<ID=16,length=90354753,assembly=b37>");
	contig_vector.push_back("##contig=<ID=17,length=81195210,assembly=b37>");
	contig_vector.push_back("##contig=<ID=18,length=78077248,assembly=b37>");
	contig_vector.push_back("##contig=<ID=19,length=59128983,assembly=b37>");
	contig_vector.push_back("##contig=<ID=20,length=63025520,assembly=b37>");
	contig_vector.push_back("##contig=<ID=21,length=48129895,assembly=b37>");
	contig_vector.push_back("##contig=<ID=22,length=51304566,assembly=b37>");
	contig_vector.push_back("##contig=<ID=X,length=155270560,assembly=b37>");
	contig_vector.push_back("##contig=<ID=Y,length=59373566,assembly=b37>");
	contig_vector.push_back("##contig=<ID=MT,length=16569,assembly=b37>");
}
