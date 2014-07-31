/*
 * variant_file_output.cpp
 *
 *  Created on: Aug 28, 2009
 *      Author: Adam Auton
 *      ($Revision: 249 $)
 */
#include "variant_file.h"

void variant_file::output_frequency(const string &output_file_prefix, bool output_counts, bool suppress_allele_output, bool derived)
{
	// Output statistics of frequency at each site
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Frequency Statistics.");

	LOG.printLOG("Outputting Frequency Statistics...\n");
	string output_file = output_file_prefix + ".frq";
	if (output_counts)
		output_file += ".count";

	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	if (suppress_allele_output == false)
	{
		out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:";
		if (output_counts)
			out << "COUNT}" << endl;
		else
			out << "FREQ}" << endl;
	}
	else
	{
		if (output_counts)
			out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}" << endl;
		else
			out << "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}" << endl;
	}

	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	unsigned int aa_idx = 0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);

		if (derived)
			e->parse_basic_entry(true, false, true);
		else
			e->parse_basic_entry(true);
		e->parse_genotype_entries(true);
		N_alleles = e->get_N_alleles();

		if (derived)
		{
			string AA = e->get_INFO_value("AA");
			std::transform(AA.begin(), AA.end(), AA.begin(), ::toupper);	// Comment this out if only want high quality sites.
			if ((AA == "?") || (AA == "."))
			{
				LOG.one_off_warning("\tWarning: Cannot output derived allele frequencies without Ancestral Alleles (AA)");
				continue;
			}
			else
			{
				bool found = false;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (AA == e->get_allele(ui))
					{
						aa_idx = ui;
						found = true;
						break;
					}
				}
				if (found == false)
				{
					LOG.one_off_warning("\tWarning: Ancestral allele does not match any SNP allele.");
					continue;
				}
			}
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << N_alleles << "\t" << N_non_missing_chr;
		if (output_counts)
		{
			if (suppress_allele_output == false)
			{
				out << "\t" << e->get_allele(aa_idx) << ":" << allele_counts[aa_idx];
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
						out << "\t" << e->get_allele(ui) << ":" << allele_counts[ui];
				}
				out << endl;
			}
			else
			{
				out << "\t" << allele_counts[aa_idx];
				for (unsigned ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
						out << "\t" << allele_counts[ui];
				}
				out << endl;
			}
		}
		else
		{
			double freq;
			if (suppress_allele_output == false)
			{
				freq = allele_counts[aa_idx] / (double)N_non_missing_chr;
				out << "\t" << e->get_allele(aa_idx) << ":" << freq;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
					{
						freq = allele_counts[ui] / (double)N_non_missing_chr;
						out << "\t" << e->get_allele(ui) << ":" << freq;
					}
				}
				out << endl;
			}
			else
			{
				freq = allele_counts[aa_idx] / (double)N_non_missing_chr;
				out << "\t" << freq;
				for (unsigned int ui=0; ui<N_alleles; ui++)
				{
					if (ui != aa_idx)
					{
						freq = allele_counts[ui] / (double)N_non_missing_chr;
						out << "\t" << freq;
					}
				}
				out << endl;
			}
		}
	}
	delete e;
	out.close();
}

void variant_file::output_het(const string &output_file_prefix)
{
	// Output statistics on Heterozygosity for each individual
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Heterozygosity Statistics.");
	// Following the calculations in PLINK....
	// Note this assumes Biallelic SNPs.

	LOG.printLOG("Outputting Individual Heterozygosity\n");

	string output_file = output_file_prefix + ".het";
	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF" << endl;

	// P(Homo) = F + (1-F)P(Homo by chance)
	// P(Homo by chance) = p^2+q^2 for a biallelic locus.
	// For an individual with N genotyped loci, we
	//   1. count the total observed number of loci which are homozygous (O),
	//   2. calculate the total expected number of loci homozygous by chance (E)
	// Then, using the method of moments, we have
	//    O = NF + (1-F)E
	// Which rearranges to give
	//    F = (O-E)/(N-E)

	// First, calc frequency of each site (should really move this to a subroutine)
	vector<double> freq(N_entries, 0.0);
	vector<int> allele_counts;
	vector<unsigned int> N_non_missing_chr(N_entries,0);
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tIndividual Heterozygosity: Only using biallelic SNPs.");
			continue;
		}

		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tIndividual Heterozygosity: Only using fully diploid SNPs.");
			continue;
		}

		// Frequency of non-reference allele
		e->get_allele_counts(allele_counts, N_non_missing_chr[s], include_indv, include_genotype[s]);

		if (N_non_missing_chr[s] > 0)
			freq[s] = allele_counts[1] / double(N_non_missing_chr[s]);
		else
			freq[s] = -1;
	}

	vector<int> N_sites_included(N_indv, 0);
	vector<int> N_obs_hom(N_indv, 0);
	vector<double> N_expected_hom(N_indv, 0.0);
	pair<int, int> alleles;

	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
			continue;

		e->parse_genotype_entries(true);
		if (e->is_diploid(include_indv, include_genotype[s]) == false)
			continue;

		if ((freq[s] <= numeric_limits<double>::epsilon())  || (1.0 - freq[s] <= numeric_limits<double>::epsilon()))
			continue;

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e->get_indv_GENOTYPE_ids(ui, alleles);
				if ((alleles.first != -1) && (alleles.second != -1))
				{
					N_sites_included[ui]++;
					if (alleles.first == alleles.second)
						N_obs_hom[ui]++;

					/////////////////////////
					// Expected homozygosity
					// E = 1 - (2pq . 2N/(2N-1))
					// (Using Nei's unbiased estimator)
					N_expected_hom[ui] += 1.0 - (2.0 * freq[s] * (1.0 - freq[s]) * (N_non_missing_chr[s] / (N_non_missing_chr[s] - 1.0)));
				}
			}
		}
	}

	out.setf(ios::fixed,ios::floatfield);
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		if (N_sites_included[ui] > 0)
		{
			double F = (N_obs_hom[ui] - N_expected_hom[ui]) / double(N_sites_included[ui] - N_expected_hom[ui]);
			out << indv[ui] << "\t" << N_obs_hom[ui] << "\t";
			out.precision(1);
			out << N_expected_hom[ui] << "\t";
			out.precision(5);
			out << N_sites_included[ui] << "\t" << F << endl;
		}
	}
	delete e;
	out.close();
}

void variant_file::output_hwe(const string &output_file_prefix)
{
	// Output HWE statistics for each site as described in Wigginton, Cutler, and Abecasis (2005)
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output HWE Statistics.");
	// Note this assumes Biallelic SNPs.
	LOG.printLOG("Outputting HWE statistics (but only for biallelic loci)\n");

	string output_file = output_file_prefix + ".hwe";
	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq\tP" << endl;

	/* PLINK code:
	// b11 = Nhom1, b12 = Nhet, b22 = Nhom2
	double tot = b11 + b12 + b22;
	double exp_11 = freq * freq * tot;
	double exp_12 = 2 * freq * (1-freq) * tot;
	double exp_22 = (1-freq) * (1-freq) * tot;

	double chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
			    + ( (b12-exp_12)*(b12-exp_12) ) / exp_12
			    + ( (b22-exp_22)*(b22-exp_22) ) / exp_22 ;

	p = chiprobP(chisq,1);
	*/

	double freq;
	unsigned int b11, b12, b22;
	double exp_11, exp_12, exp_22;
	double chisq;
	double tot;
	double p;
	unsigned int precision = out.precision();
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tHWE: Only using biallelic SNPs.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tHWE: Only using fully diploid SNPs.");
			continue;	// Isn't diploid
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		freq = allele_counts[0] / (double)N_non_missing_chr;
		e->get_genotype_counts(include_indv, include_genotype[s], b11, b12, b22);
		tot = b11 + b12 + b22;
		exp_11 = freq * freq * tot;
		exp_12 = 2.0 * freq * (1.0-freq) * tot;
		exp_22 = (1.0-freq) * (1.0-freq) * tot;

		chisq = ( (b11-exp_11)*(b11-exp_11) ) / exp_11
				+ ( (b12-exp_12)*(b12-exp_12) ) / exp_12
				+ ( (b22-exp_22)*(b22-exp_22) ) / exp_22;

		p = entry::SNPHWE(b12, b11, b22);
		out << e->get_CHROM() << "\t" << e->get_POS();
		out << "\t" << b11 << "/" << b12 << "/" << b22;
		out.precision(2);
		out << fixed << "\t" << exp_11 << "/" << exp_12 << "/" << exp_22;
		out.precision(precision);
		out << "\t" << chisq << "\t" << p << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_individuals_by_mean_depth(const string &output_file_prefix)
{
	// Output information regarding the mean depth for each individual
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Individuals by Mean Depth Statistics.");

	LOG.printLOG("Outputting Mean Depth by Individual\n");
	string output = output_file_prefix + ".idepth";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Individual Depth Output File: " + output, 2);
	out << "INDV\tN_SITES\tMEAN_DEPTH" << endl;
	vector<double> depth_sum(N_indv, 0.0);
	vector<int> count(N_indv, 0);
	int depth;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e->parse_genotype_entry(ui, false, false, true);
				depth = e->get_indv_DEPTH(ui);
				if (depth >= 0)
				{
					depth_sum[ui] += depth;
					count[ui]++;
				}
			}
		}
	}

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		double mean_depth = depth_sum[ui] / count[ui];
		out << indv[ui] << "\t" << count[ui] << "\t" << mean_depth << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_SNP_density(const string &output_file_prefix, int bin_size)
{
	// Output SNP density (technically variant density)
	if (bin_size <= 0)
		return;
	LOG.printLOG("Outputting SNP density\n");

	string output = output_file_prefix + ".snpden";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open SNP Density Output File: " + output, 2);

	// Find maximum position
	map<string, int> max_pos;
	map<string, int> min_pos;

	string CHROM; int POS;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			set_filepos(entry_file_locations[s]);
			read_CHROM_and_POS_only(CHROM, POS);
			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (POS > max_pos[CHROM])
					max_pos[CHROM] = POS;
			}
			else
				max_pos[CHROM] = POS;

			if (min_pos.find(CHROM) != min_pos.end())
			{
				if (POS < min_pos[CHROM])
					min_pos[CHROM] = POS;
			}
			else
				min_pos[CHROM] = POS;
		}
	}

	map<string, int>::iterator it;

	unsigned int N_bins;
	map<string, vector<int> > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int)((max_pos[CHROM] + bin_size) / double(bin_size));
		bins[CHROM].resize(N_bins, 0);
	}

	unsigned int idx;
	double C = 1.0 / double(bin_size);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			set_filepos(entry_file_locations[s]);
			read_CHROM_and_POS_only(CHROM, POS);
			idx = (unsigned int)(POS * C);
			bins[CHROM][idx]++;
		}
	}

	out << "CHROM\tBIN_START\tSNP_COUNT\tVARIANTS/KB" << endl;
	double sum1=0.0, sum2=0.0;
	int bin_tot;
	C = 1000.0 / bin_size;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		bool output = false;
		CHROM = (*it).first;
		sum2 += (max_pos[CHROM] - min_pos[CHROM]);
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			bin_tot = bins[CHROM][s];
			sum1 += bin_tot;
			if (bin_tot > 0)
				output = true;
			if (output == true)
				out << CHROM << "\t" << s*bin_size << "\t" << bin_tot << "\t" << bin_tot * C << endl;
		}
	}
	delete e;
	out.close();

	double mean_SNP_density = sum1 / sum2 * 1000;
	LOG.printLOG("Mean SNP density: " + output_log::dbl2str(mean_SNP_density, 5) + " variants / kb\n");
}

void variant_file::output_missingness(const string &output_file_prefix)
{
	// Output missingness by individual and site
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Missingness Statistics.");

	LOG.printLOG("Outputting Site and Individual Missingness\n");
	string output1 = output_file_prefix + ".imiss";
	ofstream out1(output1.c_str());
	if (!out1.is_open())
		LOG.error("Could not open Individual Missingness Output File: " + output1, 3);

	string output2 = output_file_prefix + ".lmiss";
	ofstream out2(output2.c_str());
	if (!out2.is_open())
		LOG.error("Could not open Site Missingness Output File: " + output2, 4);

	out1 << "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS" << endl;
	unsigned int ui, s;
	vector<unsigned int> indv_N_missing(N_indv, 0), indv_N_tot(N_indv, 0);
	vector<unsigned int> indv_N_geno_filtered(N_indv, 0);
	unsigned int site_N_missing, site_N_tot, site_N_geno_filtered;
	pair<int, int> alleles;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);

	out2 << "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS" << endl;
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry();

		site_N_missing = 0;
		site_N_tot = 0;
		site_N_geno_filtered = 0;
		for (ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (include_genotype[s][ui] == false)
			{
				site_N_geno_filtered++;
				indv_N_geno_filtered[ui]++;
				continue;
			}

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);
			if (alleles.first == -1)
			{
				site_N_missing++;
				indv_N_missing[ui]++;
			}
			indv_N_tot[ui]++;

			if (alleles.second == -1)
			{
				site_N_missing++;
			}
			site_N_tot+=2;

			if ((alleles.second == -1) && (e->get_indv_PHASE(ui) == '|'))
			{	// Phased missing genotypes indicate haploid genome
				site_N_tot--;
				site_N_missing--;
			}
		}
		out2 << e->get_CHROM() << "\t" << e->get_POS() << "\t" << site_N_tot << "\t" << site_N_geno_filtered << "\t";
		out2 << site_N_missing << "\t" << double(site_N_missing) / double(site_N_tot) << endl;
	}

	for (ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		out1 << indv[ui] << "\t" << indv_N_tot[ui] << "\t";
		out1 << indv_N_geno_filtered[ui] << "\t" << indv_N_missing[ui] << "\t";
		out1 << indv_N_missing[ui] / double(indv_N_tot[ui]) << endl;
	}

	delete e;
	out2.close();
	out1.close();
}

void variant_file::calc_hap_r2(entry *e, entry *e2, const vector<bool> &include_geno1, const vector<bool> &include_geno2, double &r2, double &D, double &Dprime, int &chr_count)
{
	double x11=0, x12=0, x21=0, x22=0;
	double X=0, X2=0, Y=0, Y2=0, XY=0;
	double sx, sy;
	double rel_x11, p1, p2, q1, q2, Dmax;
	double var1, var2, cov12;
	chr_count = 0;
	pair<int, int> geno1, geno2;
	int allele1, allele2;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == false) || (include_geno1[ui] == false) || (include_geno2[ui] == false))
			continue;

		e->get_indv_GENOTYPE_ids(ui, geno1);

		e2->parse_genotype_entry(ui, true);
		e2->get_indv_GENOTYPE_ids(ui, geno2);

//		if ((e->get_indv_ploidy(ui) != 2) || (e2->get_indv_ploidy(ui) != 2))
//		{
//			LOG.one_off_warning("\tLD: Only using diploid individuals.");
//			continue;
//		}

		if ((e->get_indv_PHASE(ui) != '|') || (e2->get_indv_PHASE(ui) != '|'))
			LOG.error("Require phased haplotypes for r^2 calculation (use --phased)\n");

		for (unsigned int c=0; c<2; c++)
		{
			if (c==0)
			{
				allele1 = geno1.first;
				allele2 = geno2.first;
			}
			else
			{
				allele1 = geno1.second;
				allele2 = geno2.second;
			}

			if ((allele1 == -1) || (allele2 == -1))
				continue;

			if (allele1 == 0 && allele2 == 0){
			  x11++;
			} else if (allele1 == 0 && allele2 != 0){
			  x12++;
			} else if (allele1 != 0 && allele2 == 0){
			  x21++;
			} else { // (allele1 !=0 && allele2 != 0)
			  x22++;
			}

			sx=0, sy=0;
			if (allele1 == 0)
				sx += 1;

			if (allele2 == 0)
				sy += 1;

			X += sx; Y += sy;
			XY += sx*sy;
			sx *= sx; sy *= sy;
			X2 += sx;
			Y2 += sy;

			chr_count++;
		}
	}

	rel_x11 = x11/double(chr_count);
	p1 = (x11 + x12)/double(chr_count);
	p2 = (x21 + x22)/double(chr_count);
	q1 = (x11 + x21)/double(chr_count);
	q2 = (x12 + x22)/double(chr_count);
	D = rel_x11 - p1*q1;
	if (D < 0)
		Dmax = min(p1*q1,p2*q2);
	else
		Dmax = min(p1*q2,p2*q1);
	Dprime = D/Dmax;

	X /= chr_count; X2 /= chr_count;
	Y /= chr_count;	Y2 /= chr_count;
	XY /= chr_count;

	var1 = X2 - X*X;
	var2 = Y2 - Y*Y;
	cov12 = XY - X*Y;

	r2 = cov12 * cov12 / (var1 * var2);
}

// Calculate r2 for either haplotypes or genotypes using the em algorithm...
void variant_file::calc_r2_em(entry *e, entry *e2, const vector<bool> &include_geno1, const vector<bool> &include_geno2, double &r2, int &indv_count)
{
	r2 = 0;
	indv_count = 0;
	pair<int, int> geno1, geno2;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == false) || (include_geno1[ui] == false) || (include_geno2[ui] == false))
			continue;

		e->get_indv_GENOTYPE_ids(ui, geno1);
		e2->parse_genotype_entry(ui, true);
		e2->get_indv_GENOTYPE_ids(ui, geno2);
		// TODO... not yet implemented...!
		LOG.error("Not yet implmented!\n");
	}

}

void variant_file::calc_geno_r2(entry *e, entry *e2, const vector<bool> &include_geno1, const vector<bool> &include_geno2, double &r2, int &indv_count)
{
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");
	double X=0, X2=0, Y=0, Y2=0, XY=0;
	double sx, sy;
	indv_count = 0;
	pair<int, int> geno1, geno2;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == false) || (include_geno1[ui] == false) || (include_geno2[ui] == false))
			continue;

		e->get_indv_GENOTYPE_ids(ui, geno1);

		e2->parse_genotype_entry(ui, true);
		e2->get_indv_GENOTYPE_ids(ui, geno2);

		if ((e->get_indv_ploidy(ui) != 2) || (e2->get_indv_ploidy(ui) != 2))
		{
			LOG.one_off_warning("\tgenoLD: Only using diploid individuals.");
			continue;
		}

		if ((geno1.first == -1) || (geno1.second == -1))
			continue;

		if ((geno2.first == -1) || (geno2.second == -1))
			continue;

		sx=0, sy=0;
		if (geno1.first == geno1.second)
		{
			if (geno1.first == 0)
			{
				sx = 2;
			}
		}
		else
			sx = 1;

		if (geno2.first == geno2.second)
		{
			if (geno2.first == 0)
			{
				sy = 2;
			}
		}
		else
			sy = 1;

		X += sx; Y += sy;
		XY += sx*sy;
		sx *= sx; sy *= sy;
		X2 += sx; Y2 += sy;

		indv_count++;
	}

	X /= indv_count; X2 /= indv_count;
	Y /= indv_count; Y2 /= indv_count;
	XY /= indv_count;

	double var1 = X2 - X*X;
	double var2 = Y2 - Y*Y;
	double cov12 = XY - X*Y;

	r2 = cov12 * cov12 / (var1 * var2);
}

void variant_file::calc_geno_chisq(entry *e, entry *e2, const vector<bool> &include_geno1, const vector<bool> &include_geno2, double &chisq, double &dof, double &pval, int &indv_count)
{
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	int N0 = e->get_N_alleles();
	int N1 = e2->get_N_alleles();

	int N_genotypes0 = N0 * (N0+1) / 2;
	int N_genotypes1 = N1 * (N1+1) / 2;

	vector<vector<double> > observed(N_genotypes0, vector<double>(N_genotypes1,0));
	indv_count = 0;
	pair<int, int> geno1, geno2;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ((include_indv[ui] == false) || (include_geno1[ui] == false) || (include_geno2[ui] == false))
			continue;

		e->get_indv_GENOTYPE_ids(ui, geno1);

		e2->parse_genotype_entry(ui, true);
		e2->get_indv_GENOTYPE_ids(ui, geno2);

		if ((e->get_indv_ploidy(ui) != 2) || (e2->get_indv_ploidy(ui) != 2))
		{
			LOG.one_off_warning("\tgenoLD: Only using diploid individuals.");
			continue;
		}

		if ((geno1.first == -1) || (geno1.second == -1))
			continue;

		if ((geno2.first == -1) || (geno2.second == -1))
			continue;

		map<pair<int, int>, int> idx_lookup1;
		int count = 0;
		for (int uj=0; uj<N0; uj++)
		{
			for (int uk=uj; uk<N0; uk++)
			{
				idx_lookup1[make_pair(uj,uk)] = count;
				idx_lookup1[make_pair(uk,uj)] = count;
				count++;
			}
		}
		map<pair<int, int>, int> idx_lookup2;
		count = 0;
		for (int uj=0; uj<N1; uj++)
		{
			for (int uk=uj; uk<N1; uk++)
			{
				idx_lookup2[make_pair(uj,uk)] = count;
				idx_lookup2[make_pair(uk,uj)] = count;
				count++;
			}
		}

		int idx1 = idx_lookup1[geno1];
		int idx2 = idx_lookup2[geno2];

		observed[idx1][idx2]++;
		indv_count++;
	}

	vector<vector<double> > expected(N_genotypes0, vector<double>(N_genotypes1,0));
	vector<double> row_tot(N_genotypes0, 0);
	vector<double> col_tot(N_genotypes1, 0);
	double tot=0;
	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			row_tot[ui] += observed[ui][uj];
			col_tot[uj] += observed[ui][uj];
			tot += observed[ui][uj];
		}
	}

	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			expected[ui][uj] = row_tot[ui] * col_tot[uj] / tot;
		}
	}
/*
	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			cout << " " << observed[ui][uj];
		}
		cout << endl;
	}
	cout << endl << endl;

	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			cout << " " << expected[ui][uj];
		}
		cout << endl;
	}
	*/

	chisq = 0;
	for (int ui=0; ui<N_genotypes0; ui++)
	{
		for (int uj=0; uj<N_genotypes1; uj++)
		{
			if ((row_tot[ui] > 0) && (col_tot[uj] > 0))	// Don't use incomplete cases
				chisq += pow(observed[ui][uj] - expected[ui][uj], 2) / expected[ui][uj];
		}
	}
	int n_col=0, n_row=0;
	for (int ui=0; ui<N_genotypes0; ui++)
		if (row_tot[ui] > 0)
			n_row++;
	for (int ui=0; ui<N_genotypes1; ui++)
		if (col_tot[ui] > 0)
			n_col++;

	dof = (n_row-1) * (n_col-1);
	//cout << chisq << " " << dof << endl;
	pval = 1.0-gammp(dof/2, chisq/2);
}

void variant_file::output_haplotype_r2(const string &output_file_prefix, int snp_window_size, int snp_window_min, int bp_window_size, int bp_window_min, double min_r2)
{
	// Output pairwise LD statistics, using traditional r^2. Requires phased haplotypes.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;

	LOG.printLOG("Outputting Pairwise LD (phased bi-allelic only)\n");
	string output = output_file_prefix + ".hap.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime" << endl;

	double r2, D, Dprime;
	int chr_count;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = s+skip; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (int(s2 - s) > snp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e->get_CHROM() != e2->get_CHROM())
			{
				s2 = N_entries;	// No need to go any further (assuming SNPs are sorted)
				continue;
			}

			if ((e2->get_POS() - e->get_POS()) < bp_window_min)
					continue;

			if ((e2->get_POS() - e->get_POS()) > bp_window_size)
			{
				s2 = N_entries;	// No need to go any further (assuming SNPs are sorted)
				continue;
			}

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tLD: Only using biallelic variants.");
				continue;
			}

			calc_hap_r2(e, e2, include_genotype[s], include_genotype[s2], r2, D, Dprime, chr_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_POS() << "\t" << chr_count << "\t" << r2 << "\t" << D << "\t" << Dprime << "\t" << endl;
		}
	}
	delete e;
	delete e2;
	out.close();
}

void variant_file::output_genotype_r2(const string &output_file_prefix, int snp_window_size, int snp_window_min, int bp_window_size, int bp_window_min, double min_r2)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;

	LOG.printLOG("Outputting Pairwise LD (bi-allelic only)\n");
	string output = output_file_prefix + ".geno.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR\tPOS1\tPOS2\tN_INDV\tR^2" << endl;

	double r2;
	int indv_count;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tgenoLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = s+skip; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (int(s2 - s) > snp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tgenoLD: Only using biallelic variants.");
				continue;	// Isn't biallelic
			}

			if (e->get_CHROM() != e2->get_CHROM())
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			if ((e2->get_POS() - e->get_POS()) < bp_window_min)
				continue;

			if ((e2->get_POS() - e->get_POS()) > bp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			calc_geno_r2(e, e2, include_genotype[s], include_genotype[s2], r2, indv_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_POS() << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	delete e;
	delete e2;
	out.close();
}

void variant_file::output_genotype_chisq(const string &output_file_prefix, int snp_window_size, int snp_window_min, int bp_window_size, int bp_window_min, double min_pval)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;

	LOG.printLOG("Outputting Pairwise LD\n");
	string output = output_file_prefix + ".geno.chisq";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR\tPOS1\tPOS2\tN_INDV\tCHI^2\tDOF\tPVAL" << endl;

	double chisq, dof, pval;
	int indv_count;
	unsigned int skip = (unsigned int)max((int)1, snp_window_min);
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);

		for (s2 = s+skip; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (int(s2 - s) > snp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e->get_CHROM() != e2->get_CHROM())
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			if ((e2->get_POS() - e->get_POS()) < bp_window_min)
				continue;

			if ((e2->get_POS() - e->get_POS()) > bp_window_size)
			{
				s2 = N_entries;	// SNPs sorted, so no need to go any further
				continue;
			}

			calc_geno_chisq(e, e2, include_genotype[s], include_genotype[s2], chisq, dof, pval, indv_count);

			if (min_pval > 0)
				if ((pval < min_pval) | (pval != pval))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_POS() << "\t" << indv_count << "\t" << chisq << "\t" << dof << "\t" << pval << endl;
		}
	}
	delete e;
	delete e2;
	out.close();
}

void variant_file::output_interchromosomal_genotype_r2(const string &output_file_prefix, double min_r2)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;

	LOG.printLOG("Outputting Interchromosomal Pairwise Genotype LD (bi-allelic only)\n");
	string output = output_file_prefix + ".interchrom.geno.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2" << endl;

	int indv_count;
	double r2;
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = s+1; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
				continue;	// Isn't biallelic
			}

			if (e->get_CHROM() == e2->get_CHROM())
				continue;

			calc_geno_r2(e, e2, include_genotype[s], include_genotype[s2], r2, indv_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_CHROM() << "\t" << e2->get_POS() << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	delete e;
	delete e2;
	out.close();
}

void variant_file::output_interchromosomal_haplotype_r2(const string &output_file_prefix, double min_r2)
{
	// Output pairwise LD statistics, using genotype r^2. This is the same formula as used by PLINK, and is basically the squared
	// correlation coefficient between genotypes numbered as 0, 1, 2.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	unsigned int s, s2;

	LOG.printLOG("Outputting Interchromosomal Pairwise LD (bi-allelic only)\n");
	string output = output_file_prefix + ".interchrom.hap.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2" << endl;

	double D, Dprime;
	int chr_count;
	double r2;
	vector<char> variant_line, variant_line2;
	entry *e, *e2;
	e = get_entry_object(N_indv);
	e2 = get_entry_object(N_indv);
	for (s=0; s<(N_entries-1); s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = s+1; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
				continue;	// Isn't biallelic
			}

			if (e->get_CHROM() == e2->get_CHROM())
				continue;

			calc_hap_r2(e, e2, include_genotype[s], include_genotype[s2], r2, D, Dprime, chr_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_CHROM() << "\t" << e2->get_POS() << "\t" << chr_count << "\t" << r2 << endl;
		}
	}
	delete e;
	delete e2;
	out.close();
}

void variant_file::output_haplotype_r2_of_SNP_list_vs_all_others(const string &output_file_prefix, const string &positions_file, double min_r2)
{
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
			LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	LOG.printLOG("Outputting haplotype pairwise LD (bi-allelic only) for a set of SNPs verses all others.\n");

	vector< set<int > > keep_positions;
	map<string, int> chr_to_idx;
	string line;
	stringstream ss;
	string chr;
	int pos1, idx;
	unsigned int N_chr=0;

	ifstream BED(positions_file.c_str());
	if (!BED.is_open())
		LOG.error("Could not open Positions file: " + positions_file);
		// Skip header
	BED.ignore(numeric_limits<streamsize>::max(), '\n');
	while (!BED.eof())
	{
		getline(BED, line);
		if (line[0] == '#')
			continue;

		ss.clear();
		ss.str(line);
		ss >> chr >> pos1;

		if (chr_to_idx.find(chr) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[chr] = (N_chr-1);
			keep_positions.resize(N_chr);
		}

		idx = chr_to_idx[chr];
		keep_positions[idx].insert(pos1);
	}
	BED.close();

	unsigned int s, s2;
	string output = output_file_prefix + ".list.hap.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2" << endl;

	double D, Dprime;
	int chr_count;
	double r2;
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		e->get_CHROM(chr);
		if (chr_to_idx.find(chr) == chr_to_idx.end())
			continue;

		idx = chr_to_idx[chr];
		pos1 = e->get_POS();

		if (keep_positions[idx].find(pos1) == keep_positions[idx].end())
			continue;

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = 0; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;

			if (s == s2)
				continue;

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
				continue;	// Isn't biallelic
			}

			calc_hap_r2(e, e2, include_genotype[s], include_genotype[s2], r2, D, Dprime, chr_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_CHROM() << "\t" << e2->get_POS() << "\t" << chr_count << "\t" << r2 << endl;
		}
	}
	out.close();
}

void variant_file::output_genotype_r2_of_SNP_list_vs_all_others(const string &output_file_prefix, const string &positions_file, double min_r2)
{
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
			LOG.error("Require Genotypes in VCF file in order to output LD Statistics.");

	LOG.printLOG("Outputting genotype pairwise LD (bi-allelic only) for a set of SNPs verses all others.\n");

	vector< set<int > > keep_positions;
	map<string, int> chr_to_idx;
	string line;
	stringstream ss;
	string chr;
	int pos1, idx;
	unsigned int N_chr=0;

	ifstream BED(positions_file.c_str());
	if (!BED.is_open())
		LOG.error("Could not open Positions file: " + positions_file);
		// Skip header
	BED.ignore(numeric_limits<streamsize>::max(), '\n');
	while (!BED.eof())
	{
		getline(BED, line);
		if (line[0] == '#')
			continue;

		ss.clear();
		ss.str(line);
		ss >> chr >> pos1;

		if (chr_to_idx.find(chr) == chr_to_idx.end())
		{
			N_chr++;
			chr_to_idx[chr] = (N_chr-1);
			keep_positions.resize(N_chr);
		}

		idx = chr_to_idx[chr];
		keep_positions[idx].insert(pos1);
	}
	BED.close();

	unsigned int s, s2;
	string output = output_file_prefix + ".list.geno.ld";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open LD Output File: " + output, 3);

	out << "CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2" << endl;

	int indv_count;
	double r2;
	vector<char> variant_line, variant_line2;
	entry *e = get_entry_object(N_indv);
	entry *e2 = get_entry_object(N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		e->get_CHROM(chr);
		if (chr_to_idx.find(chr) == chr_to_idx.end())
			continue;

		idx = chr_to_idx[chr];
		pos1 = e->get_POS();

		if (keep_positions[idx].find(pos1) == keep_positions[idx].end())
			continue;

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
			continue;	// Isn't biallelic
		}

		e->parse_genotype_entries(true);

		for (s2 = 0; s2<N_entries; s2++)
		{
			if (include_entry[s2] == false)
				continue;
			if (s == s2)
				continue;

			get_entry(s2, variant_line2);
			e2->reset(variant_line2);
			e2->parse_basic_entry(true);

			if (e2->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tinterchromLD: Only using biallelic variants.");
				continue;	// Isn't biallelic
			}

			calc_geno_r2(e, e2, include_genotype[s], include_genotype[s2], r2, indv_count);

			if (min_r2 > 0)
				if ((r2 < min_r2) | (r2 != r2))
					continue;

			out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e2->get_CHROM() << "\t" << e2->get_POS() << "\t" << indv_count << "\t" << r2 << endl;
		}
	}
	out.close();
}

void variant_file::output_singletons(const string &output_file_prefix)
{
	// Locate and output singletons (and private doubletons)
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Singletons.");

	LOG.printLOG("Outputting Singleton Locations\n");
	string output = output_file_prefix + ".singletons";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Singleton Output File: " + output, 3);

	out << "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV" << endl;

	unsigned int ui;
	int a;
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	pair<int, int> geno;
	string allele;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);
		e->parse_genotype_entries(true);

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		N_alleles = e->get_N_alleles();

		for (a=0; a<(signed)N_alleles; a++)
		{
			if (allele_counts[a] == 1)
			{	// Singleton
				for (ui=0; ui<N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e->get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) || (geno.second == a))
					{
						e->get_allele(a, allele);
						out << e->get_CHROM() << "\t" << e->get_POS() << "\tS\t" << allele << "\t" << indv[ui] << endl;
						ui=N_indv;
						break;
					}
				}
			}
			else if (allele_counts[a] == 2)
			{	// Possible doubleton
				for (ui=0; ui<N_indv; ui++)
				{
					if (include_indv[ui] == false)
						continue;
					e->get_indv_GENOTYPE_ids(ui, geno);
					if ((geno.first == a) && (geno.second == a))
					{
						e->get_allele(a, allele);
						out << e->get_CHROM() << "\t" << e->get_POS() << "\tD\t" << allele << "\t" << indv[ui] << endl;
						ui=N_indv;
						break;
					}
				}
			}
		}
	}
	delete e;
	out.close();
}

void variant_file::output_genotype_depth(const string &output_file_prefix)
{
	// Output genotype depth in tab-delimited format.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Genotype Depth Statistics.");

	LOG.printLOG("Outputting Depth for Each Genotype\n");
	string output = output_file_prefix + ".gdepth";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Genotype Depth Output File: " + output, 7);

	out << "CHROM\tPOS";
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		out << "\t" << indv[ui];
	}
	out << endl;

	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry();

		out << e->get_CHROM() << "\t" << e->get_POS();

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			if (include_genotype[s][ui] == true)
			{
				e->parse_genotype_entry(ui, false, false, true);
				out << "\t" << e->get_indv_DEPTH(ui);
			}
			else
				out << "\t-1";
		}
		out << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_FILTER_summary(const string &output_file_prefix)
{
	// Output a summary of sites in various FILTER categories.
	LOG.printLOG("Outputting Filter Summary (for bi-allelic loci only)\n");

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0;
	model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2;
	model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4;
	model_to_idx["GT"] = 5;
	string FILTER;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);

	map<string, pair<int, int> > FILTER_to_TsTv;
	map<string, int > FILTER_to_Nsites;
	map<string, int >::iterator FILTER_to_Nsites_it;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true, true);

		string model = e->get_REF() + e->get_ALT_allele(0);
		sort(model.begin(), model.end());

		FILTER = e->get_FILTER();
		FILTER_to_Nsites[FILTER]++;
		if (model_to_idx.find(model) != model_to_idx.end())
		{
			switch (model_to_idx[model])
			{
			case 1:
			case 4:
				FILTER_to_TsTv[FILTER].first++;
				break;
			case 0:
			case 2:
			case 3:
			case 5:
				FILTER_to_TsTv[FILTER].second++;
				break;
			default:
				// Don't count this snp towards Ts/Tv
				break;
			}
		}
	}

	vector<pair<int, string > > count_to_FILTER;
	for ( FILTER_to_Nsites_it=FILTER_to_Nsites.begin() ; FILTER_to_Nsites_it != FILTER_to_Nsites.end(); ++FILTER_to_Nsites_it )
	{
		FILTER = (*FILTER_to_Nsites_it).first;
		int Nsites = (*FILTER_to_Nsites_it).second;

		count_to_FILTER.push_back(make_pair(Nsites, FILTER));
	}

	sort(count_to_FILTER.begin(), count_to_FILTER.end());

	string output = output_file_prefix + ".FILTER.summary";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Filter Summary Output File: " + output, 7);

	out << "FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv" << endl;

	for (int i=count_to_FILTER.size()-1; i > -1; i--)
	{
		FILTER = count_to_FILTER[i].second;
		int Ts = FILTER_to_TsTv[FILTER].first;
		int Tv = FILTER_to_TsTv[FILTER].second;
		int Nsites = FILTER_to_Nsites[FILTER];
		out << FILTER << "\t" << Nsites << "\t";
		out << Ts << "\t" << Tv << "\t" << double(Ts)/Tv << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_TsTv(const string &output_file_prefix, int bin_size)
{
	// Output Ts/Tv ratios in bins of a given size.
	LOG.printLOG("Outputting Ts/Tv in bins of " + output_log::int2str(bin_size) + "bp\n");

	map<string, unsigned int> model_to_idx;
	model_to_idx["AC"] = 0;
	model_to_idx["AG"] = 1;
	model_to_idx["AT"] = 2;
	model_to_idx["CG"] = 3;
	model_to_idx["CT"] = 4;
	model_to_idx["GT"] = 5;

	map<string, int> max_pos;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e->get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e->get_POS();
			}
			else
				max_pos[CHROM] = e->get_POS();
		}
	}

	map<string, int>::iterator it;

	unsigned int N_bins;
	map<string, vector<int> > Ts_counts;
	map<string, vector<int> > Tv_counts;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int)((max_pos[CHROM] + bin_size) / double(bin_size));
		Ts_counts[CHROM].resize(N_bins, 0);
		Tv_counts[CHROM].resize(N_bins, 0);
	}

	vector<unsigned int> model_counts(6,0);
	double C = 1.0 / double(bin_size);
	unsigned int idx;

	string model;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (!e->is_biallelic_SNP())
			continue;

		model = e->get_REF() + e->get_ALT_allele(0);
		sort(model.begin(), model.end());

		CHROM = e->get_CHROM();
		idx = (unsigned int)(e->get_POS() * C);

		if (model_to_idx.find(model) != model_to_idx.end())
		{
			model_counts[model_to_idx[model]]++;
			switch (model_to_idx[model])
			{
			case 1:
			case 4:
				Ts_counts[CHROM][idx]++;
				break;
			case 0:
			case 2:
			case 3:
			case 5:
				Tv_counts[CHROM][idx]++;
				break;
			default:
				LOG.error("Unknown idx\n");
				break;
			}
		}
		else
			LOG.warning("Unknown model type. Not a SNP? " + CHROM + ":" + output_log::int2str(e->get_POS()) +"\n");
	}

	string output = output_file_prefix + ".TsTv";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open TsTv Output File: " + output, 7);

	out << "CHROM\tBinStart\tSNP_count\tTs/Tv" << endl;
	double ratio;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		for (unsigned int s=0; s<Ts_counts[CHROM].size(); s++)
		{
			ratio = 0.0;
			if (Tv_counts[CHROM][s] != 0)
				ratio = double(Ts_counts[CHROM][s]) / Tv_counts[CHROM][s];
			out << CHROM << "\t" << s*bin_size << "\t" << Ts_counts[CHROM][s]+Tv_counts[CHROM][s] << "\t" << ratio << endl;
		}
	}
	out.close();

	output = output_file_prefix + ".TsTv.summary";
	out.open(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open TsTv Summary Output File: " + output, 7);

	out << "MODEL\tCOUNT" << endl;
	out << "AC\t" << model_counts[0] << endl;
	out << "AG\t" << model_counts[1] << endl;
	out << "AT\t" << model_counts[2] << endl;
	out << "CG\t" << model_counts[3] << endl;
	out << "CT\t" << model_counts[4] << endl;
	out << "GT\t" << model_counts[5] << endl;
	unsigned int Ts = model_counts[1] + model_counts[4];
	unsigned int Tv = model_counts[0] + model_counts[2] + model_counts[3] + model_counts[5];
	out << "Ts\t" << Ts << endl;
	out << "Tv\t" << Tv << endl;

	LOG.printLOG("Ts/Tv ratio: " + output_log::dbl2str(double(Ts)/Tv, 4) + "\n");

	delete e;
	out.close();
}

void variant_file::output_TsTv_by_count(const string &output_file_prefix)
{
	// Output Ts/Tv ratios in bins of a given size.
	LOG.printLOG("Outputting Ts/Tv by Alternative Allele Count\n");
	vector<unsigned int> Ts_counts, Tv_counts;
	unsigned int N_kept_indv = N_kept_individuals();
	Ts_counts.resize(2*N_kept_indv);
	Tv_counts.resize(2*N_kept_indv);

	string model;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	map<string, unsigned int> model_to_Ts_or_Tv;
	model_to_Ts_or_Tv["AC"] = 1;
	model_to_Ts_or_Tv["CA"] = 1;
	model_to_Ts_or_Tv["AG"] = 0;	// Ts
	model_to_Ts_or_Tv["GA"] = 0;	// Ts
	model_to_Ts_or_Tv["AT"] = 1;
	model_to_Ts_or_Tv["TA"] = 1;
	model_to_Ts_or_Tv["CG"] = 1;
	model_to_Ts_or_Tv["GC"] = 1;
	model_to_Ts_or_Tv["CT"] = 0;	// Ts
	model_to_Ts_or_Tv["TC"] = 0;	// Ts
	model_to_Ts_or_Tv["GT"] = 1;
	model_to_Ts_or_Tv["TG"] = 1;
	unsigned int idx;
	vector<int> allele_counts;
	unsigned int allele_count;
	unsigned int N_included_indv;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true);

			if (!e->is_biallelic_SNP())
				continue;

			e->parse_genotype_entries(true);
			e->get_allele_counts(allele_counts, N_included_indv, include_indv, include_genotype[s]);
			allele_count = allele_counts[1];

			model = e->get_REF() + e->get_ALT_allele(0);
			if (model_to_Ts_or_Tv.find(model) != model_to_Ts_or_Tv.end())
			{
				idx = model_to_Ts_or_Tv[model];
				if (idx == 0) // Ts
					Ts_counts[allele_count]++;
				else if (idx == 1) // Tv;
					Tv_counts[allele_count]++;
				else
					LOG.error("Unknown model type\n");
			}
			else
				LOG.warning("Unknown model type. Not a SNP? " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()) +"\n");
		}
	}

	string output = output_file_prefix + ".TsTv.count";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open TsTv by Count Output File: " + output, 7);

	double ratio;
	out << "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv" << endl;
	for (unsigned int ui=0; ui<2*N_kept_indv; ui++)
	{
		ratio = double(Ts_counts[ui]) / Tv_counts[ui];
		out << ui << "\t" << Ts_counts[ui] << "\t" << Tv_counts[ui] << "\t" << ratio << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_TsTv_by_quality(const string &output_file_prefix)
{
	// Output Ts/Tv ratios in bins of a given size.
	LOG.printLOG("Outputting Ts/Tv By Quality\n");
	map<double, pair<unsigned int, unsigned int> > TsTv_counts;
	double max_qual = -numeric_limits<double>::max(), min_qual=numeric_limits<double>::max();

	string model;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	map<string, unsigned int> model_to_Ts_or_Tv;
	model_to_Ts_or_Tv["AC"] = 1;
	model_to_Ts_or_Tv["CA"] = 1;
	model_to_Ts_or_Tv["AG"] = 0;	// Ts
	model_to_Ts_or_Tv["GA"] = 0;	// Ts
	model_to_Ts_or_Tv["AT"] = 1;
	model_to_Ts_or_Tv["TA"] = 1;
	model_to_Ts_or_Tv["CG"] = 1;
	model_to_Ts_or_Tv["GC"] = 1;
	model_to_Ts_or_Tv["CT"] = 0;	// Ts
	model_to_Ts_or_Tv["TC"] = 0;	// Ts
	model_to_Ts_or_Tv["GT"] = 1;
	model_to_Ts_or_Tv["TG"] = 1;
	unsigned int idx;
	double QUAL;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true);

			if (!e->is_biallelic_SNP())
				continue;

			QUAL = e->get_QUAL();
			if (QUAL > max_qual)
				max_qual = QUAL;
			if (QUAL < min_qual)
				min_qual = QUAL;

			model = e->get_REF() + e->get_ALT_allele(0);;
			if (model_to_Ts_or_Tv.find(model) != model_to_Ts_or_Tv.end())
			{
				idx = model_to_Ts_or_Tv[model];
				if (idx == 0) // Ts
				{
					TsTv_counts[QUAL].first++;
				}
				else if (idx == 1) // Tv;
					TsTv_counts[QUAL].second++;
				else
					LOG.error("Unknown model type\n");
			}
			else
				LOG.warning("Unknown model type. Not a SNP? " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()) +"\n");
		}
	}

	string output = output_file_prefix + ".TsTv.qual";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open TsTv by Count Output File: " + output, 7);

	out << "QUAL_THRESHOLD";
	out << "\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD";
	out << "\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD" << endl;

	unsigned int N_TsTv = TsTv_counts.size();

	vector<double> Ts_sum_below(N_TsTv+1, 0.0), Tv_sum_below(N_TsTv+1, 0.0);
	vector<double> QUAL_vector(N_TsTv+1, 0.0);
	QUAL_vector[0] = min_qual;
	QUAL_vector[N_TsTv] = max_qual;
	idx = 1;
	for (map<double, pair<unsigned int, unsigned int> >::iterator it=TsTv_counts.begin(); it != TsTv_counts.end(); ++it)
	{
		QUAL = (it->first);
		double Ts = (it->second).first;
		double Tv = (it->second).second;
		Ts_sum_below[idx] = Ts_sum_below[idx-1]+Ts;
		Tv_sum_below[idx] = Tv_sum_below[idx-1]+Tv;
		QUAL_vector[idx-1] = QUAL;
		idx++;
	}
	QUAL_vector[N_TsTv] = max_qual;

	vector<double> Ts_sum_above(N_TsTv+1, 0.0), Tv_sum_above(N_TsTv+1, 0.0);
	idx = N_TsTv;
	for (map<double, pair<unsigned int, unsigned int> >::reverse_iterator it=TsTv_counts.rbegin(); it != TsTv_counts.rend(); ++it)
	{
		QUAL = (it->first);
		double Ts = (it->second).first;
		double Tv = (it->second).second;
		Ts_sum_above[idx] = Ts_sum_above[idx+1]+Ts;
		Tv_sum_above[idx] = Tv_sum_above[idx+1]+Tv;
		idx--;
	}

	double Ts_sum, Tv_sum, ratio;
	for (unsigned int ui=1; ui<(N_TsTv+1); ui++)
	{
		QUAL = QUAL_vector[ui-1];
		out << QUAL;
		Ts_sum = Ts_sum_below[ui-1]; Tv_sum = Tv_sum_below[ui-1];
		ratio = Ts_sum / Tv_sum;
		out << "\t" << Ts_sum << "\t" << Tv_sum << "\t" << ratio;
		Ts_sum = Ts_sum_above[ui+1]; Tv_sum = Tv_sum_above[ui+1];
		ratio = Ts_sum / Tv_sum;
		out << "\t" << Ts_sum << "\t" << Tv_sum << "\t" << ratio;
		out << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_site_quality(const string &output_file_prefix)
{
	// Output per-site quality information.
	LOG.printLOG("Outputting Quality for Each Site\n");
	string output = output_file_prefix + ".lqual";

	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Site Depth Output File: " + output, 7);

	out << "CHROM\tPOS\tQUAL" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry();

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << e->get_QUAL() << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_site_depth(const string &output_file_prefix, bool output_mean)
{
	// Output per-site depth information
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Site Depth Statistics.");

	LOG.printLOG("Outputting Depth for Each Site\n");
	string output = output_file_prefix + ".ldepth";
	if (output_mean)
		output += ".mean";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Site Depth Output File: " + output, 7);

	out << "CHROM\tPOS\t";
	if (output_mean)
		out << "MEAN_DEPTH\tVAR_DEPTH" << endl;
	else
		out << "SUM_DEPTH\tSUMSQ_DEPTH" << endl;

	int depth;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry();

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t";

		unsigned int sum=0;
		unsigned int sumsq=0;
		unsigned int n=0;
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			if (include_genotype[s][ui] == false)
				continue;

			e->parse_genotype_entry(ui, false, false, true);
			depth = e->get_indv_DEPTH(ui);
			if (depth >= 0)
			{
				sum += depth;
				sumsq += (depth*depth);
				n++;
			}
		}

		if (output_mean)
		{
			double mean = double(sum) / n;
			double var = ((double(sumsq) / n) - (mean*mean)) * double(n) / double(n-1);
			out << mean << "\t" << var << endl;
		}
		else
			out << sum << "\t" << sumsq << endl;
	}
	delete e;
	out.close();
}

void variant_file::output_hapmap_fst(const string &output_file_prefix, const vector<string> &indv_files)
{
	// Calculate Fst using individuals in one (rather than two VCF files)
	// Calculate, and output, Fst using the formula outlined in HapMap I
	// Namely:
	// Fst = 1 - (Pi_within / Pi_combined)
	// where
	// Pi_within = sum_j(nchoosek(n_j,2) * sum_i(2*n_ij * x_ij * (1-x_ij) / (n_ij -1))) / sum_j(nchoosek(n_j,2))
	// and
	// Pi_between = sum_i(2*n_i*x_i*(1-x_i) / (n_i - 1))
	// where j is the population index, and i is the SNP index

	if (indv_files.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting HapMap-style Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = indv_files.size();
	indvs_in_pops.resize(N_pops, vector<bool>(N_indv, false));
	vector<bool> all_indv(N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(indv_files[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + indv_files[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	string output = output_file_prefix + ".hapmap.fst";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Fst Output File: " + output, 7);

	out << "CHROM\tPOS\tHAPMAP_FST" << endl;

	entry *e = get_entry_object(N_indv);
	vector<char> variant_line;
	vector<int> allele_counts1;
	double Fst_tot_num=0.0, Fst_tot_denom=0.0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tFst: Only using biallelic sites.");
			continue;
		}

		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		unsigned int N_chr;
		e->get_allele_counts(allele_counts1, N_chr, all_indv, include_genotype[s]);
		double count_all = allele_counts1[1];
		double N_chr_all = N_chr;

		if ((count_all == 0) || (count_all == N_chr_all))
			continue;	// No polymorphism

		vector<double> counts(N_pops, 0);
		vector<double> pop_N_chr(N_pops, 0);
		vector<double> pop_N_choose_2(N_pops, 0);
		for (unsigned int p=0; p<N_pops; p++)
		{
			e->get_allele_counts(allele_counts1, N_chr, indvs_in_pops[p], include_genotype[s]);
			counts[p] = allele_counts1[1];
			pop_N_chr[p] = N_chr;
			pop_N_choose_2[p] = N_chr * (N_chr-1.0) / 2.0;
		}

		double Fst_SNP = 0;
		double f;
		double sum1=0.0;
		for (unsigned int p=0; p<N_pops; p++)
		{
			f = counts[p] / pop_N_chr[p];
			Fst_SNP += 2.0*pop_N_choose_2[p]*(pop_N_chr[p]/(pop_N_chr[p]-1.0))*f*(1.0-f);
			sum1 += pop_N_choose_2[p];
		}
		Fst_SNP /= sum1;
		Fst_tot_num += Fst_SNP;
		f = count_all / N_chr_all;
		double tmp = (2.0*(N_chr_all / (N_chr_all-1.0))*f*(1.0-f));
		Fst_SNP /= tmp;
		Fst_tot_denom += tmp;
		Fst_SNP = 1.0 - Fst_SNP;
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << Fst_SNP << endl;
	}
	double Fst_tot = 1.0 - (Fst_tot_num / Fst_tot_denom);
	LOG.printLOG("HapMap-style Fst = " + output_log::dbl2str(Fst_tot, 6) + "\n");

	delete e;
	out.close();
}

void variant_file::output_weir_and_cockerham_fst(const string &output_file_prefix, const vector<string> &indv_files)
{	// Implements the bi-allelic version of Weir and Cockerham's Fst
	if (indv_files.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting Weir and Cockerham Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = indv_files.size();
	indvs_in_pops.resize(N_pops, vector<bool>(N_indv, false));
	vector<bool> all_indv(N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(indv_files[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + indv_files[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	string output = output_file_prefix + ".weir.fst";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Fst Output File: " + output, 7);

	out << "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST" << endl;

	entry *e = get_entry_object(N_indv);
	vector<char> variant_line;

	double snp_Fst;
	double sum1=0.0, sum2 = 0.0;
	double sum3=0.0, count = 0.0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tFst: Only using biallelic sites.");
			continue;
		}

		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tFst: Only using diploid sites.");
			continue;
		}

		vector<unsigned int> n;
		n.resize(N_pops, 0);
		vector<double> p;
		p.resize(N_pops, 0);
		double nbar = 0.0, pbar=0.0, hbar=0.0;
		double ssqr=0.0;
		double sum_nsqr = 0.0;
		double n_sum = 0.0;
		unsigned int N_hom1, N_het, N_hom2;

		for (unsigned int j=0; j<N_pops; j++)
		{
			e->get_genotype_counts(indvs_in_pops[j], include_genotype[s], N_hom1, N_het, N_hom2);
			n[j] = N_hom1 + N_het + N_hom2;
			hbar += N_het;
			p[j] = N_het + 2*N_hom2;
			nbar += n[j];
			pbar += p[j];
			sum_nsqr += (n[j] * n[j]);

			p[j] /= (2.0*n[j]);	// diploid
		}

		n_sum = nbar;
		nbar /= double(N_pops);
		hbar /= n_sum;
		pbar /= (n_sum*2.0);	// diploid

		for (unsigned int j=0; j<N_pops; j++)
			ssqr += (n[j]*(p[j] - pbar)*(p[j] - pbar));
		ssqr /= ((N_pops-1.0)*nbar);
		double nc = (n_sum - (sum_nsqr / n_sum)) / (N_pops - 1.0);
		double C2 = N_pops * (1.0 - (nc / nbar));
		double S1 = ssqr - ((pbar * (1.0 - pbar) - ((N_pops-1.0) * ssqr / N_pops) - (hbar/4.0)) / (nbar - 1.0));
		double S2 = pbar * (1.0-pbar) * (1.0 - (nbar*C2/(N_pops*(nbar-1.0))));
		S2 += ((ssqr/N_pops)*(1.0 + ((N_pops-1.0)*nbar*C2/(N_pops*(nbar-1.0)))));
		S2 += (C2/(N_pops*(nbar-1.0)))*hbar/4.0;

		snp_Fst = S1 / S2;
		/*
		// Calculate Fit and Fis
		double S3 = nc * hbar / nbar;
		double FIT = 1.0 - (S3 / (2*S2));	// FIT
		double FIS = (FIT - snp_Fst) / (1 - snp_Fst);
		cout << FIT << " " << FIS << endl;
		*/

		if ((S2 != 0) && (!isnan(S1)) && (!isnan(S2)) && (!isnan(snp_Fst)))
		{
			sum1 += S1;
			sum2 += S2;
			sum3 += snp_Fst;
			count++;
		}
		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << snp_Fst << endl;
	}

	double weighted_Fst = sum1 / sum2;
	double mean_Fst = sum3 / count;

	LOG.printLOG("Weir and Cockerham mean Fst estimate: " + output_log::dbl2str(mean_Fst, 5) + "\n");
	LOG.printLOG("Weir and Cockerham weighted Fst estimate: " + output_log::dbl2str(weighted_Fst, 5) + "\n");
	delete e;
}

void variant_file::output_windowed_weir_and_cockerham_fst(const string &output_file_prefix, const vector<string> &indv_files, int fst_window_size, int fst_window_step)
{
	if (fst_window_size <= 0)
		return;

	if ((fst_window_step <= 0) || (fst_window_step > fst_window_size))
		fst_window_step = fst_window_size;

	if (indv_files.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting Windowed Weir and Cockerham Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = indv_files.size();
	indvs_in_pops.resize(N_pops, vector<bool>(N_indv, false));
	vector<bool> all_indv(N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(indv_files[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + indv_files[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	// Find maximum position on each chromosome
	map<string, int> max_pos;
	map<string, int>::iterator it;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e->get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e->get_POS();
			}
			else
				max_pos[CHROM] = e->get_POS();
		}
	}

	// Calculate number of bins for each chromosome and allocate memory for them.
	// Each bin is a vector with four entries:
	// N_variant_sites: Number of sites in a window that have VCF entries
	// N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
	// N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
	// N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference allele
	unsigned int N_bins;
	const vector< double > empty_vector(4, 0);	// sum1, sum2, sum3, count
	map<string, vector< vector< double > > > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int) ceil( (max_pos[CHROM]+1) / double(fst_window_step));
		bins[CHROM].resize(N_bins, empty_vector);
	}

	double snp_Fst;
	double sum1=0.0, sum2 = 0.0;
	double sum3=0.0, count = 0.0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tFst: Only using biallelic sites.");
			continue;
		}

		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tFst: Only using diploid sites.");
			continue;
		}

		vector<unsigned int> n;
		n.resize(N_pops, 0);
		vector<double> p;
		p.resize(N_pops, 0);
		double nbar = 0.0, pbar=0.0, hbar=0.0;
		double ssqr=0.0;
		double sum_nsqr = 0.0;
		double n_sum = 0.0;
		unsigned int N_hom1, N_het, N_hom2;

		for (unsigned int j=0; j<N_pops; j++)
		{
			e->get_genotype_counts(indvs_in_pops[j], include_genotype[s], N_hom1, N_het, N_hom2);

			n[j] = N_hom1 + N_het + N_hom2;
			hbar += N_het;
			p[j] = N_het + 2*N_hom2;
			nbar += n[j];
			pbar += p[j];
			sum_nsqr += (n[j] * n[j]);

			p[j] /= (2.0*n[j]);	// diploid
		}

		n_sum = nbar;
		nbar /= double(N_pops);
		hbar /= n_sum;
		pbar /= (n_sum*2.0);	// diploid

		for (unsigned int j=0; j<N_pops; j++)
			ssqr += (n[j]*(p[j] - pbar)*(p[j] - pbar));
		ssqr /= ((N_pops-1.0)*nbar);

		double nc = (n_sum - (sum_nsqr / n_sum)) / (N_pops - 1.0);
		double C2 = N_pops * (1.0 - (nc / nbar));

		// diploid
		double S1 = ssqr - ((pbar * (1.0 - pbar) - ((N_pops-1.0) * ssqr / N_pops) - (hbar/4.0)) / (nbar - 1.0));
		double S2 = pbar * (1.0-pbar) * (1.0 - (nbar*C2/(N_pops*(nbar-1.0))));
		S2 += ((ssqr/N_pops)*(1.0 + ((N_pops-1.0)*nbar*C2/(N_pops*(nbar-1.0)))));
		S2 += (C2/(N_pops*(nbar-1.0)))*hbar/4.0;

		snp_Fst = S1 / S2;
		if ((S2 != 0) && (!isnan(S1)) && (!isnan(S2)) && (!isnan(snp_Fst)))
		{
			// Place the results into bins
			int pos = (int)e->get_POS();
			CHROM = e->get_CHROM();
			int first = (int) ceil((pos - fst_window_size)/double(fst_window_step));
			if (first < 0)
				first = 0;
			int last = (int) ceil(pos/double(fst_window_step));
			for(int idx = first; idx < last; idx++)
			{
				bins[CHROM][idx][0] += S1;
				bins[CHROM][idx][1] += S2;
				bins[CHROM][idx][2] += snp_Fst;
				bins[CHROM][idx][3]++;
			}
			sum1 += S1;
			sum2 += S2;
			sum3 += snp_Fst;
			count++;
		}
	}

	double weighted_Fst = sum1 / sum2;
	double mean_Fst = sum3 / count;

	LOG.printLOG("Weir and Cockerham mean Fst estimate: " + output_log::dbl2str(mean_Fst, 5) + "\n");
	LOG.printLOG("Weir and Cockerham weighted Fst estimate: " + output_log::dbl2str(weighted_Fst, 5) + "\n");

	string output = output_file_prefix + ".windowed.weir.fst";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Fst Output File: " + output, 7);

	out << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST" << endl;

	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			if ((bins[CHROM][s][1] != 0) && (!isnan(bins[CHROM][s][0])) && (!isnan(bins[CHROM][s][1])) && (bins[CHROM][s][3] > 0))
			{
				double weighted_Fst = bins[CHROM][s][0] / bins[CHROM][s][1];
				double mean_Fst = bins[CHROM][s][2] / bins[CHROM][s][3];

				out << CHROM << "\t"
				<< s*fst_window_step + 1 << "\t"
				<< (s*fst_window_step + fst_window_size) << "\t"
				<< bins[CHROM][s][3] << "\t"
				<< weighted_Fst << "\t" << mean_Fst << endl;
			}
		}
	}
	out.close();
	delete e;
}

void variant_file::output_windowed_hapmap_fst(const string &output_file_prefix, const vector<string> &indv_files, int fst_window_size, int fst_window_step)
{
	if (fst_window_size <= 0)
		return;

	if ((fst_window_step <= 0) || (fst_window_step > fst_window_size))
		fst_window_step = fst_window_size;

	if (indv_files.size() == 1)
	{
		LOG.printLOG("Require at least two populations to estimate Fst. Skipping\n");
		return;
	}

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Fst statistics.");

	LOG.printLOG("Outputting Windowed HapMap Fst estimates.\n");

	// First, read in the relevant files.
	vector< vector<bool> > indvs_in_pops;
	unsigned int N_pops = indv_files.size();
	indvs_in_pops.resize(N_pops, vector<bool>(N_indv, false));
	vector<bool> all_indv(N_indv,false);
	map<string, int> indv_to_idx;
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
			indv_to_idx[indv[ui]] = ui;
	for (unsigned int ui=0; ui<N_pops; ui++)
	{
		ifstream indv_file(indv_files[ui].c_str());
		if (!indv_file.is_open())
			LOG.error("Could not open Individual file: " + indv_files[ui]);
		string line;
		string tmp_indv;
		stringstream ss;
		while (!indv_file.eof())
		{
			getline(indv_file, line);
			ss.str(line);
			ss >> tmp_indv;
			if (indv_to_idx.find(tmp_indv) != indv_to_idx.end())
			{
				indvs_in_pops[ui][indv_to_idx[tmp_indv]]=true;
				all_indv[indv_to_idx[tmp_indv]]=true;
			}
			ss.clear();
		}
		indv_file.close();
	}

	// Find maximum position on each chromosome
	map<string, int> max_pos;
	map<string, int>::iterator it;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e->get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e->get_POS();
			}
			else
				max_pos[CHROM] = e->get_POS();
		}
	}

	// Calculate number of bins for each chromosome and allocate memory for them.
	// Each bin is a vector with four entries:
	// N_variant_sites: Number of sites in a window that have VCF entries
	// N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
	// N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
	// N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference allele
	unsigned int N_bins;
	const vector< double > empty_vector(4, 0);	// sum1, sum2, sum3, count
	map<string, vector< vector< double > > > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int) ceil( (max_pos[CHROM]+1) / double(fst_window_step));
		bins[CHROM].resize(N_bins, empty_vector);
	}

	vector<int> allele_counts1;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		if (e->get_N_alleles() != 2)
		{
			LOG.one_off_warning("\tFst: Only using biallelic sites.");
			continue;
		}

		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		unsigned int N_chr;
		e->get_allele_counts(allele_counts1, N_chr, all_indv, include_genotype[s]);
		double count_all = allele_counts1[1];
		double N_chr_all = N_chr;

		if ((count_all == 0) || (count_all == N_chr_all))
			continue;	// No polymorphism

		vector<double> counts(N_pops, 0);
		vector<double> pop_N_chr(N_pops, 0);
		vector<double> pop_N_choose_2(N_pops, 0);
		for (unsigned int p=0; p<N_pops; p++)
		{
			e->get_allele_counts(allele_counts1, N_chr, indvs_in_pops[p], include_genotype[s]);
			counts[p] = allele_counts1[1];
			pop_N_chr[p] = N_chr;
			pop_N_choose_2[p] = N_chr * (N_chr-1.0) / 2.0;
		}

		double Fst_SNP = 0;
		double f;
		double sum1=0.0;
		for (unsigned int p=0; p<N_pops; p++)
		{
			f = counts[p] / pop_N_chr[p];
			Fst_SNP += 2.0*pop_N_choose_2[p]*(pop_N_chr[p]/(pop_N_chr[p]-1.0))*f*(1.0-f);
			sum1 += pop_N_choose_2[p];
		}
		Fst_SNP /= sum1;
		double Fst_num = Fst_SNP;
		f = count_all / N_chr_all;
		double tmp = (2.0*(N_chr_all / (N_chr_all-1.0))*f*(1.0-f));
		Fst_SNP /= tmp;
		Fst_SNP = 1.0 - Fst_SNP;

		int pos = (int)e->get_POS();
		CHROM = e->get_CHROM();
		int first = (int) ceil((pos - fst_window_size)/double(fst_window_step));
		if (first < 0)
			first = 0;
		int last = (int) ceil(pos/double(fst_window_step));
		for(int idx = first; idx < last; idx++)
		{
			bins[CHROM][idx][0] += Fst_num;
			bins[CHROM][idx][1] += tmp;
			bins[CHROM][idx][2] += Fst_SNP;
			bins[CHROM][idx][3]++;
		}
	}

	string output = output_file_prefix + ".windowed.hapmap.fst";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Fst Output File: " + output, 7);

	out << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST" << endl;

	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			// if ((denominator_total != 0) && (!isnan(numerator)) && (!isnan(denominator_total)) && (!isnan(snp_Fst)))
			if ((bins[CHROM][s][1] != 0) && (!isnan(bins[CHROM][s][0])) && (!isnan(bins[CHROM][s][1])) && (bins[CHROM][s][3] > 0))
			{
				double weighted_Fst = 1.0 - (bins[CHROM][s][0] / bins[CHROM][s][1]);
				double mean_Fst = bins[CHROM][s][2] / bins[CHROM][s][3];

				out << CHROM << "\t"
				<< s*fst_window_step + 1 << "\t"
				<< (s*fst_window_step + fst_window_size) << "\t"
				<< bins[CHROM][s][3] << "\t"
				<< weighted_Fst << "\t" << mean_Fst << endl;
			}
		}
	}
	out.close();
	delete e;
}

void variant_file::output_per_site_nucleotide_diversity(const string &output_file_prefix)
{
	// Output nucleotide diversity, calculated on a per-site basis.
	// Pi = average number of pairwise differences
	// Assumes a constant distance of 1 between all possible mutations
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Nucleotide Diversity Statistics.");

	LOG.printLOG("Outputting Per-Site Nucleotide Diversity Statistics...\n");
	string output_file = output_file_prefix + ".sites.pi";

	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "CHROM\tPOS\tPI" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	vector<int> allele_counts;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		e->parse_full_entry(true);
		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tsitePi: Only using fully diploid sites.");
			continue;
		}

		unsigned int N_non_missing_chr;
		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		unsigned int total_alleles = std::accumulate(allele_counts.begin(), allele_counts.end(), 0);

		unsigned int N_alleles = e->get_N_alleles();
		int mismatches = 0;
		for(unsigned int allele = 0; allele < N_alleles; allele++)
		{
			int other_alleles_count = (total_alleles - allele_counts[allele]);
			mismatches += (allele_counts[allele] * other_alleles_count);
		}

		int pairs = (total_alleles * (total_alleles - 1));
		double pi = (mismatches/static_cast<double>(pairs));

		out << e->get_CHROM() << "\t" << e->get_POS() << "\t" << pi << endl;
	}
	delete e;
}

// Output Tajima's D
// Carlson et al. Genome Res (2005)
void variant_file::output_Tajima_D(const string &output_file_prefix, int window_size)
{
	if (window_size <= 0)
		return;

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Tajima's D Statistic.");

	LOG.printLOG("Outputting Tajima's D Statistic...\n");
	string output_file = output_file_prefix + ".Tajima.D";

	double a1=0.0, a2=0.0, b1, b2, c1, c2, e1, e2;
	unsigned int n = N_kept_individuals()*2;
	if (n < 2)
		LOG.error("Require at least two chromosomes!");

	for (unsigned int ui=1; ui<n; ui++)
	{
		a1 += 1.0 / double(ui);
		a2 += 1.0 / double(ui * ui);
	}
	b1 = double(n+1) / 3.0 / double(n-1);
	b2 = 2.0 * double(n*n + n + 3) / 9.0 / double(n) / double(n-1);
	c1 = b1 - (1.0 / a1);
	c2 = b2 - (double(n+2)/double(a1*n)) + (a2/a1/a1);
	e1 = c1 / a1;
	e2 = c2 / ((a1*a1) + a2);

	// Find maximum position
	map<string, int> max_pos;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e->get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e->get_POS();
			}
			else
				max_pos[CHROM] = e->get_POS();
		}
	}

	map<string, int>::iterator it;
	unsigned int N_bins;
	map<string, vector< pair<int, double> > > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int)((max_pos[CHROM] + window_size) / double(window_size));
		bins[CHROM].resize(N_bins, make_pair(0,0));
	}

	unsigned int idx;
	double C = 1.0 / double(window_size);
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned int N_alleles;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tTajimaD: Only using bialleleic sites.");
			continue;
		}

		CHROM = e->get_CHROM();
		idx = (unsigned int)(e->get_POS() * C);
		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tTajimaD: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);

		double p = double(allele_counts[0]) / N_non_missing_chr;
		if ((p > 0.0) && (p < 1.0))
		{
			bins[CHROM][idx].first++;
			bins[CHROM][idx].second += p * (1.0-p);
		}
	}

	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "CHROM\tBIN_START\tN_SNPS\tTajimaD" << endl;

	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		bool output = false;
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			int S = bins[CHROM][s].first;
			double D = 0;
			if (S > 1)
			{
				double pi = 2.0*bins[CHROM][s].second*n/double(n-1);
				double tw = double(S) / a1;
				double var = (e1*S) + e2*S*(S-1);
				D = (pi - tw) / sqrt(var);
				output = true;
			}
			if (S > 0)
				output = true;
			if (output == true)
				out << CHROM << "\t" << s*window_size << "\t" << bins[CHROM][s].first << "\t" << D << endl;
		}
	}
	delete e;
	out.close();
}

void variant_file::output_windowed_nucleotide_diversity(const string &output_file_prefix, int window_size, int window_step)
{
	// Output nucleotide diversity, as calculated in windows.
	// Average number of pairwise differences in windows.
	if (window_size <= 0)
		return;

	if ((window_step <= 0) || (window_step > window_size))
		window_step = window_size;

	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Nucleotide Diversity Statistics.");

	LOG.printLOG("Outputting Windowed Nucleotide Diversity Statistics...\n");
	string output_file = output_file_prefix + ".windowed.pi";

	// Find maximum position on each chromosome
	map<string, int> max_pos;
	map<string, int>::iterator it;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();

			if (max_pos.find(CHROM) != max_pos.end())
			{
				if (e->get_POS() > max_pos[CHROM])
					max_pos[CHROM] = e->get_POS();
			}
			else
				max_pos[CHROM] = e->get_POS();
		}
	}

	// Calculate number of bins for each chromosome and allocate memory for them.
	// Each bin is a vector with four entries:
	// N_variant_sites: Number of sites in a window that have VCF entries
	// N_variant_site_pairs: Number of possible pairwise mismatches at polymorphic sites within a window
	// N_mismatches: Number of actual pairwise mismatches at polymorphic sites within a window
	// N_polymorphic_sites: number of sites within a window where there is at least 1 sample that is polymorphic with respect to the reference allele
	unsigned int N_bins;
	const unsigned int N_variant_sites = 0;
	const unsigned int N_variant_site_pairs = 1;
	const unsigned int N_mismatches = 2;
	const unsigned int N_polymorphic_sites = 3;
	const vector< unsigned long > empty_vector(4, 0);
	map<string, vector< vector< unsigned long> > > bins;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		N_bins = (unsigned int) ceil( (max_pos[CHROM]+1) / double(window_step));
		bins[CHROM].resize(N_bins, empty_vector);
	}

	// Count polymorphic sites and pairwise mismatches
	vector<int> allele_counts;
	unsigned int N_non_missing_chr;
	unsigned long N_comparisons;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		CHROM = e->get_CHROM();

		e->parse_genotype_entries(true);

		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\twindowPi: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);

		unsigned int N_site_mismatches = 0;
		for (vector<int>::iterator ac = allele_counts.begin(); ac != allele_counts.end(); ++ac)
		{
			N_site_mismatches += (*ac * (N_non_missing_chr - *ac));
		}

		if (N_site_mismatches == 0)
				continue;	// Site is actually fixed.

		// Place the counts into bins
		int pos = (int)e->get_POS();
		int first = (int) ceil((pos - window_size)/double(window_step));
		if (first < 0)
			first = 0;
		int last = (int) ceil(pos/double(window_step));
		N_comparisons = N_non_missing_chr * (N_non_missing_chr - 1);
		for(int idx = first; idx < last; idx++)
		{
			bins[CHROM][idx][N_variant_sites]++;
			bins[CHROM][idx][N_variant_site_pairs] += N_comparisons;
			bins[CHROM][idx][N_mismatches] += N_site_mismatches;

			if(allele_counts[0] < (signed)N_non_missing_chr)
				bins[CHROM][idx][N_polymorphic_sites]++;
		}
	}

	// Calculate and print nucleotide diversity statistics
	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI" << endl;

	unsigned long N_monomorphic_sites = 0;
	int N_kept_chr = 2*N_kept_individuals();
	N_comparisons = (N_kept_chr * (N_kept_chr - 1)); 	// Number of pairwise comparisons at a monomorphic site
	unsigned long N_pairs = 0; 								// Number of pairwise comparisons within a window
	double pi = 0;
	for (it=max_pos.begin(); it != max_pos.end(); ++it)
	{
		CHROM = (*it).first;
		for (unsigned int s=0; s<bins[CHROM].size(); s++)
		{
			if( (bins[CHROM][s][N_polymorphic_sites] > 0) || (bins[CHROM][s][N_mismatches] > 0) )
			{
				// This number can be slightly off for the last bin since the
				// window size can go off the end of the chromosome.
				N_monomorphic_sites = window_size - bins[CHROM][s][N_variant_sites];

				// The total number of possible pairwise comparisons is the sum of
				// pairwise comparisons at polymorphic sites and pairwise
				// comparisons at monomorphic sites.
				N_pairs = bins[CHROM][s][N_variant_site_pairs] + (N_monomorphic_sites * N_comparisons);

				pi = bins[CHROM][s][N_mismatches] / double(N_pairs);
				out << CHROM << "\t"
				    << s*window_step + 1 << "\t"
				    << (s*window_step + window_size) << "\t"
				    << bins[CHROM][s][N_polymorphic_sites] << "\t"
				    << pi << endl;
			}
		}
	}
	delete e;
	out.close();
}

void variant_file::output_kept_and_removed_sites(const string &output_file_prefix)
{
	// Output lists of sites that have been filtered (or not).
	LOG.printLOG("Outputting Kept and Removed Sites...\n");
	string output_file1 = output_file_prefix + ".kept.sites";
	string output_file2 = output_file_prefix + ".removed.sites";

	string CHROM;
	vector<char> variant_line;
	int POS;
	entry *e = get_entry_object(N_indv);

	ofstream out1(output_file1.c_str());
	if (!out1.is_open()) LOG.error("Could not open output file: " + output_file1, 12);
	out1 << "CHROM\tPOS" << endl;

	ofstream out2(output_file2.c_str());
	if (!out2.is_open()) LOG.error("Could not open output file: " + output_file2, 12);
	out2 << "CHROM\tPOS" << endl;

	for (unsigned int s=0; s<N_entries; s++)
	{
		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry();
		POS = e->get_POS();
		CHROM = e->get_CHROM();
		if (include_entry[s] == true)
			out1 << CHROM << "\t" << POS << endl;
		else
			out2 << CHROM << "\t" << POS << endl;
	}
	delete e;
	out1.close();
	out2.close();
}

void variant_file::output_LROH(const string &output_file_prefix)
{
	// Detect and output Long Runs of Homozygosity, following the method
	// developed by Adam Boyko, and described in Auton et al., Genome Research, 2009
	// (Although using Forward-backwards algorithm in place of Viterbi).
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output LROH.");

	LOG.printLOG("Outputting Long Runs of Homozygosity (Experimental)... \n");
	string output_file = output_file_prefix + ".LROH";

	unsigned int nGen=4;				// Number of generations since common ancestry
	double genotype_error_rate = 0.01;	// Assumed genotype error rate
	double p_auto_prior = 0.05;			// Prior probability of being in autozygous state
	double p_auto_threshold = 0.99;		// Threshold for reporting autozygous region
	int min_SNPs=0;						// Threshold for reporting autozygous region

	string CHROM;
	vector<char> variant_line;
	int POS;
	entry *e = get_entry_object(N_indv);
	pair<int, int> alleles;
	vector<unsigned int> s_vector;
	vector<pair<double, double> > p_emission;
	vector<vector<double> > p_trans;

	ofstream out(output_file.c_str());
	if (!out.is_open()) LOG.error("Could not open output file: " + output_file, 12);
	out << "CHROM\tAUTO_START\tAUTO_END\tN_VARIANTS\tINDV" << endl;

	// TODO - refactor this so that Entries loop is on the outside.
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;

		LOG.printLOG("\t" + indv[ui] + "\n");

		int last_POS = -1;
		s_vector.resize(0);	p_emission.resize(0); p_trans.resize(0);

		for (unsigned int s=0; s<N_entries; s++)
		{
			if ((include_entry[s] == false) || (include_genotype[s][ui] == false))
				continue;

			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true);

			if (e->get_N_alleles() != 2)
			{
				LOG.one_off_warning("\tLROH: Only using bialleleic sites.");
				continue;	// TODO: Probably could do without this...
			}

			POS = e->get_POS();

			e->parse_genotype_entry(ui, true);
			e->get_indv_GENOTYPE_ids(ui, alleles);

			if (e->get_indv_ploidy(ui) != 2)
			{
				LOG.one_off_warning("\tLROH: Only using diploid sites.");
				continue;
			}

			if ((alleles.first == -1) || (alleles.second == -1))
				continue;

			unsigned int X = alleles.first + alleles.second;

			// Calculate heterozyogosity of this site.
			// TODO: Would be better to do this once, but for simplicity, do it for each individual.
			unsigned int N_genotypes = 0;
			unsigned int N_hets = 0;
			for (unsigned int uj=0; uj<N_indv; uj++)
			{
				if ((include_indv[uj] == false) || (include_genotype[s][ui] == false))
					continue;

				e->parse_genotype_entry(uj, true);
				e->get_indv_GENOTYPE_ids(uj, alleles);
				if ((alleles.first != -1) && (alleles.second != -1))
				{
					N_genotypes++;
					if (alleles.first != alleles.second)
						N_hets++;
				}
			}
			double h = N_hets / double(N_genotypes);
			double p_emission_given_nonauto;
			double p_emission_given_auto;

			if (X == 1)
			{	// Heterozygote
				p_emission_given_nonauto = h;
				p_emission_given_auto = genotype_error_rate;
				p_emission.push_back(make_pair(p_emission_given_auto, p_emission_given_nonauto));
			}
			else
			{	// Homozygote
				p_emission_given_nonauto = 1.0-h;
				p_emission_given_auto = 1.0-genotype_error_rate;
				p_emission.push_back(make_pair(p_emission_given_auto, p_emission_given_nonauto));
			}

			double r = 0;
			if (last_POS > 0)
			{	// Assume 1cM/Mb.
				r = (POS - last_POS) / 1000000.0 / 100.0; // Morgans
			}

			double e = (1.0 - exp(-2.0*nGen*r));
			double p_trans_auto_to_nonauto = (1.0 - p_auto_prior) * e;	//A[1]
			double p_trans_nonauto_to_auto = p_auto_prior * e;	//A[2]
			double p_trans_auto_to_auto = 1.0 - p_trans_nonauto_to_auto;	//A[0]
			double p_trans_nonauto_to_nonauto = 1.0 - p_trans_auto_to_nonauto; // A[3]
			vector<double> A(4);
			A[0] = p_trans_auto_to_auto;
			A[1] = p_trans_auto_to_nonauto;
			A[2] = p_trans_nonauto_to_auto;
			A[3] = p_trans_nonauto_to_nonauto;

			s_vector.push_back(s);

			p_trans.push_back(A);
			last_POS = POS;
		}

		// Forward-backward algorithm
		int N_obs = (int)p_emission.size();
		if (N_obs == 0)
			continue;

		vector<vector<double> > alpha(N_obs, vector<double>(2,0));
		vector<vector<double> > beta(N_obs, vector<double>(2,0));

		alpha[0][0] = p_emission[0].first;
		alpha[0][1] = p_emission[0].second;
		for (int i=1; i<N_obs; i++)
		{
			alpha[i][0] = alpha[i-1][0] * p_trans[i-1][0] * p_emission[i].first;
			alpha[i][0] += alpha[i-1][1] * p_trans[i-1][2] * p_emission[i].first;

			alpha[i][1] = alpha[i-1][1] * p_trans[i-1][3] * p_emission[i].second;
			alpha[i][1] += alpha[i-1][0] * p_trans[i-1][1] * p_emission[i].second;

			while (alpha[i][0] + alpha[i][1] < 1e-20)
			{	// Renormalise to prevent underflow
				alpha[i][0] *= 1e20;
				alpha[i][1] *= 1e20;
			}
		}

		beta[N_obs-1][0] = 1.0;
		beta[N_obs-1][1] = 1.0;
		for (int i=N_obs-2; i>=0; i--)
		{
			beta[i][0] = beta[i+1][0] * p_trans[i][0] * p_emission[i].first;
			beta[i][0] += beta[i+1][1] * p_trans[i][2] * p_emission[i].first;

			beta[i][1] = beta[i+1][1] * p_trans[i][3] * p_emission[i].second;
			beta[i][1] += beta[i+1][0] * p_trans[i][1] * p_emission[i].second;

			while (beta[i][0] + beta[i][1] < 1e-20)
			{	// Renormalise to prevent underflow
				beta[i][0] *= 1e20;
				beta[i][1] *= 1e20;
			}
		}

		// Calculate probability of each site being autozygous
		vector<double> p_auto(N_obs);
		for (int i=0; i<N_obs; i++)
			p_auto[i] = alpha[i][0] * beta[i][0] / (alpha[i][0] * beta[i][0] + alpha[i][1] * beta[i][1]);

		// Generate output
		// TODO: Would be good to report actual limits of homozygosity
		// (i.e. extend regions out until first heterozygote),
		// as opposed to regions with p>threshold.
		// TODO: Also would be good to report heterozygotic SNPs found in homozygotic regions.
		bool in_auto=false;
		int start_pos=0, end_pos=0;
		int N_SNPs = 0;
		for (int i=0; i<N_obs; i++)
		{
			if (p_auto[i] > p_auto_threshold)
			{
				if (in_auto == false)
				{	// Start of autozygous region
					unsigned int s = s_vector[i];
					get_entry(s, variant_line);
					e->reset(variant_line);
					e->parse_basic_entry(true);
					CHROM = e->get_CHROM();
					start_pos = e->get_POS();
				}
				N_SNPs++;
				in_auto = true;
			}
			else
			{
				if (in_auto == true)
				{	// end of autozygous region
					unsigned int s = s_vector[i];
					get_entry(s, variant_line);
					e->reset(variant_line);
					e->parse_basic_entry(true);
					end_pos = e->get_POS();
					if (N_SNPs >= min_SNPs)
						out << CHROM << "\t" << start_pos << "\t" << end_pos << "\t" << N_SNPs << "\t" << indv[ui] << endl;
				}
				in_auto = false;
				N_SNPs = 0;
			}
		}
		if (in_auto == true)
		{	// Report final region if needed
			unsigned int s = s_vector[N_obs-1];
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true);
			end_pos = e->get_POS();
			if (N_SNPs >= min_SNPs)
				out << CHROM << "\t" << start_pos << "\t" << end_pos << "\t" << N_SNPs << "\t" << indv[ui] << endl;
		}
	}
	delete e;
	out.close();
}

void variant_file::output_indv_relatedness(const string &output_file_prefix)
{
	// Calculate and output a relatedness statistic based on the method of
	// Yang et al, 2010 (doi:10.1038/ng.608). Specifically, calculate the
	// unadjusted Ajk statistic (equation 6 of paper).
	// Expectation of Ajk is zero for individuals within a populations, and
	// one for an individual with themselves.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to output Individual Relatedness.");

	LOG.printLOG("Outputting Individual Relatedness\n");
	string output = output_file_prefix + ".relatedness";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Individual Relatedness Output File: " + output, 2);
	out << "INDV1\tINDV2\tRELATEDNESS" << endl;

	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	vector<int> allele_counts;
	unsigned int N_alleles, N_non_missing_chr;
	double freq;
	pair<int, int> geno_id;
	vector<vector<double> > Ajk(N_indv, vector<double>(N_indv, 0.0));
	vector<vector<double> > N_sites(N_indv, vector<double>(N_indv, 0.0));

	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);

		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();

		if (N_alleles != 2)
		{
			LOG.one_off_warning("\tRelatedness: Only using biallelic sites.");
			continue;	// Only use biallelic loci
		}

		e->parse_genotype_entries(true);
		if (e->is_diploid(include_indv, include_genotype[s]) == false)
		{
			LOG.one_off_warning("\tRelatedness: Only using fully diploid sites.");
			continue;
		}

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

		if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
			continue;

		vector<double> x(N_indv, -1.0);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			x[ui] = geno_id.first + geno_id.second;
		}

		double div = 1.0/(2.0*freq*(1.0-freq));
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if ((include_indv[ui] == false) || (include_genotype[s][ui] == false) || (x[ui] < 0))
				continue;
			Ajk[ui][ui] += (x[ui]*x[ui] - (1 + 2.0*freq)*x[ui] + 2.0*freq*freq) * div;
			N_sites[ui][ui]++;
			for (unsigned int uj=(ui+1); uj<N_indv; uj++)
			{
				if ((include_indv[uj] == false) || (include_genotype[s][uj] == false) || (x[uj] < 0))
					continue;
				Ajk[ui][uj] += (x[ui] - 2.0*freq) * (x[uj] - 2.0*freq) * div;
				N_sites[ui][uj]++;
			}
		}
	}

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		Ajk[ui][ui] = 1.0 + (Ajk[ui][ui] / N_sites[ui][ui]);
		out << indv[ui] << "\t" << indv[ui] << "\t" << Ajk[ui][ui] << endl;
		for (unsigned int uj=(ui+1); uj<N_indv; uj++)
		{
			if (include_indv[uj] == false)
				continue;
			Ajk[ui][uj] /= N_sites[ui][uj];
			out << indv[ui] << "\t" << indv[uj] << "\t" << Ajk[ui][uj] << endl;
		}
	}
	delete e;
	out.close();
}

void variant_file::output_PCA(const string &output_file_prefix, bool use_normalisation, int SNP_loadings_N_PCs)
{
#ifndef VCFTOOLS_PCA
	use_normalisation = true;
	SNP_loadings_N_PCs = -1;
	string out = "Cannot run PCA analysis. Vcftools has been compiled without PCA enabled (requires LAPACK).";
	LOG.error(out);
#else
	// Output PCA, following method of Patterson, Price and Reich 2006.
	if ((has_genotypes == false) | (N_kept_individuals() == 0))
		LOG.error("Require Genotypes in VCF file in order to perform PCA.");

	if (use_normalisation)
		LOG.printLOG("Outputting Principal Component Analysis (with normalisation)\n");
	else
		LOG.printLOG("Outputting Principal Component Analysis (without normalisation)\n");
	string output = output_file_prefix + ".pca";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Principal Component Analysis Output File: " + output, 2);

	unsigned int N_indvs = N_kept_individuals();
	unsigned int N_sites = N_kept_sites();

	if (N_indvs >= N_sites)
		LOG.error("PCA computation requires that there are more sites than individuals.");

	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	pair<int, int> geno_id;
	double x, freq;
	vector<int> allele_counts;
	unsigned int N_alleles, N_non_missing_chr;

	// Store list of included individuals
	vector<string> included_indvs(N_indvs);
	unsigned int ui_prime = 0;
	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if (include_indv[ui] == false)
			continue;
		included_indvs[ui_prime] = indv[ui];
		ui_prime++;
	}

	// Potentially uses a lot of memory. Should issue a warning about this.
	double **M = new double*[N_indvs];	// m rows = indv
	for (unsigned int ui=0; ui<N_indvs; ui++)
		M[ui] = new double[N_sites];	// n columns

	// Populate M
	unsigned int s_prime = 0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s]==false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);

		e->parse_basic_entry(true);
		N_alleles = e->get_N_alleles();
		if (N_alleles != 2)
			LOG.error("PCA only works for biallelic sites.");

		e->parse_genotype_entries(true);
		if (e->is_diploid(include_indv, include_genotype[s]) == false)
			LOG.error("PCA only works for fully diploid sites. Non-diploid site at " + e->get_CHROM() + ":" + output_log::int2str(e->get_POS()));

		e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
		freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

		if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
			continue;

		double mu = freq*2.0;
		double div = 1.0 / sqrt(freq * (1.0-freq));

		ui_prime = 0;
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;

			e->get_indv_GENOTYPE_ids(ui, geno_id);
			x = geno_id.first + geno_id.second;
			if (x > -1)
			{
				if (use_normalisation == true)
					M[ui_prime][s_prime] = (x - mu) * div;
				else
					M[ui_prime][s_prime] = (x - mu);
			}
			ui_prime++;
		}
		s_prime++;
	}

	// Now construct X = (1/n)MM'.
	double **X = new double *[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		X[ui] = new double[N_indvs];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] = 0;

	// Only populate one half of matrix
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=ui; uj<N_indvs; uj++)
			for (unsigned int s=0; s<N_sites; s++)
				X[ui][uj] += M[ui][s] * M[uj][s];

	delete [] M;

	// Populate other half
	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<ui; uj++)
			X[ui][uj] = X[uj][ui];

	for (unsigned int ui=0; ui<N_indvs; ui++)
		for (unsigned int uj=0; uj<N_indvs; uj++)
			X[ui][uj] /= N_sites;

	double *Er = new double[N_indvs];
	double *Ei = new double[N_indvs];
	double **Evecs = new double*[N_indvs];
	for (unsigned int ui=0; ui<N_indvs; ui++)
		Evecs[ui] = new double[N_indvs];

	// Call LAPACK routine to calculate eigenvectors and eigenvalues
	dgeev(X, N_indvs, Er, Ei, Evecs);

	// Check there are no complex eigenvalues.
	for (unsigned int ui=0; ui<N_indvs; ui++)
		if (Ei[ui] != 0)
			LOG.error("Complex eigenvalue.");

	// Output results
	out << "INDV";
	for (unsigned int ui=0; ui<N_indvs; ui++)
		out << "\tEIG_" << ui;
	out << endl;

	out << "EIGENVALUE";
	for (unsigned int ui=0; ui<N_indvs; ui++)
		out << "\t" << Er[ui];
	out << endl;

	// Output eigenvectors (as columns)
	for (unsigned int ui=0; ui<N_indvs; ui++)
	{
		out << included_indvs[ui];
		for (unsigned int uj=0; uj<N_indvs; uj++)
			out << "\t" << Evecs[ui][uj];
		out << endl;
	}

	out.close();

	if (SNP_loadings_N_PCs > 0)
	{	// Output SNP loadings
		LOG.printLOG("Outputting " + output_log::int2str(SNP_loadings_N_PCs) + " SNP loadings\n");
		output = output_file_prefix + ".pca.loadings";
		out.open(output.c_str());
		if (!out.good())
			LOG.error("Could not open Principal Component SNP Loading Output File: " + output, 2);
		out << "CHROM\tPOS";
		for (unsigned int ui=0; ui<(unsigned int)SNP_loadings_N_PCs; ui++)
			out << "\tGAMMA_" << ui;
		out << endl;

		for (unsigned int s=0; s<N_entries; s++)
		{
			if (include_entry[s]==false)
				continue;

			get_entry(s, variant_line);
			e->reset(variant_line);

			e->parse_basic_entry(true);
			N_alleles = e->get_N_alleles();
			if (N_alleles != 2)
				LOG.error("PCA only works for biallelic sites.");

			e->parse_genotype_entries(true);
			if (e->is_diploid(include_indv, include_genotype[s]) == false)
				LOG.error("PCA only works for fully diploid sites.");

			e->get_allele_counts(allele_counts, N_non_missing_chr, include_indv, include_genotype[s]);
			freq = allele_counts[1] / (double)N_non_missing_chr;	// Alt allele frequency

			if ((freq <= numeric_limits<double>::epsilon()) || (freq >= (1.0-numeric_limits<double>::epsilon())))
				continue;

			vector<double> gamma(SNP_loadings_N_PCs, 0.0);
			vector<double> a_sum(SNP_loadings_N_PCs, 0.0);

			ui_prime = 0;
			for (unsigned int ui=0; ui<N_indv; ui++)
			{
				if (include_indv[ui] == false)
					continue;

				e->get_indv_GENOTYPE_ids(ui, geno_id);
				x = geno_id.first + geno_id.second;
				if (x > -1)
				{
					for (unsigned int uj=0; uj<(unsigned int)SNP_loadings_N_PCs; uj++)
					{
						gamma[uj] += (x * Evecs[ui_prime][uj]);
						a_sum[uj] += (Evecs[ui_prime][uj]*Evecs[ui_prime][uj]);
					}
				}
				ui_prime++;
			}

			out << e->get_CHROM() << "\t" << e->get_POS();
			for (unsigned int uj=0; uj<(unsigned int)SNP_loadings_N_PCs; uj++)
				out << "\t" << gamma[uj] / a_sum[uj];
			out << endl;
		}
		out.close();
	}
	delete e;
	delete [] Er;
	delete [] Ei;
	delete [] Evecs;
	delete [] X;
#endif
}

void variant_file::output_indel_hist(const string &output_file_prefix)
{
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	string allele;
	unsigned int ref_len, N_alleles;
	int indel_len, smallest_len, largest_len, snp_count;
	vector<int> s_vector;

	string output = output_file_prefix + ".indel.hist";
	ofstream out(output.c_str());
	if (!out.is_open())
		LOG.error("Could not open Indel Hist File: " + output, 7);

	LOG.printLOG("Outputting Indel Histogram\n");
	out << "LENGTH\tCOUNT\tPRCT" << endl;
	largest_len = 0;
	smallest_len = 0;
	snp_count = 0;
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == false)
			continue;

		get_entry(s, variant_line);
		e->reset(variant_line);
		e->parse_basic_entry(true);

		allele = e->get_REF();
		ref_len = allele.size();
		N_alleles = e->get_N_alleles();
		if (e->is_SNP() )
			snp_count++;

		for (unsigned int ui=1; ui<N_alleles; ui++)
		{
			e->get_allele(ui, allele);
			if (allele.size() != ref_len)
			{
				indel_len = allele.size() - ref_len;
				s_vector.push_back (indel_len);
				if (indel_len > largest_len)
					largest_len = indel_len;
				else if (indel_len < smallest_len)
					smallest_len = indel_len;

			}
		}
	}
	double total = s_vector.size() + snp_count;
	double pct;
	for (int i=smallest_len; i<=largest_len; i++)
	{
		int icount = (int) count (s_vector.begin(), s_vector.end(), i);
		if (icount > 0)
		{
			pct = 100.0*icount/total;
			out << i << "\t" << icount << "\t" << pct << endl;

		}
		else if ((i == 0) and (snp_count>0))
		{
			pct = 100.0*snp_count/total;
			out << i << "\t" << snp_count << "\t" << pct << endl;
		}
	}
	out.close();
}

