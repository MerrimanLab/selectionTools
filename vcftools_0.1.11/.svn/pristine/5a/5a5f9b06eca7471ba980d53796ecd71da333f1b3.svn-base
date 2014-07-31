/*
 * variant_file_diff.cpp
 *
 *  Created on: Oct 30, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "variant_file.h"

void variant_file::return_site_union(variant_file &file2, map<pair<string, int>, pair<int, int> > &CHROMPOS_to_filepos_pair)
{
	unsigned int s;
	int POS;
	string CHROM;
	vector<char> variant_line;
	entry *e = get_entry_object(N_indv);
	entry *e2 = file2.get_entry_object(file2.N_indv);
	for (s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry();

			CHROM = e->get_CHROM();
			POS = e->get_POS();

			CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)] = make_pair<int,int>(s, -1);
		}
	}
	for (s=0; s<file2.N_entries; s++)
	{
		if (file2.include_entry[s] == true)
		{
			file2.get_entry(s, variant_line);
			e2->reset(variant_line);
			e2->parse_basic_entry();

			CHROM = e2->get_CHROM();
			POS = e2->get_POS();

			if (CHROMPOS_to_filepos_pair.find(make_pair<string,int>(CHROM, POS)) != CHROMPOS_to_filepos_pair.end())
			{
				CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)].second = s;
			}
			else
			{
				CHROMPOS_to_filepos_pair[make_pair<string,int>(CHROM, POS)] = make_pair<int,int>(-1, s);
			}
		}
	}
	delete e;
	delete e2;
}

void variant_file::return_indv_union(variant_file &file2, map<string, pair< int, int> > &combined_individuals, const string &indv_ID_map_file)
{
	map<string, string> indv_map;
	bool use_map = false;
	if (indv_ID_map_file != "")
	{
		LOG.printLOG("Reading individual mapping file. ");
		ifstream map(indv_ID_map_file.c_str());
		if (!map.is_open())
			LOG.error("Could not open map file: " + indv_ID_map_file);
		while (!map.eof())
		{
			string indv1, indv2;
			map >> indv1 >> indv2;
			map.ignore(numeric_limits<streamsize>::max(), '\n');
			if ((indv1 != "") && (indv1.substr(0,1) != "#"))
			{
				indv_map[indv1] = indv2;
			}
		}
		map.close();
		use_map = true;
		LOG.printLOG("Read " + LOG.int2str(indv_map.size()) + " entries.\n");
	}

	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui] == true)
		{
			combined_individuals[indv[ui]] = make_pair<int,int>(ui, -1);
		}

	for (unsigned int ui=0; ui<file2.N_indv; ui++)
		if (file2.include_indv[ui] == true)
		{
			string indv_id = file2.indv[ui];
			if (use_map == true)
			{
				if (indv_map.find(indv_id) != indv_map.end())
					indv_id = indv_map[indv_id];
			}
			if (combined_individuals.find(indv_id) != combined_individuals.end())
			{
				combined_individuals[indv_id].second = ui;
				//cout << combined_individuals[indv_id].first << " " << combined_individuals[indv_id].second << endl;
			}
			else
				combined_individuals[indv_id] = make_pair<int,int>(-1, ui);
		}
}

void variant_file::output_sites_in_files(const string &output_file_prefix, variant_file &diff_variant_file)
{
	LOG.printLOG("Comparing sites in VCF files...\n");

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_variant_file, CHROMPOS_to_filepos_pair);

	vector<char> variant_line;
	string CHROM;
	int POS;

	string output_file = output_file_prefix + ".diff.sites_in_files";
	ofstream sites_in_files(output_file.c_str());
	sites_in_files << "CHROM\tPOS\tIN_FILE\tREF\tALT1\tALT2" << endl;

	int s1, s2;
	int N_common_SNPs = 0, N_SNPs_file1_only=0, N_SNPs_file2_only=0;
	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it!=CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		entry *e1 = get_entry_object(N_indv);
		entry *e2 = diff_variant_file.get_entry_object(diff_variant_file.N_indv);

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_entry(s1, variant_line);
			e1->reset(variant_line);
		}

		if (s2 != -1)
		{
			diff_variant_file.get_entry(s2, variant_line);
			e2->reset(variant_line);
		}

		e1->parse_basic_entry(true);
		e2->parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1->get_REF();
		string REF2 = e2->get_REF();

		if ((REF == "N") || (REF == ".") || (REF == "") )
			REF = REF2;
		if ((REF2 == "N") || (REF2 == ".") || (REF2 == "") )
			REF2 = REF;

		if ((REF != REF2) && (REF2 != "N") && (REF != "N") && (REF != ".") && (REF2 != ".") && (REF != "") && (REF2 != ""))
		{
			LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
			continue;
		}

		sites_in_files << CHROM << "\t" << POS << "\t";
		if ((s1 != -1) && (s2 != -1))
		{
			N_common_SNPs++;
			sites_in_files << "B";
		}
		else if ((s1 != -1) && (s2 == -1))
		{
			N_SNPs_file1_only++;
			sites_in_files << "1";
		}
		else if ((s1 == -1) && (s2 != -1))
		{
			N_SNPs_file2_only++;
			sites_in_files << "2";
		}
		else
			LOG.error("SNP in neither file!?");

		sites_in_files << "\t" << REF << "\t" << e1->get_ALT() << "\t" << e2->get_ALT() << endl;

		delete e1;
		delete e2;
	}

	sites_in_files.close();

	LOG.printLOG("Found " + output_log::int2str(N_common_SNPs) + " SNPs common to both files.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file1_only) + " SNPs only in main file.\n");
	LOG.printLOG("Found " + output_log::int2str(N_SNPs_file2_only) + " SNPs only in second file.\n");
}

void variant_file::output_indv_in_files(const string &output_file_prefix, variant_file &diff_variant_file, const string &indv_ID_map_file)
{
	LOG.printLOG("Comparing individuals in VCF files...\n");

	string output_file = output_file_prefix + ".diff.indv_in_files";

	ofstream out(output_file.c_str());
	if (!out.is_open())
		LOG.error("Could not open Indv Differences File: " + output_file, 3);
	out << "INDV\tFILES" << endl;

	// Build a list of individuals contained in each file
	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, indv_ID_map_file);

	unsigned int N_combined_indv = combined_individuals.size();
	unsigned int N[3]={0,0,0};
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		if ((combined_individuals_it->second.first != -1) && (combined_individuals_it->second.second != -1))
		{
			N[0]++;
			out << combined_individuals_it->first << "\tB" << endl;
		}
		else if (combined_individuals_it->second.first != -1)
		{
			N[1]++;
			out << combined_individuals_it->first << "\t1" << endl;
		}
		else if (combined_individuals_it->second.second != -1)
		{
			N[2]++;
			out << combined_individuals_it->first << "\t2" << endl;
		}
		else
			LOG.error("Unhandled case");
	}
	out.close();

	LOG.printLOG("N_combined_individuals:\t" + output_log::int2str(N_combined_indv) + "\n");
	LOG.printLOG("N_individuals_common_to_both_files:\t" + output_log::int2str(N[0]) + "\n");
	LOG.printLOG("N_individuals_unique_to_file1:\t" + output_log::int2str(N[1]) + "\n");
	LOG.printLOG("N_individuals_unique_to_file2:\t" + output_log::int2str(N[2]) + "\n");
}

void variant_file::output_discordance_by_indv(const string &output_file_prefix, variant_file &diff_variant_file, const string &indv_ID_map_file)
{
	LOG.printLOG("Outputting Discordance By Individual...\n");

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_variant_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, indv_ID_map_file);

	map<string, pair<int, int> > indv_sums;

	string CHROM;
	vector<char> variant_line;
	int POS;
	int s1, s2, indv1, indv2;

	entry * e1 = get_entry_object(N_indv);
	entry * e2 = diff_variant_file.get_entry_object(diff_variant_file.N_indv);
	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_entry(s1, variant_line);
			e1->reset(variant_line);
		}

		if (s2 != -1)
		{
			diff_variant_file.get_entry(s2, variant_line);
			e2->reset(variant_line);
		}

		e1->parse_basic_entry(true);
		e2->parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1->get_REF();
		string REF2 = e2->get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if ((REF.size() != REF2.size()) || ((REF != REF2) && (REF2 != "N") && (REF != "N")))
		{
			LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
			continue;
		}

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1->get_ALT();
		ALT2 = e2->get_ALT();

		bool alleles_match = (ALT == ALT2) && (REF == REF2);
		e1->parse_full_entry(true);
		e1->parse_genotype_entries(true);

		e2->parse_full_entry(true);
		e2->parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<int,int> geno_ids1, geno_ids2;
		pair<string, string> missing_genotype(".",".");
		pair<int, int> missing_id(-1,-1);

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			if (alleles_match)
			{	// Alleles match, so can compare ids instead of strings
				e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
				e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

				if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
				{
					indv_sums[combined_individuals_it->first].first++;
					if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
						((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
					{	// Match
						// Don't do anything
					}
					else
					{	// Mismatch
						indv_sums[combined_individuals_it->first].second++;
					}
				}
				else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
				{	// Both missing
					// Don't do anything.
				}
				else if (geno_ids1 != missing_id)
				{	// Genotype 1 is not missing, genotype 2 is.
					// Don't do anything.
				}
				else if (geno_ids2 != missing_id)
				{	// Genotype 2 is not missing, genotype 1 is.
					// Don't do anything.
				}
				else
					LOG.error("Unknown condition");
			}
			else
			{	// Alleles don't match, so need to be more careful and compare strings
				e1->get_indv_GENOTYPE_strings(indv1, genotype1);
				e2->get_indv_GENOTYPE_strings(indv2, genotype2);

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{	// No missing data
					indv_sums[combined_individuals_it->first].first++;
					if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
						((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
					{	// Match
						// Don't do anything
					}
					else
					{	// Mismatch
						indv_sums[combined_individuals_it->first].second++;
					}
				}
				else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
				{	// Both missing
					// Don't do anything
				}
				else if (genotype1 != missing_genotype)
				{	// Genotype 1 is not missing, genotype 2 is.
					// Don't do anything
				}
				else if (genotype2 != missing_genotype)
				{	// Genotype 2 is not missing, genotype 1 is.
					// Don't do anything
				}
				else
					LOG.error("Unknown condition");
			}
		}
	}

	string output_file = output_file_prefix + ".diff.indv";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		LOG.error("Could not open Sites Differences File: " + output_file, 3);
	out << "INDV\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	int N, N_discord;
	double discordance;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		out << combined_individuals_it->first;
		N = indv_sums[combined_individuals_it->first].first;
		N_discord = indv_sums[combined_individuals_it->first].second;
		discordance = N_discord / double(N);
		out << "\t" << N << "\t" << N_discord << "\t" << discordance << endl;
	}

	delete e1;
	delete e2;
	out.close();
}

void variant_file::output_discordance_by_site(const string &output_file_prefix, variant_file &diff_variant_file, const string &indv_ID_map_file)
{
	LOG.printLOG("Outputting Discordance By Site...\n");

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_variant_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, indv_ID_map_file);

	string CHROM;
	vector<char> variant_line;
	int POS;
	int s1, s2, indv1, indv2;
	entry * e1 = get_entry_object(N_indv);
	entry * e2 = diff_variant_file.get_entry_object(diff_variant_file.N_indv);

	string output_file = output_file_prefix + ".diff.sites";
	ofstream diffsites(output_file.c_str());
	if (!diffsites.is_open())
		LOG.error("Could not open Sites Differences File: " + output_file, 3);
	diffsites << "CHROM\tPOS\tFILES\tMATCHING_ALLELES\tN_COMMON_CALLED\tN_DISCORD\tDISCORDANCE" << endl;

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		diffsites << CHROM << "\t" << POS;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		bool data_in_both = true;
		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_entry(s1, variant_line);
			e1->reset(variant_line);
		}
		else
			data_in_both = false;

		if (s2 != -1)
		{
			diff_variant_file.get_entry(s2, variant_line);
			e2->reset(variant_line);
		}
		else
			data_in_both = false;

		if (data_in_both)
			diffsites << "\tB";
		else if ((s1 != -1) && (s2 == -1))
			diffsites << "\t1";
		else if ((s1 == -1) && (s2 != -1))
			diffsites << "\t2";
		else
			LOG.error("Unhandled condition");

		e1->parse_basic_entry(true);
		e2->parse_basic_entry(true);

		// Set the reference to the non-missing entry (if available)
		string REF = e1->get_REF();
		string REF2 = e2->get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if ((REF.size() != REF2.size()) || ((REF != REF2) && (REF2 != "N") && (REF != "N")))
		{
			LOG.one_off_warning("Non-matching REF. Skipping all such sites.");
			continue;
		}

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1->get_ALT();
		ALT2 = e2->get_ALT();

		bool alleles_match = ((ALT == ALT2) && (REF == REF2));
		diffsites << "\t" << alleles_match;

		e1->parse_full_entry(true);
		e1->parse_genotype_entries(true);

		e2->parse_full_entry(true);
		e2->parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<int,int> geno_ids1, geno_ids2;
		pair<string, string> missing_genotype(".",".");
		pair<int, int> missing_id(-1,-1);

		unsigned int N_common_called=0;	// Number of genotypes called in both files
		unsigned int N_missing_1=0, N_missing_2=0;
		unsigned int N_discord=0;
		unsigned int N_concord_non_missing=0;

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			if (alleles_match)
			{	// Alleles match, so can compare ids instead of strings
				e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
				e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

				if ((geno_ids1 != missing_id) && (geno_ids2 != missing_id))
				{
					N_common_called++;
					if (((geno_ids1.first == geno_ids2.first) && (geno_ids1.second == geno_ids2.second)) ||
						((geno_ids1.first == geno_ids2.second) && (geno_ids1.second == geno_ids2.first)) )
					{	// Match
						N_concord_non_missing++;
					}
					else
					{	// Mismatch
						N_discord++;
					}
				}
				else if ((geno_ids1 == missing_id) && (geno_ids2 == missing_id))
				{	// Both missing
					N_missing_1++; N_missing_2++;
				}
				else if (geno_ids1 != missing_id)
				{	// Genotype 1 is not missing, genotype 2 is.
					N_missing_2++;
				}
				else if (geno_ids2 != missing_id)
				{	// Genotype 2 is not missing, genotype 1 is.
					N_missing_1++;
				}
				else
					LOG.error("Unknown condition");
			}
			else
			{	// Alleles don't match, so need to be more careful and compare strings
				e1->get_indv_GENOTYPE_strings(indv1, genotype1);
				e2->get_indv_GENOTYPE_strings(indv2, genotype2);

				if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
				{	// No missing data
					N_common_called++;
					if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
						((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
					{	// Match
						N_concord_non_missing++;
					}
					else
					{	// Mismatch
						N_discord++;
					}
				}
				else if ((genotype1 == missing_genotype) && (genotype2 == missing_genotype))
				{	// Both missing
					N_missing_1++; N_missing_2++;
				}
				else if (genotype1 != missing_genotype)
				{	// Genotype 1 is not missing, genotype 2 is.
					N_missing_2++;
				}
				else if (genotype2 != missing_genotype)
				{	// Genotype 2 is not missing, genotype 1 is.
					N_missing_1++;
				}
				else
					LOG.error("Unknown condition");
			}
		}
		double discordance = N_discord / double(N_common_called);
		diffsites << "\t" << N_common_called << "\t" << N_discord << "\t" << discordance;
		diffsites << endl;
	}
	delete e1;
	delete e2;
	diffsites.close();
}

void variant_file::output_discordance_matrix(const string &output_file_prefix, variant_file &diff_variant_file, const string &indv_ID_map_file)
{
	LOG.printLOG("Outputting Discordance Matrix\n\tFor bi-allelic loci, called in both files, with matching alleles only...\n");

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_variant_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, indv_ID_map_file);

	vector<char> variant_line;
	int s1, s2, indv1, indv2;
	entry * e1 = get_entry_object(N_indv);
	entry * e2 = diff_variant_file.get_entry_object(diff_variant_file.N_indv);

	vector<vector<int> > discordance_matrix(4, vector<int>(4, 0));

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_entry(s1, variant_line);
			e1->reset(variant_line);
		}

		if (s2 != -1)
		{
			diff_variant_file.get_entry(s2, variant_line);
			e2->reset(variant_line);
		}

		e1->parse_basic_entry(true);
		e2->parse_basic_entry(true);

		if ((e1->get_N_alleles() != 2) || (e2->get_N_alleles() != 2))
			continue;

		// Set the reference to the non-missing entry (if available)
		string REF = e1->get_REF();
		string REF2 = e2->get_REF();
		if (REF == "N")
			REF = REF2;
		if (REF2 == "N")
			REF2 = REF;

		if (REF.size() != REF2.size())
			continue;

		if ((REF != REF2) && (REF2 != "N") && (REF != "N"))
			continue;

		// Do the alternative alleles match?
		string ALT, ALT2;
		ALT = e1->get_ALT();
		ALT2 = e2->get_ALT();

		bool alleles_match = (ALT == ALT2) && (REF == REF2);
		if (alleles_match == false)
			continue;

		e1->parse_full_entry(true);
		e1->parse_genotype_entries(true);

		e2->parse_full_entry(true);
		e2->parse_genotype_entries(true);

		pair<int,int> geno_ids1, geno_ids2;
		int N1, N2;

		for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			// Alleles match, so can compare ids instead of strings
			e1->get_indv_GENOTYPE_ids(indv1, geno_ids1);
			e2->get_indv_GENOTYPE_ids(indv2, geno_ids2);

			if (((geno_ids1.first != -1) && (geno_ids1.second == -1)) ||
				((geno_ids2.first != -1) && (geno_ids2.second == -1)))
			{	// Haploid
				LOG.one_off_warning("***Warning: Haploid chromosomes not counted!***");
				continue;
			}

			N1 = geno_ids1.first + geno_ids1.second;
			N2 = geno_ids2.first + geno_ids2.second;

			if ((N1 == -1) || (N1 < -2) || (N1 > 2))
				LOG.error("Unhandled case");
			if ((N2 == -1) || (N2 < -2) || (N2 > 2))
				LOG.error("Unhandled case");

			if (N1 == -2)
				N1 = 3;

			if (N2 == -2)
				N2 = 3;

			discordance_matrix[N1][N2]++;
		}
	}

	string output_file = output_file_prefix + ".diff.discordance_matrix";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		LOG.error("Could not open Discordance Matrix File: " + output_file, 3);

	out << "-\tN_0/0_file1\tN_0/1_file1\tN_1/1_file1\tN_./._file1" << endl;
	out << "N_0/0_file2\t" << discordance_matrix[0][0] << "\t" << discordance_matrix[1][0] << "\t" << discordance_matrix[2][0] << "\t" << discordance_matrix[3][0] << endl;
	out << "N_0/1_file2\t" << discordance_matrix[0][1] << "\t" << discordance_matrix[1][1] << "\t" << discordance_matrix[2][1] << "\t" << discordance_matrix[3][1] << endl;
	out << "N_1/1_file2\t" << discordance_matrix[0][2] << "\t" << discordance_matrix[1][2] << "\t" << discordance_matrix[2][2] << "\t" << discordance_matrix[3][2] << endl;
	out << "N_./._file2\t" << discordance_matrix[0][3] << "\t" << discordance_matrix[1][3] << "\t" << discordance_matrix[2][3] << "\t" << discordance_matrix[3][3] << endl;
	out.close();
	delete e1;
	delete e2;
}

void variant_file::output_switch_error(const string &output_file_prefix, variant_file &diff_variant_file, const string &indv_ID_map_file)
{
	LOG.printLOG("Outputting Phase Switch Errors...\n");

	map<pair<string, int>, pair<int, int> > CHROMPOS_to_filepos_pair;
	map<pair<string, int>, pair<int, int> >::iterator CHROMPOS_to_filepos_pair_it;
	return_site_union(diff_variant_file, CHROMPOS_to_filepos_pair);

	map<string, pair< int, int> > combined_individuals;
	map<string, pair< int, int> >::iterator combined_individuals_it;
	return_indv_union(diff_variant_file, combined_individuals, indv_ID_map_file);

	string CHROM;
	vector<char> variant_line;
	int POS;
	int s1, s2, indv1, indv2;

	entry * e1 = get_entry_object(N_indv);
	entry * e2 = diff_variant_file.get_entry_object(diff_variant_file.N_indv);

	string output_file = output_file_prefix + ".diff.switch";
	ofstream switcherror(output_file.c_str());
	if (!switcherror.is_open())
		LOG.error("Could not open Switch Error file: " + output_file, 4);
	switcherror << "CHROM\tPOS\tINDV" << endl;

	unsigned int N_combined_indv = combined_individuals.size();
	vector<int> N_phased_het_sites(N_combined_indv, 0);
	vector<int> N_switch_errors(N_combined_indv, 0);

	pair<string, string> missing_genotype(".",".");
	vector<pair<string, string> > prev_geno_file1(N_combined_indv, missing_genotype);
	vector<pair<string, string> > prev_geno_file2(N_combined_indv, missing_genotype);
	pair<string, string> file1_hap1, file1_hap2, file2_hap1;

	for (CHROMPOS_to_filepos_pair_it=CHROMPOS_to_filepos_pair.begin(); CHROMPOS_to_filepos_pair_it != CHROMPOS_to_filepos_pair.end(); ++CHROMPOS_to_filepos_pair_it)
	{
		CHROM = CHROMPOS_to_filepos_pair_it->first.first;
		POS = CHROMPOS_to_filepos_pair_it->first.second;

		s1 = CHROMPOS_to_filepos_pair_it->second.first;
		s2 = CHROMPOS_to_filepos_pair_it->second.second;

		// Read entries from file (if available)
		if (s1 != -1)
		{
			get_entry(s1, variant_line);
			e1->reset(variant_line);
		}

		if (s2 != -1)
		{
			diff_variant_file.get_entry(s2, variant_line);
			e2->reset(variant_line);
		}

		e1->parse_basic_entry(true);
		e2->parse_basic_entry(true);

		e1->parse_full_entry(true);
		e1->parse_genotype_entries(true);

		e2->parse_full_entry(true);
		e2->parse_genotype_entries(true);

		pair<string, string> genotype1, genotype2;
		pair<string, string> missing_genotype(".",".");

		unsigned int N_common_called=0;	// Number of genotypes called in both files
		unsigned int indv_count=0;

		// Bug fix applied (#3354189) - July 5th 2011
		for (combined_individuals_it=combined_individuals.begin();
				combined_individuals_it!=combined_individuals.end();
				++combined_individuals_it, indv_count++)
		{
			indv1 = combined_individuals_it->second.first;
			indv2 = combined_individuals_it->second.second;

			if ((indv1 == -1) || (indv2 == -1))
				continue;	// Individual not found in one of the files

			e1->get_indv_GENOTYPE_strings(indv1, genotype1);
			e2->get_indv_GENOTYPE_strings(indv2, genotype2);

			if ((genotype1 != missing_genotype) && (genotype2 != missing_genotype))
			{	// No missing data
				N_common_called++;
				if (((genotype1.first == genotype2.first) && (genotype1.second == genotype2.second)) ||
					((genotype1.first == genotype2.second) && (genotype1.second == genotype2.first)) )
				{	// Have a matching genotypes in files 1 and 2
					if (genotype1.first != genotype1.second)
					{	// It's a heterozgote
						char phase1, phase2;
						phase1 = e1->get_indv_PHASE(indv1);
						phase2 = e2->get_indv_PHASE(indv2);
						if ((phase1 == '|') && (phase2 == '|'))
						{	// Calculate Phasing error (switch error)
							N_phased_het_sites[indv_count]++;
							file1_hap1 = make_pair<string,string>(prev_geno_file1[indv_count].first, genotype1.first);
							file1_hap2 = make_pair<string,string>(prev_geno_file1[indv_count].second, genotype1.second);
							file2_hap1 = make_pair<string,string>(prev_geno_file2[indv_count].first, genotype2.first);

							if ((file2_hap1 != file1_hap1) && (file2_hap1 != file1_hap2))
							{	// Must be a switch error
								string indv_id;
								N_switch_errors[indv_count]++;
								if (indv1 != -1)
									indv_id = indv[indv1];
								else
									indv_id = diff_variant_file.indv[indv2];
								switcherror << CHROM << "\t" << POS << "\t" << indv_id << endl;
							}
							prev_geno_file1[indv_count] = genotype1;
							prev_geno_file2[indv_count] = genotype2;
						}
					}
				}
			}
		}
	}
	switcherror.close();

	output_file = output_file_prefix + ".diff.indv.switch";
	ofstream idiscord(output_file.c_str());
	if (!idiscord.is_open())
		LOG.error("Could not open Individual Discordance File: " + output_file, 3);

	idiscord << "INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH" << endl;
	unsigned int indv_count=0;
	double switch_error;
	string indv_id;
	for (combined_individuals_it=combined_individuals.begin(); combined_individuals_it!=combined_individuals.end(); ++combined_individuals_it)
	{
		indv1 = combined_individuals_it->second.first;
		indv2 = combined_individuals_it->second.second;

		if (indv1 != -1)
			indv_id = indv[indv1];
		else
			indv_id = diff_variant_file.indv[indv2];

		if (N_phased_het_sites[indv_count] > 0)
			switch_error = double(N_switch_errors[indv_count]) / N_phased_het_sites[indv_count];
		else
			switch_error = 0;
		idiscord << indv_id << "\t" << N_phased_het_sites[indv_count] << "\t" << N_switch_errors[indv_count] << "\t" << switch_error << endl;

		indv_count++;
	}
	delete e1;
	delete e2;
	idiscord.close();
}

