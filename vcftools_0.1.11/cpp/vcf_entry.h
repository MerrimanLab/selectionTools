/*
 * vcf_entry.h
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#ifndef VCF_ENTRY_H_
#define VCF_ENTRY_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <cassert>

#include "entry.h"
#include "output_log.h"

extern output_log LOG;

using namespace std;

class vcf_entry : public entry
{
public:
	vcf_entry(const unsigned int N_indv);
	vcf_entry(const unsigned int N_indv, const vector<char> &data_line);
	virtual ~vcf_entry();

	static string convert_line;

	void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false);
	void parse_full_entry(bool parse_FORMAT=true);
	void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false);
	void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false);
	void parse_FORMAT();

	void reset(const vector<char> &data_line);
	void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out);

	void set_ALT(const string &in);
	void set_QUAL(const double in);
	void set_FILTER(const string &FILTER_str);
	void set_FORMAT(const string &in);
	void set_INFO(const string &INFO_str);

	void add_FORMAT_entry(const string &in, unsigned int pos);

	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const string &in);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase);
	void set_indv_GENOTYPE_alleles(unsigned int indv, const pair<string, string> &in);
	void set_indv_GENOTYPE_alleles(unsigned int indv, char a1, char a2);
	void set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in);
	void set_indv_PHASE(unsigned int indv, char in);
	void set_indv_GQUALITY(unsigned int indv, double in);
	void set_indv_DEPTH(unsigned int indv, int in);
	void set_indv_GFILTER(unsigned int indv, const string &in);

	void add_indv_GFILTER(unsigned int indv, const string &in);

	static int add_INFO_descriptor(const string &in, unsigned int index);
	static int add_FILTER_descriptor(const string &in, unsigned int index);
	static int add_FORMAT_descriptor(const string &in, unsigned int index);
	static void add_CONTIG_descriptor(const string &in, unsigned int index);

	void print(ostream &out);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype);
	void print_bcf(BGZF* out);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype);

	void filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth);
	void filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality);
	void filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all = false);

private:
	string ALT_str, FILTER_str, INFO_str, FORMAT_str, QUAL_str;
	vector<string> GENOTYPE_str;

	static map<string, Field_description> INFO_map;
	static map<string, unsigned int> INFO_reverse_map;
	static map<string, string> FILTER_map;
	static map<string, unsigned int> FILTER_reverse_map;
	static map<string, Field_description> FORMAT_map;
	static map<string, unsigned int> FORMAT_reverse_map;
	static map<string, unsigned int> CONTIG_map;
};

#endif /* VCF_ENTRY_H_ */
