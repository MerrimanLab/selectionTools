/*
 * bcf_entry.h
 *
 *  Created on: Sep 20, 2012
 *      Author: Anthony Marcketta
 *      ($Revision: 1 $)
 */

#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <stdint.h>

#include "output_log.h"
#include "entry.h"
#include "header.h"

extern output_log LOG;

class bcf_entry : public entry {
public:
	bcf_entry(const unsigned int N_indv, const header &header_obj, const vector<char> &line);
	bcf_entry(const unsigned int N_indv, const header &header_obj);
	~bcf_entry();

	map<int, Field_description> INFO_map;
	map<int, Field_description> FILTER_map;
	map<int, Field_description> FORMAT_map;
	map<int, Field_description> CONTIG_map;
	map<string, int> CONTIG_reverse_map;
	map<string, int> FILTER_reverse_map;
	map<string, int> INFO_reverse_map;
	map<string, int> FORMAT_reverse_map;

	header entry_header;

	unsigned int N_samples;
	unsigned int N_info;
	unsigned int N_format;
	unsigned int L_shared;
	unsigned int L_indiv;
	unsigned int line_pos;

	void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false);
	void parse_full_entry(bool parse_FORMAT=true);
	void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false);
	void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false);

	void set_ALT(const int n_allele);
	void set_ALT(const string &in);
	void set_QUAL(const float &in);
	void set_FILTER();
	void set_FORMAT();
	void set_INFO();
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const vector<char> &in);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase);
	void set_indv_GENOTYPE_and_PHASE(unsigned int indv, const unsigned int &pos, const unsigned int &size);
	void set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in);
	void set_indv_GQUALITY(unsigned int indv, const vector<char> &in);
	void set_indv_GQUALITY(unsigned int indv, const float &in);
	void set_indv_DEPTH(unsigned int indv, const vector<char> &in);
	void set_indv_DEPTH(unsigned int indv, int in);
	void set_indv_GFILTER(unsigned int indv, const string &in);
	void set_indv_GFILTER(unsigned int indv, const vector<char> &in);
	void set_indv_PHASE(unsigned int indv, char in);
	void set_indv_GENOTYPE_alleles(unsigned int indv, char a1, char a2);
	void set_indv_GENOTYPE_alleles(unsigned int indv, const pair<int, int> &in);
	void reset(const vector<char> &data_line);

	void add_FORMAT_entry(const string &in, const unsigned int &fmt_key, const unsigned int &pos, const unsigned int &line_pos, const unsigned int &type, const unsigned int &size);
	void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out);
	void read_indv_generic_entry(unsigned int indv, const int &idx, string &out);
	void write_out( const char * filename, const bool stream );
	void read_all_entries(string &out, const vector<bool> &include_indv, const vector<bool> &include_genotype);

	void filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality);
	void filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth);
	void filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all = false);

	void print(ostream &out);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype);
	void print_bcf(BGZF* out);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype);

	static int add_INFO_descriptor(const string &in, int index);
	static int add_FILTER_descriptor(const string &in, int index);
	static int add_FORMAT_descriptor(const string &in, int index);
	static void add_CONTIG_descriptor(const string &in, int index);

private:
	vector<char> line;

	vector<char> INFO_str, QUAL_str;
	vector<int> FILTER_str;
	vector<string> ALT_str;

	unsigned int INFO_pos, FILTER_pos, ALT_pos, FORMAT_pos;
};
