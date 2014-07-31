/*
 * entry.h
 *
 *  Created on: Dec 12, 2012
 *      Author: amarcketta
 */
#ifndef ENTRY_H_
#define ENTRY_H_

#include <cstring>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <vector>
#include <stdint.h>

#include <cassert>
#include "bgzf.h"
#include "output_log.h"

using namespace std;
extern output_log LOG;

enum Type_enum {Integer=0, Float=1, Character=2, String=3, Flag=4};

class Field_description
{
public:
	string ID;
	int N_entries;
	string N_entries_str;
	string Type_str;
	Type_enum Type;
	string Description;
	string Length;
	string Assembly;

	Field_description() : ID(""), N_entries(0), Type(Integer), Description("") {};
	~Field_description() {};
};

class entry
{
public:
	unsigned int N_indv;

	virtual void parse_basic_entry(bool parse_ALT=false, bool parse_FILTER=false, bool parse_INFO=false) = 0;
	virtual void parse_full_entry(bool parse_FORMAT=true) = 0;
	virtual void parse_genotype_entry(unsigned int indv, bool GT=false, bool GQ=false, bool DP=false, bool FT=false) = 0;
	virtual void parse_genotype_entries(bool GT=false, bool GQ=false, bool DP=false, bool FT=false) = 0;

	virtual void reset(const vector<char> &data_line) = 0;

	string get_CHROM() const;
	void get_CHROM(string &out) const;
	int get_POS() const;
	string get_ID() const;
	string get_REF() const;
	string get_ALT() const;
	string get_ALT_allele(int allele_num) const;
	void get_allele(int allele_num, string &out) const;
	string get_allele(int allele_num) const;
	void get_alleles_vector(vector<string> &out) const;
	string get_FILTER() const;
	void get_FILTER_vector(vector<string> &out) const;
	double get_QUAL() const;
	string get_INFO(const set<string> &INFO_to_keep, bool keep_all_INFO=false) const;
	string get_INFO_value(const string &key) const;
	string get_FORMAT() const;

	void get_indv_GENOTYPE_ids(unsigned int indv, pair<int, int> &out) const;
	void get_indv_GENOTYPE_strings(unsigned int indv, pair<string, string> &out) const;
	char get_indv_PHASE(unsigned int indv) const;
	double get_indv_GQUALITY(unsigned int indv) const;
	int get_indv_DEPTH(unsigned int indv) const;
	void get_indv_GFILTER(unsigned int indv, string &out) const;
	void get_indv_GFILTER_vector(unsigned int indv, vector<string> &out) const;
	int get_indv_ploidy(unsigned int indv) const;

	bool is_SNP() const;
	bool is_biallelic_SNP() const;
	bool is_diploid(const vector<bool> &include_indv, const vector<bool> &include_genotype) const;
	virtual void read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out) = 0;
	bool FORMAT_id_exists(const string &FORMAT_id);

	void get_allele_counts(vector<int> &out, unsigned int &N_non_missing_chr_out, const vector<bool> &include_indv, const vector<bool> &include_genotype) const;
	void get_genotype_counts(const vector<bool> &include_indv, const vector<bool> &include_genotype, unsigned int &out_N_hom1, unsigned int &out_N_het, unsigned int &out_N_hom2) const;
	unsigned int get_N_alleles() const;
	unsigned int get_N_chr(const vector<bool> &include_indv, const vector<bool> &include_genotype) const;

	void get_POS_binary(vector<char> &out) const;
	void get_ID_binary(vector<char> &out);
	void get_rlen(vector<char> &out) const;
	void get_QUAL_binary(vector<char> &out) const;
	void get_n_allele_info(vector<char> &out) const;
	void get_n_fmt_sample(vector<char> &out) const;
	void get_ALLELES_binary(vector<char> &out);
	vector<pair<string, string> > get_INFO_vector(const set<string> &INFO_to_keep, bool keep_all_INFO=false) const;
	void get_FORMAT_binary(vector<char> &out) const;

	string get_typed_string( unsigned int * line_position, const vector<char>& line );
	void get_type(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size);
	vector<int> get_int_vector(unsigned int * line_position, const vector<char>& line);
	int get_typed_int(unsigned int * line_position, const vector<char>& line, unsigned int &type, unsigned int &size);
	float get_typed_float(unsigned int * line_position, const vector<char>& line);
	vector<float> get_typed_float_vector(unsigned int * line_position, const vector<char>& line);
	void get_number(uint32_t &out, unsigned int * line_position, const vector<char>& line);
	void decode_genotype(int8_t in, int &GT, bool &phased);

	void make_typed_string(vector<char> &out, const string &in, bool typed);
	void make_typed_int(vector<char> &out, const int &in, bool typed);
	void make_int(vector<char> &out, const int &in, int type);
	void make_typed_int_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_int_vector(vector<char> &out, const string &in, int number = -1);
	void make_typed_int_vector(vector<char> &out, const vector<int> &in);
	void make_typed_float_vector(vector<char> &out, const string &in, int number = -1);
	void make_typed_float_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_string_vector(vector<char> &out, const vector<string> &in, int number = -1);
	void make_typed_GT_vector(vector<char> &out, vector<string> &in);
	void make_type_size(vector<char> &out, const unsigned int &type, const unsigned int &size);
	void encode_genotype(vector<char> &out, string &in, int exp_size);

	void copy_object(vector<char> &out, int &position, const vector<char> &in);
	void skip_section(unsigned int *line_position, const vector<char> &line);
	bool check_missing(unsigned int line_position, const unsigned int type, const vector<char> &line);

	void set_CHROM(const string &in);
	void set_POS(const int in);
	void set_ID(const string &in);
	void set_REF(const string &in);

	void add_ALT_allele(const string &in);
	void add_FILTER_entry(const string &in);

	virtual void print(ostream &out) = 0;
	virtual void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO=false) = 0;
	virtual void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype) = 0;
	virtual void print_bcf(BGZF* out) = 0;
	virtual void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO=false) = 0;
	virtual void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype) = 0;

	virtual void filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth) = 0;
	virtual void filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality) = 0;
	virtual void filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all = false) = 0;

	static double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);
	static void tokenize(const string &in, char token, vector<string> &out);

	static int str2int(const string &in, const int missing_value=-1);
	static double str2double(const string &in, const double missing_value=-1.0);

	static string int2str(const int in, const int missing_value=-1);
	static string double2str(const double in, const double missing_value=-1.0);
	static inline bool is_big_endian() { long one= 1; return !(*((char *)(&one))); };

protected:
	istringstream data_stream;

	bool basic_parsed;
	bool fully_parsed;
	bool parsed_ALT;
	bool parsed_FILTER;
	bool parsed_INFO;
	bool parsed_FORMAT;
	bool parsed_FORMAT_binary;

	string CHROM;
	int POS;
	string ID;
	string REF;
	vector<string> ALT;
	double QUAL;
	vector<string> FILTER;
	bool passed_filters;
	vector<pair<string, string> > INFO;
	vector<string> FORMAT;
	vector<char> FORMAT_binary;
	int N_INFO_removed;
	int N_FORMAT_removed;

	vector< pair<int,int> > GENOTYPE;
	vector<int> ploidy;
	vector<char>  PHASE;
	vector<double> GQUALITY;
	vector<int>   DEPTH;
	vector< vector<string> > GFILTER;

	vector<bool> parsed_GT;
	vector<bool> parsed_GQ;
	vector<bool> parsed_DP;
	vector<bool> parsed_FT;

	map<string, unsigned int> FORMAT_to_idx;
	int GT_idx;
	int GQ_idx;
	int DP_idx;
	int FT_idx;

	vector<unsigned int> FORMAT_positions, FORMAT_types, FORMAT_sizes, FORMAT_skip, FORMAT_keys;
};

#endif /* ENTRY_H_ */
