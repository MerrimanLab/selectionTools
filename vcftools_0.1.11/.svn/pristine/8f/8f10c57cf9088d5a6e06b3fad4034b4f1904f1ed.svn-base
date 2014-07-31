/*
 * vcf_file.h
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#ifndef VCF_FILE_H_
#define VCF_FILE_H_

#include "output_log.h"
#include "vcf_entry.h"
#include "parameters.h"
#include "variant_file.h"


extern output_log LOG;

using namespace std;

class vcf_file : public variant_file
{
public:
	vcf_file(const string &filename, bool compressed, const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index=false);
	vcf_file();

	void get_entry(unsigned int entry_num, vector<char> &out);
	entry* get_entry_object(unsigned int N_indv);

	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO);
	void print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO);
	void print_bcf(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO=false, bool stream=false);

	inline char peek();
	streampos get_filepos();
	void set_filepos(streampos &filepos);

protected:
	~vcf_file();

private:
	ifstream vcf_in;
	gzFile gzvcf_in;
	char *gz_readbuffer;
	unsigned int gzMAX_LINE_LEN;
	unsigned int contig_index;
	void open();
	void close();
	bool eof();
	inline void read_line(string &out);
	inline void read_line(vector<char> &out);
	inline void read_CHROM_only(string &CHROM);
	void read_CHROM_and_POS_only(string &CHROM, int &POS);
	inline int read_CHROM_and_POS_and_skip_remainder_of_line(string &CHROM, int &POS);
	void parse_header(const string &line);
	void parse_meta(const string &line, unsigned int &line_index);
	void scan_file(const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index=false);
};

#endif /* VCF_FILE_H_ */
