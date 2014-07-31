/*
 * vcf_file.h
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#ifndef BCF_FILE_H_
#define BCF_FILE_H_

#include "output_log.h"
#include "parameters.h"
#include "variant_file.h"
#include "bgzf.h"

extern output_log LOG;
using namespace std;

class bcf_file : public variant_file
{
public:
	bcf_file(const string &filename, const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index=false, bool gatk=false);

	void get_entry(unsigned int entry_num, vector<char> &out);
	entry* get_entry_object(unsigned int N_indv);

	void print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO);
	void print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO=false);
	void print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO);
	void print_bcf(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO=false, bool stream=false);

protected:
	~bcf_file();

private:
	gzFile bcf_infile_bgzf;
	FILE *bcf_infile;
	bool file_open;
	bool is_BGZF;
	bool is_GATK;
	bool big_endian;
	unsigned int gzMAX_LINE_LEN;

	int read(void *buffer, unsigned int len, size_t size);
	void read_header(bool skip_meta=false);
	void open();
	void close();
	bool eof();

	inline void read_CHROM_only(string &CHROM);
	void read_CHROM_and_POS_only(string &CHROM, int &POS);
	inline int read_CHROM_and_POS_and_skip_remainder_of_line(string &CHROM, int &POS);

	streampos get_filepos();
	void set_filepos(streampos &filepos);
	streampos get_eof();

	void scan_file(const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index=false);

};

#endif /* BCF_FILE_H_ */
