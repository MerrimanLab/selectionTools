/*
 * bcf_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "bcf_file.h"

bcf_file::bcf_file(const string &fname, const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index, bool gatk)
{
	filename = fname; has_body = false;
	has_file_format = false; has_header = false;
	has_meta = false; file_open = false;
	is_BGZF = false; bcf_format = true;
	big_endian = is_big_endian();
	is_GATK = gatk;
	header_obj = header();

	open();
	scan_file(chrs_to_keep, exclude_chrs, force_write_index);
}

bcf_file::~bcf_file()
{
	close();
}

void bcf_file::open()
{
	int ret;
	char magic[5];
	char test_str[3] = {'B','C','F'};

	if (filename.substr(filename.size()-4) == ".vcf")
			LOG.error("Filename ends in '.vcf'. Shouldn't you be using --vcf?\n");

	if (filename.substr(filename.size()-7) == ".bcf")
			LOG.error("Filename ends in '.vcf.gz'. Shouldn't you be using --gzvcf?\n");

	ret = bgzf_is_bgzf(filename.c_str());

	if (ret == 1)
		is_BGZF = true;
	else
		is_BGZF = false;

	bcf_infile_bgzf = NULL;	bcf_infile = NULL;
	if (is_BGZF){
		gzMAX_LINE_LEN = 1024*1024;
		bcf_infile_bgzf = gzopen(filename.c_str(), "rb");
		if (bcf_infile_bgzf == NULL)
			LOG.error("Could not open BGZF BCF file: " + filename, 0);
		#ifdef ZLIB_VERNUM
			string tmp(ZLIB_VERSION);
			LOG.printLOG("Using zlib version: " + tmp + "\n");
			#if (ZLIB_VERNUM >= 0x1240)
				ret = gzbuffer(bcf_infile_bgzf, gzMAX_LINE_LEN); // Included in zlib v1.2.4 and makes things MUCH faster
			#else
				LOG.printLOG("Versions of zlib >= 1.2.4 will be *much* faster when reading compressed BCF files.\n");
			#endif
		#endif
	}
	else
		bcf_infile = fopen(filename.c_str(), "r");
	if ((bcf_infile == NULL) && (bcf_infile_bgzf == NULL))
	{
		LOG.error("Could not open BCF file\n");
	}

	read(magic, 5, 1);

	if ( (magic[0] != test_str[0]) || (magic[1] != test_str[1]) || (magic[2] != test_str[2]) )
		LOG.error("Does not appear to be a BCF file\n");

	if ( ((int)magic[3] != 2) || ((int)magic[4] != 1 ) )
	{
		stringstream tmp_stream;
		tmp_stream << "File version number: " << (int)magic[3] << "." << (int)magic[4] << "\n";
		LOG.printLOG( tmp_stream.str() );
		LOG.error("VCFtools is currently only compatible with BCFv2.1\n");
	}
}

streampos bcf_file::get_filepos()
{
	if (!is_BGZF)
		return ftell(bcf_infile);
	else
		return gztell(bcf_infile_bgzf);
}

void bcf_file::set_filepos(streampos &filepos)
{
	if (!is_BGZF)
		fseek(bcf_infile, filepos, SEEK_SET);
	else{
		gzseek(bcf_infile_bgzf, filepos, SEEK_SET);
	}
}

void bcf_file::close()
{
	if (!is_BGZF)
		fclose(bcf_infile);
	else{
		gzclose(bcf_infile_bgzf);
	}
}

void bcf_file::scan_file(const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index)
{
	bool filter_by_chr = (chrs_to_keep.size() != 0);
	bool exclude_by_chr = (exclude_chrs.size() != 0);
	string index_filename = filename + ".bcfidx";
	bool could_read_index_file = false;

	N_entries = 0;
	if (force_write_index == false)
		could_read_index_file = read_index_file(index_filename);

	string CHROM, last_CHROM="";
	streampos filepos, endpos;

	if (could_read_index_file == false)
	{
		int POS, last_POS = -1;
		char magic[5];

		endpos = get_eof();
		read(magic, 5, 1);
		read_header();

		if ((has_header == false) || (has_meta == false))
			LOG.error("No header or meta information. Invalid file: " + filename);

		filepos = get_filepos();
		if (!is_BGZF)
		{
			while( filepos<endpos )
			{
				read_CHROM_and_POS_and_skip_remainder_of_line(CHROM, POS);

				if (last_CHROM != CHROM)
				{
					LOG.printLOG("\tScanning Chromosome: " + CHROM + "\n");
					last_CHROM = CHROM;
				}
				if (POS == last_POS)
					LOG.one_off_warning("\tWarning - file contains entries with the same position. These entries will be processed separately.\n");

				last_POS = POS;

				entry_file_locations.push_back(filepos);
				N_entries++;
				filepos = get_filepos();
			}
		}
		else
		{
			int ret = 1;
			while( (ret!=0) and (ret!=-1) )
			{
				ret = read_CHROM_and_POS_and_skip_remainder_of_line(CHROM, POS);

				if ((ret==0) or (ret==-1))
					break;

				if (last_CHROM != CHROM)
				{
					LOG.printLOG("\tScanning Chromosome: " + CHROM + "\n");
					last_CHROM = CHROM;
				}
				if (POS == last_POS)
					LOG.one_off_warning("\tWarning - file contains entries with the same position. These entries will be processed separately.\n");

				last_POS = POS;

				entry_file_locations.push_back(filepos);
				N_entries++;
				filepos = get_filepos();
			}
		}
	write_index_file(index_filename);
	}
	else
		read_header(true);

	if (header_obj.CONTIG_map.empty())
		LOG.error("BCF file does not contain any contig declarations.\n", 0);

	LOG.printLOG("File contains " + output_log::int2str(N_entries) + " entries and " + output_log::int2str(N_indv) + " individuals.\n");

	if ((exclude_by_chr == true) || (filter_by_chr == true))
	{
		unsigned int N_found_required_chr = chrs_to_keep.size();

		LOG.printLOG("Filtering by chromosome.\n");
		for (unsigned int ui=0; ui<N_entries; ui++)
		{
			if ((filter_by_chr == true) && (N_found_required_chr == 0))
			{
				LOG.printLOG("Skipping Remainder.\n");
				entry_file_locations.erase(entry_file_locations.begin()+ui, entry_file_locations.end());
				break;
			}

			set_filepos(entry_file_locations[ui]);
			read_CHROM_only(CHROM);

			if (last_CHROM != CHROM)
			{
				LOG.printLOG("\tChromosome: " + CHROM + "\n");
				if ((filter_by_chr == true) && (chrs_to_keep.find(last_CHROM) != chrs_to_keep.end()))
					N_found_required_chr--;

				last_CHROM = CHROM;
			}
			if ((exclude_by_chr == true) && (exclude_chrs.find(CHROM) != exclude_chrs.end()))
			{
				entry_file_locations[ui] = -1;
				continue;
			}
			if ((filter_by_chr == true) && (chrs_to_keep.find(CHROM) == chrs_to_keep.end()))
			{
				entry_file_locations[ui] = -1;
				continue;
			}
		}
		sort(entry_file_locations.begin(), entry_file_locations.end());
		while((entry_file_locations.size() > 0) && (entry_file_locations[0] < 0))
			entry_file_locations.pop_front();

		N_entries = entry_file_locations.size();
		LOG.printLOG("Keeping " + output_log::int2str(N_entries) + " entries on specified chromosomes.\n");
	}

	include_indv.clear();
	include_indv.resize(N_indv, true);
	include_entry.clear();
	include_entry.resize(N_entries, true);
	include_genotype.clear();
	include_genotype.resize(N_entries, vector<bool>(N_indv, true));

}

void bcf_file::read_CHROM_only(string &CHROM)
{
	int32_t chrom_int[3];

	read(&chrom_int[0], 3, sizeof(int32_t) );
	CHROM = header_obj.CONTIG_map[chrom_int[2]].ID;
}

void bcf_file::read_CHROM_and_POS_only(string &CHROM, int &POS)
{
	int32_t chrom_int[4];//, pos_int;

	read(&chrom_int[0], 4, sizeof(int32_t) );
	CHROM = header_obj.CONTIG_map[chrom_int[2]].ID;
	POS = chrom_int[3] + (int32_t)1;
}

int bcf_file::read_CHROM_and_POS_and_skip_remainder_of_line(string &CHROM, int &POS)
{
	int ret;
	int32_t chrom_int[4];

	ret = read(&chrom_int[0], 4, sizeof(int32_t) );
	if (ret != 4*sizeof(int32_t)) return 0;
	CHROM = header_obj.CONTIG_map[chrom_int[2]].ID;
	POS = chrom_int[3] + (int32_t)1;

	size_t forward = chrom_int[0] + chrom_int[1] - 2*sizeof(int32_t);
	char whole_line[forward];

	ret = read(&whole_line, 1, forward);
	return ret;
}

void bcf_file::get_entry(unsigned int entry_num, vector<char> &out)
{
	uint32_t size_int[2];
	int read_size = 0;
	set_filepos(entry_file_locations[entry_num]);

	read(&size_int[0], 2, sizeof(uint32_t) );
	read_size = size_int[0] + size_int[1];

	out.resize(read_size+2*sizeof(uint32_t));
	memcpy(&out[0], size_int, 2*sizeof(uint32_t));
	read(&out[2*sizeof(uint32_t)], 1, read_size);
}

entry* bcf_file::get_entry_object(unsigned int N_indv)
{
	return new bcf_entry(N_indv, header_obj);
}

int bcf_file::read(void *buffer, unsigned int len, size_t size)
{
	int ret;
	if (is_BGZF)
		ret = gzread(bcf_infile_bgzf, buffer, size*len);
	else
		ret = fread(buffer, 1, size*len, bcf_infile);

	if ((big_endian) && (size > 1)) // Note: don't both swapping character arrays - BCF is defined as little endian.
	{
		unsigned int ui;
		for (ui=0; ui<len; ui++)
			ByteSwap((unsigned char *)buffer+(size*ui), size);
	}
	return ret;
}

void bcf_file::read_header(bool skip_meta)
{
	uint32_t len_text;
	vector<string> headers;
	vector<string> ID_lookup;
	unsigned int N_header_indv = 0;

	read(&len_text, 1, sizeof(uint32_t));
	char *header_array = new char[(unsigned int)len_text];
	read(header_array, len_text, 1);
	string header(header_array);
	delete [] header_array;

	int contig_count = 0;
	header.erase( header.find_last_not_of(" \f\n\r\t\v\0" ) + 1 );

	istringstream iss(header);
	string line;
	while (getline(iss, line))
		headers.push_back(line);

	if (headers.size() == 0 )
	{
		LOG.error(" Input BCF file does not have a header.\n");
		exit(0);
	}
	else
		has_meta = true;

	// It needs to parse the header information and store in a
	// structure for later access.
	ID_lookup.resize(headers.size(), "");
	int pos_correct = 0;
	pos_correct += header_obj.add_FILTER_descriptor("ID=PASS,Description=PASS", pos_correct);

	for (unsigned int ui=0; ui<headers.size(); ui++)
	{
		if ( (skip_meta == false) and (headers[ui].substr(0, 6) != "#CHROM") )
			meta.push_back(headers[ui]);

		if (headers[ui].substr(0, 6) == "##ALT=")
		{
			if (is_GATK)
				pos_correct += 1;
		}
		else if (headers[ui].substr(0, 9) == "##FORMAT=")
		{
			string tmp = headers[ui].substr(10, headers[ui].size()-8);
			pos_correct += header_obj.add_FORMAT_descriptor(tmp, pos_correct);

		}
		else if (headers[ui].substr(0, 9) == "##FILTER=")
		{
			string tmp = headers[ui].substr(10, headers[ui].size()-8);
			pos_correct += header_obj.add_FILTER_descriptor(tmp, pos_correct);
		}
		else if (headers[ui].substr(0, 7) == "##INFO=")
		{
			string tmp = headers[ui].substr(8, headers[ui].size()-8);
			pos_correct += header_obj.add_INFO_descriptor(tmp, pos_correct);
		}
		else if (headers[ui].substr(0, 9) == "##contig=")
		{
			string tmp = headers[ui].substr(10, headers[ui].size()-8);
			header_obj.add_CONTIG_descriptor(tmp, contig_count);
			contig_count += 1;
		}
		else if (headers[ui].substr(0, 6) == "#CHROM")
		{
			vector<string> tmp;
			has_header = true;
			entry::tokenize(headers[ui],'\t',tmp);
			for (unsigned int ui = 0; ui < tmp.size(); ui++)
			{
				switch (ui)
				{
					case 0: if (tmp[ui] != "#CHROM") LOG.warning("First Header entry should be #CHROM: " + tmp[ui]); break;
					case 1: if (tmp[ui] != "POS") LOG.warning("Second Header entry should be POS: " + tmp[ui]); break;
					case 2: if (tmp[ui] != "ID") LOG.warning("Third Header entry should be ID: " + tmp[ui]); break;
					case 3: if (tmp[ui] != "REF") LOG.warning("Fourth Header entry should be REF: " + tmp[ui]); break;
					case 4: if (tmp[ui] != "ALT") LOG.warning("Fifth Header entry should be ALT: " + tmp[ui]); break;
					case 5: if (tmp[ui] != "QUAL") LOG.warning("Sixth Header entry should be QUAL: " + tmp[ui]); break;
					case 6: if (tmp[ui] != "FILTER") LOG.warning("Seventh Header entry should be FILTER: " + tmp[ui]); break;
					case 7: if (tmp[ui] != "INFO") LOG.warning("Eighth Header entry should be INFO: " + tmp[ui]); break;
					case 8:
						if (tmp[ui] != "FORMAT")
							LOG.warning("Ninth Header entry should be FORMAT: " + tmp[ui]);
						else
							has_genotypes = true;
						break;
					default:
					{
						if (ui <= 8)
							LOG.error("Incorrectly formatted header.");
						indv.push_back(tmp[ui]);
						N_header_indv++;
					}
					break;
				}
			}
			N_indv = N_header_indv;
		}
	}
}

bool bcf_file::eof()
{
	if (is_BGZF){
		return gzeof(bcf_infile_bgzf);
	}
	else
		return(feof(bcf_infile));
}

streampos bcf_file::get_eof()
{
	streampos end_pos;
	if (!is_BGZF)
	{
		fseek(bcf_infile, 0, SEEK_END);
		end_pos = get_filepos();
		fseek(bcf_infile, 0, SEEK_SET);
	}
	else
	{
		gzseek(bcf_infile_bgzf, 0, SEEK_END);
		end_pos = get_filepos();
		gzseek(bcf_infile_bgzf, 0, SEEK_SET);
	}
	return end_pos;
}

void bcf_file::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	for (unsigned int ui=0; ui<meta.size(); ui++)
		out << meta[ui] << endl;

	out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (N_indv > 0)
		out << "\tFORMAT";
	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui])
			out << "\t" << indv[ui];
	out << endl;

	vector<char> variant_line;
	entry *e = new bcf_entry(N_indv, header_obj);
	for (unsigned int s=0; s<N_entries; s++)
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true, true, true);
			e->parse_full_entry(true);
			e->parse_genotype_entries(true);
			e->print(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype[s]);
		}
	delete e;
}

void bcf_file::print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	LOG.printLOG("Outputting VCF file... ");

	string output_file = output_file_prefix + ".recode.vcf";
	ofstream out(output_file.c_str());
	if (!out.is_open())
		LOG.error("Could not open VCF Output File: " + output_file, 3);

	print(out, INFO_to_keep, keep_all_INFO);

	out.close();
	LOG.printLOG("Done\n");
}

void bcf_file::print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	string header_str;
	uint32_t len_text = 0;
	vector<char> header;

	char magic[5] = {'B','C','F','\2','\1'};
	bgzf_write(out, magic, 5);

	for (unsigned int ui=0; ui<meta.size(); ui++)
	{
		for (unsigned int uj=0; uj<meta[ui].length(); uj++)
			header.push_back( meta[ui][uj] );
		header.push_back('\n');
	}

	header_str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (N_indv > 0)
		header_str += "\tFORMAT";

	for (unsigned int ui=0; ui<N_indv; ui++)
		if (include_indv[ui])
		{
			header_str += "\t";
			header_str += indv[ui];
		}
	header_str += "\n";

	for (unsigned int ui=0; ui<header_str.length(); ui++)
		header.push_back( header_str[ui] );
	header.push_back( '\0' );
	len_text = header.size();

	bgzf_write(out, (char *)&len_text, sizeof(len_text) );
	bgzf_write(out, (char *)&header[0], len_text );
	vector<char> variant_line;
	entry * e = new bcf_entry(N_indv, header_obj);
	for (unsigned int s=0; s<N_entries; s++)
	{
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true, true, true);
			e->parse_full_entry(true);
			e->parse_genotype_entries(true);
			e->print_bcf(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype[s]);
		}
	}
	delete e;
}

void bcf_file::print_bcf(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO, bool stream)
{
	LOG.printLOG("Outputting BCF file... ");
	BGZF * out;
	if(!stream)
	{
		string output_file = output_file_prefix + ".recode.bcf";
		out = bgzf_open(output_file.c_str(), "w");
	}
	else
		out = bgzf_dopen(1, "w");

	print_bcf(out, INFO_to_keep, keep_all_INFO);

	bgzf_close(out);
	LOG.printLOG("Done\n");
}
