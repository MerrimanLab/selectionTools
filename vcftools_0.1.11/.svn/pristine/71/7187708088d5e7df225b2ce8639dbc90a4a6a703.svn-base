/*
 * vcf_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "vcf_file.h"

vcf_file::vcf_file(const string &fname, bool comp, const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index)
{
	filename = fname; compressed = comp;
	has_body = false; has_file_format = false;
	has_header = false; has_meta = false;
	bcf_format = false; has_genotypes = false;
	has_contigs = false; contig_index = 0;
	gzMAX_LINE_LEN = 0; meta.clear();

	open();
	scan_file(chrs_to_keep, exclude_chrs, force_write_index);
}

vcf_file::vcf_file()
{
	gzvcf_in = false;
	gz_readbuffer = NULL;
	compressed = false;	has_body = false;
	has_file_format = false; has_header = false;
	has_meta = false; has_genotypes = false;
	has_contigs = false; contig_index = 0;
	gzMAX_LINE_LEN = 0; meta.clear();
}

vcf_file::~vcf_file()
{
	close();
}

// Parse VCF meta information
void vcf_file::parse_meta(const string &line, unsigned int &line_index)
{
	has_meta = true;
	meta.push_back(line);
	size_t found=line.find("##fileformat=");
	if (found!=string::npos)
	{
		has_file_format = true;
		found = line.find_first_of("=");
		string version = line.substr(found+1);
		if ((version != "VCFv4.0") && (version != "VCFv4.1"))
			LOG.error("VCF version must be v4.0 or v4.1:\nYou are using version " + version);
	}

	found=line.find("##INFO=");
	if (found!=string::npos)
	{	// Found an INFO descriptor
		line_index += vcf_entry::add_INFO_descriptor(line, line_index);
	}

	found=line.find("##FILTER=");
	if (found!=string::npos)
	{	// Found a FILTER descriptor
		line_index += vcf_entry::add_FILTER_descriptor(line, line_index);
	}

	found=line.find("##FORMAT=");
	if (found!=string::npos)
	{	// Found a genotype filter descriptor
		line_index += vcf_entry::add_FORMAT_descriptor(line, line_index);
	}

//ALT FIELDS NO LONGER COUNT
//	found=line.find("##ALT=");
//	if (found!=string::npos)
//		line_index += 1;

	found=line.find("##contig=");
	if (found!=string::npos)
	{	// Found a contig descriptor
		vcf_entry::add_CONTIG_descriptor(line, contig_index);
		contig_index++;
		has_contigs = true;
	}
}

void vcf_file::parse_header(const string &line)
{
	// #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	(FORMAT	NA00001 NA00002 ... )
	if (has_header == true)
		LOG.warning("Multiple Header lines.");

	has_header = true;
	istringstream header(line);
	int count = 0;
	string tmp_str;
	unsigned int N_header_indv = 0;
	has_genotypes = false;
	while (!header.eof())
	{
		getline(header, tmp_str, '\t');
		switch (count)
		{
			case 0: if (tmp_str != "#CHROM") LOG.warning("First Header entry should be #CHROM: " + tmp_str); break;
			case 1: if (tmp_str != "POS") LOG.warning("Second Header entry should be POS: " + tmp_str); break;
			case 2: if (tmp_str != "ID") LOG.warning("Third Header entry should be ID: " + tmp_str); break;
			case 3: if (tmp_str != "REF") LOG.warning("Fourth Header entry should be REF: " + tmp_str); break;
			case 4: if (tmp_str != "ALT") LOG.warning("Fifth Header entry should be ALT: " + tmp_str); break;
			case 5: if (tmp_str != "QUAL") LOG.warning("Sixth Header entry should be QUAL: " + tmp_str); break;
			case 6: if (tmp_str != "FILTER") LOG.warning("Seventh Header entry should be FILTER: " + tmp_str); break;
			case 7: if (tmp_str != "INFO") LOG.warning("Eighth Header entry should be INFO: " + tmp_str); break;
			case 8:
				if (tmp_str != "FORMAT")
					LOG.warning("Ninth Header entry should be FORMAT: " + tmp_str);
				else
					has_genotypes = true;
				break;
			default:
			{
				if (count <= 8)
					LOG.error("Incorrectly formatted header.");
				indv.push_back(tmp_str);
				N_header_indv++;
			}
			break;
		}
		count++;
	}
	N_indv = N_header_indv;

	if ((has_genotypes == true ) && (N_indv == 0))
		LOG.warning("FORMAT field without genotypes?");
}

void vcf_file::scan_file(const set<string> &chrs_to_keep, const set<string> &exclude_chrs, bool force_write_index)
{
	bool filter_by_chr = (chrs_to_keep.size() != 0);
	bool exclude_by_chr = (exclude_chrs.size() != 0);
	string index_filename = filename + ".vcfidx";
	bool could_read_index_file = false;
	if (force_write_index == false)
		could_read_index_file = read_index_file(index_filename);
	string CHROM, last_CHROM="";
	unsigned int meta_counter = 1;

	if (could_read_index_file == false)
	{
		int POS, last_POS = -1;
		bool found_header = false;
		bool found_meta = false;

		LOG.printLOG("Building new index file.\n");
		string line, CHROM, last_CHROM = "";
		streampos filepos;
		char c;
		N_entries=0;
		N_indv = 0;

		while (!eof())
		{
			filepos = get_filepos();
			c = peek();

			if ((c == '\n') || (c == '\r'))
			{
				read_line(line);
				continue;
			}
			else if (c == EOF)
				break;

			if (c == '#')
			{
				read_line(line);
				if (line[1] == '#')
				{	// Meta information
					parse_meta(line, meta_counter);
					found_meta = true;
				}
				else
				{	// Must be header information: #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	(FORMAT	NA00001 NA00002 ... )
					parse_header(line);
					found_header = true;
				}
			}
			else
			{	// Must be a data line
				if ((found_header == false) || (found_meta == false))
					LOG.error("No header or meta information. Invalid file: " + filename);

				read_CHROM_and_POS_and_skip_remainder_of_line(CHROM, POS);

				if (POS == last_POS)
				{
					if (last_CHROM == CHROM)
						LOG.one_off_warning("\tWarning - file contains entries with the same position. These entries will be processed separately.\n");
				}
				else if (last_POS > POS)
				{
					if (last_CHROM == CHROM)
						LOG.error(" VCF file is not sorted at position " + CHROM + ":" + LOG.int2str(POS) + ".\n");
				}

				if (last_CHROM != CHROM)
				{
					LOG.printLOG("\tScanning Chromosome: " + CHROM + "\n");
					last_CHROM = CHROM;
				}

				last_POS = POS;
				entry_file_locations.push_back(filepos);
				N_entries++;
			}
		}

		if ((found_header == false) || (found_meta == false))
			LOG.error("No header or meta information. Invalid file: " + filename);

		write_index_file(index_filename);
	}
	else
	{
		vector<string> meta_lines = meta; meta.resize(0);
		meta_counter = 1;
		for (unsigned int ui=0; ui<meta_lines.size(); ui++)
			parse_meta(meta_lines[ui], meta_counter);
	}

	if (has_contigs == false)
	{
		vector<string> contig_vector;
		get_default_contigs(contig_vector);
		for(unsigned int ui=0; ui<contig_vector.size(); ui++)
			vcf_entry::add_CONTIG_descriptor(contig_vector[ui], ui);
	}

	LOG.printLOG("File contains " + output_log::int2str(N_entries) + " entries and " + output_log::int2str(N_indv) + " individuals.\n");
	has_genotypes = (N_indv > 0);

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

void vcf_file::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO)
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
	entry * e = new vcf_entry(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true, true, true);
			e->parse_full_entry(true);
			e->parse_genotype_entries(true,true,true,true);
			e->print(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype[s]);
		}
	delete e;
}

void vcf_file::print(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO)
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

void vcf_file::print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	string header_str;
	uint32_t len_text = 0;
	vector<char> header;

	char magic[5] = {'B','C','F','\2', '\1'};
	bgzf_write(out, magic, 5);

	for (unsigned int ui=0; ui<meta.size(); ui++)
	{
		for (unsigned int uj=0; uj<meta[ui].length(); uj++)
			header.push_back( meta[ui][uj] );
		header.push_back('\n');
	}

	if (has_contigs == false)
	{
		vector<string> contig_vector;
		get_default_contigs(contig_vector);
		for(unsigned int ui=0; ui<contig_vector.size(); ui++)
		{
			for(unsigned int uj=0; uj<contig_vector[ui].size(); uj++)
				header.push_back(contig_vector[ui][uj]);
			header.push_back('\n');
		}
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
	entry * e = new vcf_entry(N_indv);
	for (unsigned int s=0; s<N_entries; s++)
		if (include_entry[s] == true)
		{
			get_entry(s, variant_line);
			e->reset(variant_line);
			e->parse_basic_entry(true, true, true);
			e->parse_full_entry(true);
			e->parse_genotype_entries(true,true,true,true);
			e->print_bcf(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype[s]);
		}
	delete e;
}

void vcf_file::print_bcf(const string &output_file_prefix, const set<string> &INFO_to_keep, bool keep_all_INFO, bool stream)
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

void vcf_file::open()
{
	struct stat buf;

	int i = stat(filename.c_str(), &buf);
	if (i != 0)
	{
		perror("stat error");
		LOG.error("Can't determine file type of " + filename, 0);
	}
	if (!S_ISREG(buf.st_mode))
		LOG.error("Does not appear to be a regular file: " + filename, 0);

	if (filename.substr(filename.size()-4) == ".bcf")
		LOG.error("Filename ends in '.bcf'. Shouldn't you be using --bcf?\n");

	if (!compressed)
	{
		if (filename.substr(filename.size()-3) == ".gz")
			LOG.error("Filename ends in '.gz'. Shouldn't you be using --gzvcf or --gzdiff?\n");
		vcf_in.open(filename.c_str(), ios::in);
		if (!vcf_in.is_open())
			LOG.error("Could not open VCF file: " + filename, 0);
	}
	else
	{
		gzMAX_LINE_LEN = 1024*1024;
		gz_readbuffer = new char[gzMAX_LINE_LEN];
		gzvcf_in = gzopen(filename.c_str(), "rb");
		if (gzvcf_in == NULL)
			LOG.error("Could not open GZVCF file: " + filename, 0);
#ifdef ZLIB_VERNUM
		string tmp(ZLIB_VERSION);
		LOG.printLOG("Using zlib version: " + tmp + "\n");
	#if (ZLIB_VERNUM >= 0x1240)
		gzbuffer(gzvcf_in, gzMAX_LINE_LEN); // Included in zlib v1.2.4 and makes things MUCH faster
	#else
		LOG.printLOG("Versions of zlib >= 1.2.4 will be *much* faster when reading zipped VCF files.\n");
	#endif
#endif
	}
}

void vcf_file::close()
{
	if (!compressed)
		vcf_in.close();
	else
	{
		gzclose(gzvcf_in);
		delete [] gz_readbuffer;
	}
}

bool vcf_file::eof()
{
	bool out;
	if (!compressed)
		out = vcf_in.eof();
	else
	{
		out = gzeof(gzvcf_in);	// Returns 1 when EOF has previously been detected reading the given input stream, otherwise zero.
	}
	return out;
}

streampos vcf_file::get_filepos()
{
	if (!compressed)
		return vcf_in.tellg();
	else
	{
		return gztell(gzvcf_in);	// TODO: Type check
	}
}

void vcf_file::set_filepos(streampos &filepos)
{
	if (!compressed)
	{
		vcf_in.clear();
		vcf_in.seekg(filepos, ios::beg);
	}
	else
	{
		gzseek(gzvcf_in, filepos, SEEK_SET);
	}
}

void vcf_file::get_entry(unsigned int entry_num, vector<char> &out)
{
	set_filepos( entry_file_locations[entry_num] );
	read_line(out);
}

entry* vcf_file::get_entry_object(unsigned int N_indv)
{
	return new vcf_entry(N_indv);
}

void vcf_file::read_line(string &out)
{
	if (!compressed)
	{
		getline(vcf_in, out);
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
	else
	{
		out = "";
		bool again = true;
		while (again == true)
		{
			gzgets(gzvcf_in, gz_readbuffer, gzMAX_LINE_LEN);
			out.append(gz_readbuffer);
			if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
				again = false;
		}
		out.erase( out.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)
	}
}

void vcf_file::read_line(vector<char> &out)
{
	static string tmp;
	tmp.resize(0);
	if (!compressed)
	{
		getline(vcf_in, tmp);
		tmp.erase( tmp.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line
	}
	else
	{
		bool again = true;
		while (again == true)
		{
			gzgets(gzvcf_in, gz_readbuffer, gzMAX_LINE_LEN);
			tmp.append(gz_readbuffer);
			if (strlen(gz_readbuffer) != gzMAX_LINE_LEN-1)
				again = false;
		}
		tmp.erase( tmp.find_last_not_of(" \t\n\r") + 1);	// Trim whitespace at end of line (required in gzipped case!)
	}
//	out.assign(tmp.begin(),tmp.end());
	vector<char> tmp_char(tmp.begin(),tmp.end());
	out = tmp_char;
}

char vcf_file::peek()
{
	if (!compressed)
		return vcf_in.peek();
	else
	{
		char c = gzgetc(gzvcf_in);
		gzungetc(c, gzvcf_in);
		return c;
	}
}

int vcf_file::read_CHROM_and_POS_and_skip_remainder_of_line(string &CHROM, int &POS)
{
	if (!compressed)
	{
		getline(vcf_in, CHROM, '\t');
		vcf_in >> POS;
		vcf_in.ignore(std::numeric_limits<streamsize>::max(), '\n');
	}
	else
	{
		static string line;
		static stringstream ss;
		read_line(line);
		ss.clear(); ss.str(line);
		getline(ss, CHROM, '\t');
		ss >> POS;
	}
	return eof();
}

void vcf_file::read_CHROM_only(string &CHROM)
{	// Just read in the chromosome. Note: leaves the stream in a funny state, but is faster than reading whole line
	if (!compressed)
	{
		getline(vcf_in, CHROM, '\t');
	}
	else
	{
		CHROM = "";
		char c = gzgetc(gzvcf_in);
		while (c != '\t')
		{
			CHROM += c;
			c = gzgetc(gzvcf_in);
		}
	}
}

void vcf_file::read_CHROM_and_POS_only(string &CHROM, int &POS)
{	// Just read in the chromosome and position. Note: leaves the stream in a funny state, but is faster than reading whole line
	if (!compressed)
	{
		getline(vcf_in, CHROM, '\t');
		vcf_in >> POS;
	}
	else
	{
		CHROM = "";
		char c = gzgetc(gzvcf_in);
		while (c != '\t')
		{
			CHROM += c;
			c = gzgetc(gzvcf_in);
		}
		string tmp;
		c = gzgetc(gzvcf_in);
		while (c != '\t')
		{
			tmp += c;
			c = gzgetc(gzvcf_in);
		}
		POS = atoi(tmp.c_str());
	}
}
