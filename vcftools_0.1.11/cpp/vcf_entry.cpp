/*
 * vcf_entry.cpp
 *
 *  Created on: Aug 19, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"

map<string, Field_description> vcf_entry::INFO_map;
map<string, unsigned int> vcf_entry::INFO_reverse_map;
map<string, string> vcf_entry::FILTER_map;
map<string, unsigned int> vcf_entry::FILTER_reverse_map;
map<string, Field_description> vcf_entry::FORMAT_map;
map<string, unsigned int> vcf_entry::FORMAT_reverse_map;
map<string, unsigned int> vcf_entry::CONTIG_map;
string vcf_entry::convert_line;

vcf_entry::vcf_entry(const unsigned int n_indv, const vector<char> &line)
{
	N_indv = n_indv;
	basic_parsed = false; fully_parsed = false;
	parsed_ALT = false; parsed_FILTER = false;
	parsed_INFO = false; parsed_FORMAT = false;
	CHROM = ""; POS = -1; REF = ""; QUAL = -1;
	passed_filters = false; parsed_FORMAT_binary = false;
	N_INFO_removed = 0; N_FORMAT_removed = 0;
	parsed_GT = vector<bool>(N_indv, false); parsed_GQ = vector<bool>(N_indv, false);
	parsed_DP = vector<bool>(N_indv, false); parsed_FT = vector<bool>(N_indv, false);
	GT_idx = -1; GQ_idx = -1; DP_idx = -1; FT_idx = -1;
	FORMAT_positions.resize(n_indv); FORMAT_types.resize(n_indv); FORMAT_sizes.resize(n_indv); FORMAT_skip.resize(n_indv); FORMAT_keys.resize(n_indv);
	convert_line.clear();
	convert_line.assign(line.begin(), line.end());
	data_stream.str(convert_line);
}

// Create an empty VCF entry
vcf_entry::vcf_entry(const unsigned int n_indv)
{
	N_indv = n_indv;
	basic_parsed = false; fully_parsed = false;
	parsed_ALT = false; parsed_FILTER = false;
	parsed_INFO = false; parsed_FORMAT = false;
	CHROM = ""; POS = -1; REF = ""; QUAL = -1;
	passed_filters = false; parsed_FORMAT_binary = false;
	N_INFO_removed = 0; N_FORMAT_removed = 0;
	parsed_GT = vector<bool>(N_indv, false); parsed_GQ = vector<bool>(N_indv, false);
	parsed_DP = vector<bool>(N_indv, false); parsed_FT = vector<bool>(N_indv, false);
	GT_idx = -1; GQ_idx = -1; DP_idx = -1; FT_idx = -1;
	FORMAT_positions.resize(n_indv); FORMAT_types.resize(n_indv); FORMAT_sizes.resize(n_indv); FORMAT_skip.resize(n_indv); FORMAT_keys.resize(n_indv);
	convert_line.clear();
	data_stream.str("");
}

vcf_entry::~vcf_entry() {}

// Reset the VCF entry object with a new data line
void vcf_entry::reset(const vector<char> &data_line)
{
	basic_parsed = false;
	fully_parsed = false;
	parsed_ALT = false;
	parsed_FILTER = false;
	parsed_INFO = false;
	parsed_FORMAT = false;
	parsed_FORMAT_binary = false;

	data_stream.clear();
	convert_line.clear();
	convert_line.assign(data_line.begin(), data_line.end());
	data_stream.str(convert_line);

	fill(parsed_GT.begin(), parsed_GT.end(), false);
	fill(parsed_GQ.begin(), parsed_GQ.end(), false);
	fill(parsed_DP.begin(), parsed_DP.end(), false);
	fill(parsed_FT.begin(), parsed_FT.end(), false);

	N_INFO_removed = 0; N_FORMAT_removed = 0;
	FORMAT_positions.clear(); FORMAT_types.clear(); FORMAT_sizes.clear(); FORMAT_skip.clear(); FORMAT_keys.clear();
}

// Tokenize the basic information in a VCF data line (at the tab level)
void vcf_entry::parse_basic_entry(bool parse_ALT, bool parse_FILTER, bool parse_INFO)
{
	// The following would break on spaces too, which caused a bug :-(
	//data_stream >> CHROM >> POS >> ID >> REF >> ALT_str >> QUAL_str >> FILTER_str >> INFO_str;

	getline(data_stream, CHROM, '\t');

	getline(data_stream, ID, '\t');
	POS = atoi(ID.c_str());
	getline(data_stream, ID, '\t');
	getline(data_stream, REF, '\t');
	getline(data_stream, ALT_str, '\t');
	getline(data_stream, QUAL_str, '\t');
	getline(data_stream, FILTER_str, '\t');
	getline(data_stream, INFO_str, '\t');

	QUAL = str2double(QUAL_str);

	// Convert to uppercase for consistency
	// Note that VCF v4.1 allows mixtures of lower/upper case in REF and ALT.
	// However, the spec specifically states that tools using VCF are not required
	// to preserve the case.
	std::transform(REF.begin(), REF.end(), REF.begin(), ::toupper);
	std::transform(ALT_str.begin(), ALT_str.end(),ALT_str.begin(), ::toupper);

	parsed_ALT = false;
	parsed_FILTER = false;
	parsed_INFO = false;
	basic_parsed = true;

	if (parse_ALT)
		set_ALT(ALT_str);
	if (parse_FILTER)
		set_FILTER(FILTER_str);
	if (parse_INFO)
		set_INFO(INFO_str);
}

// Tokenize the genotype information (at the 'tab' level) in the VCF entry
void vcf_entry::parse_full_entry(bool parse_FORMAT)
{
	if (basic_parsed == false)
		parse_basic_entry();

	//data_stream >> FORMAT_str;
	getline(data_stream, FORMAT_str, '\t');

	if (parse_FORMAT)
		set_FORMAT(FORMAT_str);

	string tmpstr; tmpstr.reserve(64);
	GENOTYPE_str.resize(N_indv, tmpstr);

	for (unsigned int ui=0; ui<N_indv; ui++)
		//data_stream >> GENOTYPE_str[ui];
		getline(data_stream, GENOTYPE_str[ui], '\t');

	// The following line copies the GENOTYPE fields from the stringstream into the GENOTYPE_str vector.
	// Is actually slower than the above code.
	//copy(istream_iterator<string>(data_stream), istream_iterator<string>(), GENOTYPE_str.begin());

	fully_parsed = true;
}

// Tokenize a given genotype entry into it's component parts
void vcf_entry::parse_genotype_entry(unsigned int indv, bool GT, bool GQ, bool DP, bool FT)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT(FORMAT_str);

	static string tmpstr;
	static istringstream ss;
	ss.clear(); ss.str(GENOTYPE_str[indv]);

	int N_required = GT + GQ + DP + FT;
	int N_got = 0;

	int i=0;
	while (getline(ss, tmpstr, ':'))
	{
		if (GT && (i == GT_idx)) // (FORMAT[ui] == "GT")
		{
			set_indv_GENOTYPE_and_PHASE(indv, tmpstr);
			N_got++;
		}
		else if (GQ && (i == GQ_idx)) // (FORMAT[ui] == "GQ")
		{
			set_indv_GQUALITY(indv, str2double(tmpstr));
			N_got++;
		}
		else if (DP && (i == DP_idx)) // (FORMAT[ui] == "DP")
		{
			set_indv_DEPTH(indv, str2int(tmpstr));
			N_got++;
		}
		else if (FT && (i == FT_idx)) // (FORMAT[ui] == "FT")
		{
			set_indv_GFILTER(indv, tmpstr);
			N_got++;
		}

		if (N_got == N_required)
			break;
		i++;
	}

	// Set missing return values if requested a value, but couldn't find it
	if (GT && (parsed_GT[indv] == false))
	{
		set_indv_GENOTYPE_and_PHASE(indv, make_pair(-1,-1), '/');
	}
	if (GQ && (parsed_GQ[indv] == false))
	{
		set_indv_GQUALITY(indv, -1);
	}
	if (DP && (parsed_DP[indv] == false))
	{
		set_indv_DEPTH(indv, -1);
	}
	if (FT && (parsed_FT[indv] == false))
	{
		set_indv_GFILTER(indv, "");
	}
}

// Read the VCF entry and fully populate the object
void vcf_entry::parse_genotype_entries(bool GT, bool GQ, bool DP, bool FT)
{
	for (unsigned int ui=0; ui<N_indv; ui++)
		parse_genotype_entry(ui, GT, GQ, DP, FT);
}

void vcf_entry::parse_FORMAT()
{
	FORMAT_binary.resize(0);
	vector<char> tmp_vector;
	vector<string> tmp_split;
	vector< vector<string> > format_matrix(N_indv);
	unsigned int type, number, size, position=0;

	for (unsigned int ui=0; ui<N_indv; ui++)
	{
		if ( (GENOTYPE_str[ui] != "") and (GENOTYPE_str[ui] != ".") )
			tokenize(GENOTYPE_str[ui],':', tmp_split);
		else
			tmp_split.assign(FORMAT.size(),".");
		format_matrix[ui] = tmp_split;

		for (unsigned int uj=format_matrix[ui].size(); uj<FORMAT.size(); uj++)
			format_matrix[ui].push_back(".");
	}

	FORMAT_positions.resize(FORMAT.size()); FORMAT_types.resize(FORMAT.size()); FORMAT_skip.resize(FORMAT.size());

	for (unsigned int ui=0; ui<FORMAT.size(); ui++)
	{
		if (FORMAT_reverse_map.find(FORMAT[ui]) == FORMAT_reverse_map.end())
		{
			LOG.one_off_warning("FORMAT value " + FORMAT[ui] + " is not defined in the header and will not be encoded.");
			N_FORMAT_removed++;
			continue;
		}
		FORMAT_positions[ui] = position;
		make_typed_int(tmp_vector, FORMAT_reverse_map[ FORMAT[ui] ], true );
		FORMAT_binary.insert(FORMAT_binary.end(), tmp_vector.begin(), tmp_vector.end() );

		tmp_vector.resize(0);
		tmp_split.resize(0);

		type = FORMAT_map[ FORMAT[ui] ].Type;
		number = FORMAT_map[ FORMAT[ui] ].N_entries;

		for (unsigned int uj=0; uj<N_indv; uj++)
			tmp_split.push_back( format_matrix[uj][ui] );

		if ((int)ui == GT_idx)
			make_typed_GT_vector(tmp_vector, tmp_split );
		else if (type == Integer)
			make_typed_int_vector(tmp_vector, tmp_split, number );
		else if (type == Float)
			make_typed_float_vector(tmp_vector, tmp_split, number );
		else if ( (type == Character) or (type == String) )
			make_typed_string_vector(tmp_vector, tmp_split, number );
		else
			LOG.error("Invalid type in FORMAT definition", 0);

		position = 0;
		get_type(&position, tmp_vector, type, size);
		FORMAT_types[ui] = type;

		if ( (type == 1) || (type == 7) )
			FORMAT_skip[ui] = size*sizeof(int8_t);
		else if (type == 2)
			FORMAT_skip[ui] = size*sizeof(int16_t);
		else if ( (type == 3) || (type == 5) )
			FORMAT_skip[ui] = size*sizeof(int32_t);

		FORMAT_binary.insert(FORMAT_binary.end(), tmp_vector.begin(), tmp_vector.end() );
		tmp_vector.resize(0);
		position = FORMAT_binary.size();
	}
	parsed_FORMAT_binary = true;
}

void vcf_entry::print(ostream &out)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	set<string> INFO_to_keep;
	print(out, INFO_to_keep, false, include_indv, include_genotype);
}

void vcf_entry::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	print(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype);
}

// Output VCF entry to output stream
void vcf_entry::print(ostream &out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype)
{
	if (fully_parsed == false)
		parse_full_entry();

	out << get_CHROM() << '\t' << POS << '\t' << get_ID() << '\t' << REF << '\t' << get_ALT();

	out << '\t' << double2str(QUAL);
	out << '\t' << get_FILTER();
	if (keep_all_INFO == false)
		out << '\t' << get_INFO(INFO_to_keep);
	else
		out << '\t' << INFO_str;

	pair<int, int> genotype;
	string GFILTER_tmp;
	if (FORMAT.size() > 0)
	{
		char PHASE;
		out << '\t' << get_FORMAT();

		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (include_indv[ui] == false)
				continue;
			out << '\t';
			for (int count=0; count<(int)FORMAT.size(); count++)
			{
				if (count == GT_idx) // (FORMAT[count] == "GT")
				{
					if (count != 0)	out << ':';
					if (include_genotype[ui] == true)
					{
						get_indv_GENOTYPE_ids(ui, genotype);
						PHASE = get_indv_PHASE(ui);
						if ((genotype.first != -1) && (genotype.second != -1))
							out << int2str(genotype.first) << PHASE << int2str(genotype.second);
						else if ((PHASE == '|') && (genotype.second == -1))
							out << int2str(genotype.first);	// Handle haploid case
						else
							out << int2str(genotype.first) << PHASE << int2str(genotype.second);
					}
					else
						out << "./.";
				}
				else if (count == GQ_idx) //(FORMAT[count] == "GQ")
				{
					if (count != 0)	out << ':';
					out << double2str(get_indv_GQUALITY(ui));
				}
				else if (count == DP_idx) // (FORMAT[count] == "DP")
				{
					if (count != 0)	out << ':';
					out << int2str(get_indv_DEPTH(ui));
				}
				else if (count == FT_idx) // (FORMAT[count] == "FT")
				{
					if (count != 0)	out << ':';
					get_indv_GFILTER(ui, GFILTER_tmp);
					out << GFILTER_tmp;
				}
				else
				{	// Unknown FORMAT so just replicate original output
					if (count != 0)	out << ':';
					read_indv_generic_entry(ui, FORMAT[count], GFILTER_tmp);
					out << GFILTER_tmp;
				}
			}
		}
	}
	out << '\n';	// endl flushes the buffer, which is slow. This (should be) quicker.
}

void vcf_entry::print_bcf(BGZF* out)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	set<string> INFO_to_keep;
	print_bcf(out, INFO_to_keep, false, include_indv, include_genotype);
}

void vcf_entry::print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO)
{
	vector<bool> include_indv(N_indv, true);
	vector<bool> include_genotype(N_indv, true);
	print_bcf(out, INFO_to_keep, keep_all_INFO, include_indv, include_genotype);
}

// Output VCF entry to output stream in binary
void vcf_entry::print_bcf(BGZF* out, const set<string> &INFO_to_keep, bool keep_all_INFO, const vector<bool> &include_indv, const vector<bool> &include_genotype)
{
	if (fully_parsed == false)
		parse_full_entry();

	if (parsed_FORMAT_binary == false)
		parse_FORMAT();

	vector<char> out_vector, tmp_vector;
	out_vector.resize(8*sizeof(int32_t));
	int vector_pos = 2*sizeof(uint32_t);
	string tmp_string;
	int index;
	vector<string> filter_vector;
	vector<pair< string, string > > tmp_info;

	tmp_string = get_CHROM();
	if (tmp_string == "." or tmp_string == " " or tmp_string == "")
		LOG.error("CHROM value must be defined for all entries.",0);

	if (CONTIG_map.find(tmp_string) == CONTIG_map.end() )
		LOG.error("CHROM value " + tmp_string + " is not defined on contig dictionary.",0);

	int32_t chrom = (int32_t)CONTIG_map[tmp_string];
	memcpy(&out_vector[vector_pos], &chrom, sizeof(chrom));
	vector_pos += sizeof(chrom);

	get_POS_binary(tmp_vector);
	memcpy(&out_vector[vector_pos], &tmp_vector[0], tmp_vector.size());
	vector_pos += tmp_vector.size();
	tmp_vector.resize(0);

	get_rlen(tmp_vector);
	memcpy(&out_vector[vector_pos], &tmp_vector[0], tmp_vector.size());
	vector_pos += tmp_vector.size();
	tmp_vector.resize(0);

	get_QUAL_binary(tmp_vector);
	memcpy(&out_vector[vector_pos], &tmp_vector[0], tmp_vector.size());
	vector_pos += tmp_vector.size();
	tmp_vector.resize(0);

	get_ID_binary(tmp_vector);
	out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());
	tmp_vector.resize(0);

	get_ALLELES_binary(tmp_vector);
	out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());
	tmp_vector.resize(0);

	get_FILTER_vector(filter_vector);
	if (passed_filters == true)
		make_typed_int(tmp_vector, 0, true);
	else if (filter_vector.empty())
		make_typed_int_vector(tmp_vector, filter_vector);
	else
	{
		vector<int> index_vector;
		for(unsigned int ui=0; ui<filter_vector.size(); ui++)
		{
			if ( FILTER_reverse_map.find( filter_vector[ui] ) == FILTER_reverse_map.end() )
				LOG.one_off_warning("FILTER value " + filter_vector[ui] + " is not defined in the header and will not be encoded.");
			else
				index_vector.push_back( FILTER_reverse_map[ filter_vector[ui] ] );
		}
		make_typed_int_vector(tmp_vector, index_vector );
	}

	if (tmp_vector.empty())
		make_typed_int(tmp_vector, 0, true);

	out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());

	int map_type, number;
	tmp_info = get_INFO_vector(INFO_to_keep, keep_all_INFO);

	for(unsigned int ui=0; ui<tmp_info.size(); ui++)
	{
		if (INFO_reverse_map.find(tmp_info[ui].first) == INFO_reverse_map.end())
		{
			LOG.one_off_warning("INFO value " + tmp_info[ui].first + " is not defined in the header and will not be encoded.");
			N_INFO_removed++;
			continue;
		}

		tmp_vector.resize(0);
		index = INFO_reverse_map[ tmp_info[ui].first ];
		make_typed_int(tmp_vector, index, true);
		out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());

		tmp_vector.resize(0);
		map_type = INFO_map[ tmp_info[ui].first ].Type;
		number = INFO_map[ tmp_info[ui].first ].N_entries;

		if (map_type == Integer)
			make_typed_int_vector(tmp_vector, tmp_info[ui].second, number );
		else if (map_type == Float)
			make_typed_float_vector(tmp_vector, tmp_info[ui].second, number );
		else if ( (map_type == Character) or (map_type == String) )
			make_typed_string(tmp_vector, tmp_info[ui].second, true );
		else if (map_type == Flag)
			make_typed_int(tmp_vector, 1, true );
		else
			LOG.error("Invalid type in INFO definition", 0);

		out_vector.insert(out_vector.end(), tmp_vector.begin(), tmp_vector.end());
	}
	tmp_vector.resize(0);
	uint32_t l_shared = (uint32_t)out_vector.size() - (uint32_t)(2*sizeof(uint32_t));
	memcpy(&out_vector[0], &l_shared, sizeof(l_shared));

	tmp_vector.resize(0);
	get_n_allele_info(tmp_vector);
	memcpy(&out_vector[vector_pos], &tmp_vector[0], tmp_vector.size());
	vector_pos += tmp_vector.size();

	tmp_vector.resize(0);
	get_FORMAT_binary(tmp_vector);
	unsigned int key_size, size, type, old_pos, out_size, format_pos = 0;
	unsigned int skip_size = 0;
	for (unsigned int ui=0; ui<(FORMAT.size()-N_FORMAT_removed); ui++)
	{
		format_pos = FORMAT_positions[ui];
		old_pos = format_pos;

		get_typed_int(&format_pos, tmp_vector, type, size);
		get_type(&format_pos, tmp_vector, type, size);
		key_size = format_pos - old_pos;
		out_size = out_vector.size();
		out_vector.resize( out_size + key_size );
		memcpy(&out_vector[out_size], &tmp_vector[old_pos], key_size );

		skip_size = FORMAT_skip[ui];
		for (unsigned int uj=0; uj<N_indv; uj++)
		{
			if (include_indv[uj] == false)
				continue;

			format_pos = FORMAT_positions[ui]+key_size+skip_size*uj;
			out_size = out_vector.size();
			out_vector.resize( out_size + skip_size );
			if ( ((int)ui == GT_idx) and (include_genotype[ui] == false) )
			{
				for (int p = 0; p < (int)skip_size; p++)
				{
					if (p < ploidy[uj])
						out_vector[out_size] = (int8_t)0x00;
					else
						out_vector[out_size] = (int8_t)0x80;
					out_size++;
				}
			}
			else
			{
				memcpy(&out_vector[out_size], &tmp_vector[format_pos], skip_size);
				out_size += skip_size;
			}
		}
	}
	tmp_vector.resize(0);
	get_n_fmt_sample(tmp_vector);
	memcpy(&out_vector[vector_pos], &tmp_vector[0], tmp_vector.size());

	uint32_t l_indv = (uint32_t)out_vector.size() - l_shared - (uint32_t)(2*sizeof(uint32_t));
	memcpy(&out_vector[sizeof(l_shared)], &l_indv, sizeof(l_indv));
	bgzf_write(out, &out_vector[0], out_vector.size());
}

// Set the include_genotype flag on the basis of depth
void vcf_entry::filter_genotypes_by_depth(vector<bool> &include_genotype_out, int min_depth, int max_depth)
{
	if (fully_parsed == false)
		parse_full_entry();

	//if (FORMAT_to_idx.find("DP") != FORMAT_to_idx.end())
	if (DP_idx != -1)
	{	// Have depth info
		int depth;
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_DP[ui] == false)
				parse_genotype_entry(ui, false, false, true);
			depth = get_indv_DEPTH(ui);
			if ((depth < min_depth) || (depth > max_depth))
				include_genotype_out[ui] = false;
		}
	}
}

// Filter specific genotypes by quality
void vcf_entry::filter_genotypes_by_quality(vector<bool> &include_genotype_out, double min_genotype_quality)
{
	if (fully_parsed == false)
		parse_full_entry();

	//if (FORMAT_to_idx.find("GQ") != FORMAT_to_idx.end())
	if (GQ_idx != -1)
	{	// Have quality info
		double quality;
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_GQ[ui] == false)
				parse_genotype_entry(ui, false, true);
			quality = get_indv_GQUALITY(ui);
			if (quality < min_genotype_quality)
				include_genotype_out[ui] = false;
		}
	}
}

// Exclude genotypes with a filter flag.
void vcf_entry::filter_genotypes_by_filter_status(vector<bool> &include_genotype_out, const set<string> &filter_flags_to_remove, bool remove_all)
{
	if (fully_parsed == false)
		parse_full_entry();

	vector<string> GFILTERs;
	if (FT_idx != -1)
	{	// Have GFilter info
		include_genotype_out.resize(N_indv, true);
		for (unsigned int ui=0; ui<N_indv; ui++)
		{
			if (parsed_FT[ui] == false)
				parse_genotype_entry(ui, false, false, false, true);
			get_indv_GFILTER_vector(ui, GFILTERs);

			if (remove_all == true)
			{	// If removing all filters, only keep things with label PASS
				if (!GFILTERs.empty())
					if ((GFILTERs[0] != "PASS") && (GFILTERs[0] != "."))
						include_genotype_out[ui] = false;
			}
			else
			{	// Only removing specific filters
				for (unsigned int uj=0; uj<GFILTERs.size(); uj++)
					if (filter_flags_to_remove.find(GFILTERs[uj]) != filter_flags_to_remove.end())
							include_genotype_out[ui] = false;
			}
		}
	}
}

void vcf_entry::read_indv_generic_entry(unsigned int indv, const string &FORMAT_id, string &out)
{
	if (fully_parsed == false)
		parse_full_entry(true);

	if (parsed_FORMAT == false)
		set_FORMAT(FORMAT_str);

	out = ".";

	if (FORMAT_to_idx.find(FORMAT_id) != FORMAT_to_idx.end())
	{
		unsigned int idx = FORMAT_to_idx[FORMAT_id];
		static string tmpstr;
		static istringstream ss;
		ss.clear();
		ss.str(GENOTYPE_str[indv]);

		for (unsigned int ui=0; ui <= idx; ui++)
		{
			getline(ss, tmpstr, ':');
			if (ui == idx)
			{
				out = tmpstr;
				break;
			}
			if (!ss.good())
				break;
		}
	}
}
