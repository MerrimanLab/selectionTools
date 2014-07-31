/*
 * vcf_entry_setters.cpp
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 230 $)
 */

#include "vcf_entry.h"
#include "entry.h"

void vcf_entry::set_ALT(const string &in)
{
	istringstream ss(in);
	string tmpstr;
	ALT.resize(0);
	while(!ss.eof())
	{
		getline(ss, tmpstr, ',');
		add_ALT_allele(tmpstr);
	}
	parsed_ALT = true;
}

void vcf_entry::set_QUAL(const double in)
{
	QUAL = in;
}

void vcf_entry::set_FORMAT(const string &in)
{
	FORMAT.resize(0);
	FORMAT_to_idx.clear();

	if (in.size() > 0)
	{
		istringstream ss(in);
		string tmpstr;

		unsigned int pos=0;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ':');
			add_FORMAT_entry(tmpstr, pos);
			pos++;
		}
	}

	GT_idx = -1;
	GQ_idx = -1;
	DP_idx = -1;
	FT_idx = -1;

	if (FORMAT_to_idx.find("GT") != FORMAT_to_idx.end())
		GT_idx = FORMAT_to_idx["GT"];
	if (FORMAT_to_idx.find("GQ") != FORMAT_to_idx.end())
		GQ_idx = FORMAT_to_idx["GQ"];
	if (FORMAT_to_idx.find("DP") != FORMAT_to_idx.end())
		DP_idx = FORMAT_to_idx["DP"];
	if (FORMAT_to_idx.find("FT") != FORMAT_to_idx.end())
		FT_idx = FORMAT_to_idx["FT"];

	parsed_FORMAT = true;
}

void vcf_entry::add_FORMAT_entry(const string &in, unsigned int pos)
{
	FORMAT.push_back(in);
	FORMAT_to_idx[in] = pos;
}

// The following function reads in a genotype from a '0/1'-like string.
// Should handle haploid types to, but NOT polyploidy.
void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const string &in)
{
	ploidy.resize(N_indv);
	if ((in.size() == 3) && ((in.c_str()[1] == '/') || (in.c_str()[1] == '|')))
	{	// Fast, diploid case...
		ploidy[indv] = 2;
		set_indv_PHASE(indv, in.c_str()[1]);
		set_indv_GENOTYPE_alleles(indv, in.c_str()[0], in.c_str()[2]);
	}
	else
	{	// More complex case...
		size_t pos = in.find_first_of("/|");
		if (pos != string::npos)
		{	// autosome
			ploidy[indv] = 2;
			set_indv_PHASE(indv, in[pos]);
			set_indv_GENOTYPE_alleles(indv, make_pair(in.substr(0,pos), in.substr(pos+1)));
		}
		else
		{	// Male chrX, or chrY
			ploidy[indv] = 1;
			set_indv_PHASE(indv, '|');
			set_indv_GENOTYPE_alleles(indv, make_pair(in.substr(0,pos), "."));
		}

		// Check for polypoidy
		size_t pos2 = in.find_last_of("/|");
		if (pos != pos2)
			LOG.error("Polypolidy found, and not supported by vcftools: " + CHROM + ":" + int2str(POS));
	}

	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<int, int> &genotype, char phase)
{
	ploidy.resize(N_indv);
	ploidy[indv] = 2;
	set_indv_GENOTYPE_ids(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_and_PHASE(unsigned int indv, const pair<string, string> &genotype, char phase)
{
	ploidy.resize(N_indv);
	ploidy[indv] = 2;
	set_indv_GENOTYPE_alleles(indv, genotype);
	set_indv_PHASE(indv, phase);
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, const pair<string, string> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (in.first != ".")
		a.first = str2int(in.first);

	if (in.second != ".")
		a.second = str2int(in.second);

	GENOTYPE[indv] = a;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GENOTYPE_alleles(unsigned int indv, char a1, char a2)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));

	pair<int, int> a(-1,-1);

	if (a1 != '.')
		a.first = a1 - '0';

	if (a2 != '.')
		a.second = a2 - '0';

	GENOTYPE[indv] = a;
	parsed_GT[indv] = true;
}


void vcf_entry::set_indv_GENOTYPE_ids(unsigned int indv, const pair<int, int> &in)
{
	if (GENOTYPE.size() == 0)
		GENOTYPE.resize(N_indv, make_pair(-1,-1));
	GENOTYPE[indv] = in;
}

void vcf_entry::set_indv_PHASE(unsigned int indv, char in)
{
	if (PHASE.size() == 0)
		PHASE.resize(N_indv, '/');

	PHASE[indv] = in;
	parsed_GT[indv] = true;
}

void vcf_entry::set_indv_GQUALITY(unsigned int indv, double in)
{
	parsed_GQ[indv] = true;
	if (in == -1)
	{
		if (GQUALITY.size() > 0)
			GQUALITY[indv] = -1;
		return;
	}
	if (GQUALITY.size() == 0)
		GQUALITY.resize(N_indv, -1);

	if (in > 99)
		in = 99;
	GQUALITY[indv] = in;
}

void vcf_entry::set_indv_DEPTH(unsigned int indv, int in)
{
	parsed_DP[indv] = true;
	if (in == -1)
	{
		if (DEPTH.size() > 0)
			DEPTH[indv] = -1;
		return;
	}
	if (DEPTH.size() == 0)
		DEPTH.resize(N_indv, -1);

	DEPTH[indv] = in;
}

void vcf_entry::add_indv_GFILTER(unsigned int indv, const string &in)
{
	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	if (in != ".")
		if (find(GFILTER[indv].begin(), GFILTER[indv].end(), in) == GFILTER[indv].end())
			GFILTER[indv].push_back(in);
	parsed_FT[indv] = true;
}

void vcf_entry::set_indv_GFILTER(unsigned int indv, const string &in)
{
	parsed_FT[indv] = true;

	if (GFILTER.size() == 0)
		GFILTER.resize(N_indv);

	GFILTER[indv].resize(0);
	if ((in.size() == 0) || (in == "."))
		return;

	static istringstream ss;
	static string ith_FILTER;
	ss.clear();
	ss.str(in);
	while (!ss.eof())
	{
		getline(ss, ith_FILTER, ';');

		if ((ith_FILTER.size()==0) || (ith_FILTER == "."))
			continue;	// Don't bother storing "unfiltered" state.

		GFILTER[indv].push_back(ith_FILTER);
	}
}

void vcf_entry::set_FILTER(const string &FILTER_str)
{
	FILTER.resize(0);
	passed_filters = false;
	if (FILTER_str == "PASS")
		passed_filters = true;
	else
	{
		if (FILTER_str != ".")
		{
			istringstream ss(FILTER_str);
			string ith_FILTER;
			while (!ss.eof())
			{
				getline(ss, ith_FILTER, ';');
				FILTER.push_back(ith_FILTER);
			}
		}
	}
	sort(FILTER.begin(), FILTER.end());
	parsed_FILTER = true;
}

void vcf_entry::set_INFO(const string &INFO_str)
{
	INFO.resize(0);
	if ((INFO_str.size() > 0) && (INFO_str != "."))
	{
		istringstream ss(INFO_str);
		string tmpstr;
		while(!ss.eof())
		{
			getline(ss, tmpstr, ';');

			istringstream ss2(tmpstr);
			getline(ss2, tmpstr, '=');
			pair<string, string> INFO_entry(tmpstr, ".");

			if (!ss2.eof())
			{	// If there is a value entry, read it now
				getline(ss2, tmpstr);
				INFO_entry.second = tmpstr;
			}
			else	// Otherwise, set it equal to 1
				INFO_entry.second = "1";

			INFO.push_back(INFO_entry);
		}
	}
	parsed_INFO = true;
}

int vcf_entry::add_INFO_descriptor(const string &in, unsigned int index)
{
	size_t found=in.find("##INFO=");
	if (found!=string::npos)
	{	// Found an INFO descriptor
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		Field_description I;

		vector<string> tokens;
		tokenize(details, ',', tokens);
		if (tokens.size() < 4)
			LOG.error("Expected 4 parts in INFO definition: " + in);

		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else LOG.error("Expected ID entry as first field in INFO description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Number")
		{	// TODO - handle 'A' and 'G' categories correctly.
			if ((entry[1] == "A") || (entry[1] == "G"))
				I.N_entries = -1;	// Currently just treat as missing.
			else
				I.N_entries =  str2int(entry[1]);
		}
		else LOG.error("Expected Number entry as second field in INFO description: " + in);

		tokenize(tokens[2], '=', entry);
		if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") I.Type = Integer;
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) I.Type = Float;
			else if (entry[1] == "Character") I.Type = Character;
			else if (entry[1] == "String") I.Type = String;
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
			}
			else LOG.error("Unknown Type in INFO meta-information: " + in);
		}
		else LOG.error("Expected Type entry as third field in INFO description: " + in);

		tokenize(tokens[3], '=', entry);
		if (entry[0] == "Description")
		{
			I.Description = entry[1];
			for (unsigned int i=4; i<tokens.size(); i++)
			{
				I.Description += "; " + tokens[i];
			}
		}
		else LOG.error("Expected Description entry as fourth field in INFO description: " + in);

		if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
		{
			INFO_map[ I.ID ] = I;
			INFO_reverse_map[I.ID] = FORMAT_reverse_map[I.ID];
			return 0;
		}
		else if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
		{
			INFO_map[ I.ID ] = I;
			INFO_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
			return 0;
		}
		else if (INFO_map.find(I.ID) != INFO_map.end())
			return 0;
		else
		{
			INFO_map[I.ID] = I;
			INFO_reverse_map[I.ID] = index;
			return 1;
		}
	}
	return 0;
}

int vcf_entry::add_FILTER_descriptor(const string &in, unsigned int index)
{
	size_t found=in.find("##FILTER=");
	if (found!=string::npos)
	{
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		vector<string> tokens;
		tokenize(details, ',', tokens);

		if (tokens.size() < 2)
			LOG.error("Expected 2 parts in FILTER definition: " + in);

		string ID, Description;
		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") ID = entry[1];
		else LOG.error("Expected ID as first field in FILTER description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Description")
		{
			Description = entry[1];
			for (unsigned int i=2; i<tokens.size(); i++)
			{
				Description += "; " + tokens[i];
			}
		}
		else
			LOG.error("Expected Description as second field in FILTER description: " + in);

		if ( INFO_reverse_map.find( ID ) != INFO_reverse_map.end() )
		{
			FILTER_map[ ID ] = Description;
			FILTER_reverse_map[ID] = INFO_reverse_map[ ID ];
			return 0;
		}
		else if ( FORMAT_reverse_map.find( ID ) != FORMAT_reverse_map.end() )
		{
			FILTER_map[ ID ] = Description;
			FILTER_reverse_map[ID] = FORMAT_reverse_map[ID];
			return 0;
		}
		else if (FILTER_map.find(ID) != FILTER_map.end())
			return 0;
		else
		{
			FILTER_map[ID] = Description;
			FILTER_reverse_map[ID] = index;
			return 1;
		}
	}
	return 0;
}

int vcf_entry::add_FORMAT_descriptor(const string &in, unsigned int index)
{
	size_t found=in.find("##FORMAT=");
	if (found!=string::npos)
	{	// Found an FORMAT descriptor
		size_t found_start=in.find_first_of("<");
		size_t found_end=in.find_last_of(">");
		string details = in.substr(found_start+1, found_end-found_start-1);
		vector<string> tokens;
		tokenize(details, ',', tokens);
		Field_description I;

		if (tokens.size() < 4)
			LOG.error("Expected 4 parts in FORMAT definition: " + in);

		vector<string> entry;
		tokenize(tokens[0], '=', entry);
		if (entry[0] == "ID") I.ID = entry[1];
		else LOG.error("Expected ID entry as first field in FORMAT description: " + in);

		tokenize(tokens[1], '=', entry);
		if (entry[0] == "Number")
		{	// TODO - handle 'A' and 'G' categories correctly.
			if ((entry[1] == "A") || (entry[1] == "G"))
				I.N_entries = -1;	// Currently just treat as missing.
			else
				I.N_entries = str2int(entry[1]);
		}
		else LOG.error("Expected Number entry as second field in FORMAT description: " + in);

		tokenize(tokens[2], '=', entry);
		if (entry[0] == "Type")
		{
			if (entry[1] == "Integer") I.Type = Integer;
			else if ((entry[1] == "Float") || (entry[1] == "Numeric")) I.Type = Float;
			else if (entry[1] == "Character") I.Type = Character;
			else if (entry[1] == "String") I.Type = String;
			else if (entry[1] == "Flag")
			{
				I.Type = Flag;
				if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
			}
			else LOG.error("Unknown Type in FORMAT meta-information: " + in);
		}
		else LOG.error("Expected Type entry as third field in FORMAT description: " + in);

		tokenize(tokens[3], '=', entry);
		if (entry[0] == "Description")
		{
			I.Description = entry[1];
			for (unsigned int i=4; i<tokens.size(); i++)
			{
				I.Description += "; " + tokens[i];
			}
		}
		else LOG.error("Expected Description entry as fourth field in FORMAT description: " + in);

		if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
		{
			FORMAT_map[ I.ID ] = I;
			FORMAT_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
			return 0;
		}
		else if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
		{
			FORMAT_map[ I.ID ] = I;
			FORMAT_reverse_map[I.ID] = INFO_reverse_map[ I.ID ];
			return 0;
		}
		else if (FORMAT_map.find(I.ID) != FORMAT_map.end())
			return 0;
		else
		{
			FORMAT_map[I.ID] = I;
			FORMAT_reverse_map[I.ID] = index;
			return 1;
		}
	}
	return 0;
}

void vcf_entry::add_CONTIG_descriptor(const string &in, unsigned int index)
{
	size_t found_start=in.find_first_of("<");
	size_t found_end=in.find_last_of(">");
	string details = in.substr(found_start+1, found_end-found_start-1);

	vector<string> tokens;
	entry::tokenize(details, ',', tokens);
	Field_description I;
	bool id_found = false;
	vector<string> entry;

	for (unsigned int ui=0; ui<tokens.size(); ui++)
	{
		entry::tokenize(tokens[ui], '=', entry);
		if (entry[0] == "ID")
		{
			I.ID = entry[1];
			id_found = true;
		}
		else if (entry[0] == "length") I.Length = entry[1];
		else if (entry[0] == "assembly") I.Assembly = entry[1];
	}
	if (id_found == false)
		LOG.warning("CONTIG declaration found without ID: "+ in + "\n");

	CONTIG_map[I.ID] = index;
}
