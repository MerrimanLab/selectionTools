/*
 * header.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: amarcketta
 */

#include "header.h"
#include "entry.h"

int header::add_INFO_descriptor(const string &in, int index)
{
	Field_description I;
	vector<string> tokens;
	entry::tokenize(in, ',', tokens);

	if (tokens.size() < 4)
		LOG.error("Expected 4 parts in INFO definition: " + in);

	vector<string> entry;
	entry::tokenize(tokens[0], '=', entry);
	if (entry[0] == "ID") I.ID = entry[1];
	else LOG.error("Expected ID entry as first field in INFO description: " + in);

	entry::tokenize(tokens[1], '=', entry);
	if (entry[0] == "Number")
	{
		if ((entry[1] == "A") || (entry[1] == "G"))
		{
			I.N_entries = -1;
			I.N_entries_str = entry[1];
		}
		else{
			I.N_entries =  entry::str2int(entry[1]);
			I.N_entries_str = entry[1];
		}
	}
	else LOG.error("Expected Number entry as second field in INFO description: " + in);

	entry::tokenize(tokens[2], '=', entry);
	if (entry[0] == "Type")
	{
		if (entry[1] == "Integer") { I.Type_str = "Integer"; I.Type = Integer; }
		else if ((entry[1] == "Float") || (entry[1] == "Numeric")) {I.Type_str = "Float"; I.Type = Float;}
		else if (entry[1] == "Character") {I.Type_str = "Character"; I.Type = Character;}
		else if (entry[1] == "String") {I.Type_str = "String"; I.Type = String;}
		else if (entry[1] == "Flag")
		{
			I.Type = Flag;
			I.Type_str = "Flag";
			if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
		}
			else LOG.error("Unknown Type in INFO meta-information: " + in);
	}
		else LOG.error("Expected Type entry as third field in INFO description: " + in);

	entry::tokenize(tokens[3], '=', entry);
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
		INFO_map[ FORMAT_reverse_map[ I.ID ] ] = I;
		INFO_reverse_map[I.ID] = FORMAT_reverse_map[I.ID];
		return 0;
	}
	else if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
	{
		INFO_map[ FILTER_reverse_map[ I.ID ] ] = I;
		INFO_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
		return 0;
	}
	else if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
		return 0;
	else
	{
		INFO_map[index] = I;
		INFO_reverse_map[I.ID] = index;
		return 1;
	}
}

int header::add_FORMAT_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end-1);

	vector<string> tokens;
	entry::tokenize(details, ',', tokens);
	Field_description I;
	if (tokens.size() < 4)
		LOG.error("Expected 4 parts in FORMAT definition: " + in);

	vector<string> entry;
	entry::tokenize(tokens[0], '=', entry);
	if (entry[0] == "ID") I.ID = entry[1];
	else LOG.error("Expected ID entry as first field in FORMAT description: " + in);

	entry::tokenize(tokens[1], '=', entry);
	if (entry[0] == "Number")
	{
		if ((entry[1] == "A") || (entry[1] == "G"))
			I.N_entries = -1;
		else
			I.N_entries = entry::str2int(entry[1]);
		I.N_entries_str = entry[1];
	}
	else LOG.error("Expected Number entry as second field in FORMAT description: " + in);
	entry::tokenize(tokens[2], '=', entry);
	if (entry[0] == "Type")
	{
		if (entry[1] == "Integer") {I.Type = Integer;}
		else if ((entry[1] == "Float") || (entry[1] == "Numeric")) {I.Type = Float;}
		else if (entry[1] == "Character") {I.Type = Character;}
		else if (entry[1] == "String") {I.Type = String;}
		else if (entry[1] == "Flag")
		{
			I.Type = Flag;
			I.Type_str = "Flag";
			if (I.N_entries != 0) LOG.error("Flag Type must have 0 entries: " + in);
		}
		else LOG.error("Unknown Type in FORMAT meta-information: " + in);
	}
	else LOG.error("Expected Type entry as third field in FORMAT description: " + in);

	entry::tokenize(tokens[3], '=', entry);
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
		FORMAT_map[ FILTER_reverse_map[ I.ID ] ] = I;
		FORMAT_reverse_map[I.ID] = FILTER_reverse_map[ I.ID ];
		return 0;
	}
	else if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
	{
		FORMAT_map[ INFO_reverse_map[ I.ID ] ] = I;
		FORMAT_reverse_map[I.ID] = INFO_reverse_map[ I.ID ];
		return 0;
	}
	else if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
		return 0;
	else
	{
		FORMAT_map[ index ] = I;
		FORMAT_reverse_map[I.ID] = index;
		return 1;
	}
}

void header::add_CONTIG_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end-1);

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

	CONTIG_map[index] = I;
	CONTIG_reverse_map[I.ID] = index;
}

int header::add_FILTER_descriptor(const string &in, int index)
{
	size_t found_end=in.find_last_of(">");
	string details = in.substr(0, found_end-1);
	vector<string> tokens;
	entry::tokenize(details, ',', tokens);
	if (tokens.size() < 2)
		LOG.error("Expected 2 parts in FILTER definition: " + in);

	string Description;
	Field_description I;
	vector<string> entry;
	entry::tokenize(tokens[0], '=', entry);
	if (entry[0] == "ID") I.ID = entry[1];
	else LOG.error("Expected ID as first field in FILTER description: " + in);

	entry::tokenize(tokens[1], '=', entry);
	if (entry[0] == "Description")
	{
		Description = entry[1];
		for (unsigned int i=2; i<tokens.size(); i++)
		{
			Description += "; " + tokens[i];
		}
		I.Description = Description;
	}
	else
		LOG.error("Expected Description as second field in FILTER description: " + in);

	if ( INFO_reverse_map.find( I.ID ) != INFO_reverse_map.end() )
	{
		FILTER_map[ INFO_reverse_map[ I.ID ] ] = I;
		FILTER_reverse_map[I.ID] = INFO_reverse_map[ I.ID ];
		return 0;
	}
	else if ( FORMAT_reverse_map.find( I.ID ) != FORMAT_reverse_map.end() )
	{
		FILTER_map[ FORMAT_reverse_map[ I.ID ] ] = I;
		FILTER_reverse_map[I.ID] = FORMAT_reverse_map[ I.ID ];
		return 0;
	}
	else if ( FILTER_reverse_map.find( I.ID ) != FILTER_reverse_map.end() )
		return 0;
	else
	{
		FILTER_map[index] = I;
		FILTER_reverse_map[I.ID] = index;
		return 1;
	}
}


