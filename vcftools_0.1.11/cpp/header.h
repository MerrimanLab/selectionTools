/*
 * header.h
 *
 *  Created on: Apr 29, 2013
 *      Author: amarcketta
 */

#ifndef HEADER_H_
#define HEADER_H_

#include <cstring>
#include <string>
#include <map>
#include "entry.h"

using namespace std;

class header
{
public:
	map<int, Field_description> INFO_map;
	map<int, Field_description> FILTER_map;
	map<int, Field_description> FORMAT_map;
	map<int, Field_description> CONTIG_map;
	map<string, int> CONTIG_reverse_map;
	map<string, int> FILTER_reverse_map;
	map<string, int> INFO_reverse_map;
	map<string, int> FORMAT_reverse_map;

	header() {};
	~header() {};

	int add_INFO_descriptor(const string &in, int index);
	int add_FILTER_descriptor(const string &in, int index);
	int add_FORMAT_descriptor(const string &in, int index);
	void add_CONTIG_descriptor(const string &in, int index);
};
//
//class Field_description
//{
//public:
//	string ID;
//	int N_entries;
//	string N_entries_str;
//	string Type_str;
//	Type_enum Type;
//	string Description;
//	string Length;
//	string Assembly;
//
//	Field_description() : ID(""), N_entries(0), Type(Integer), Description("") {};
//	~Field_description() {};
//};

#endif /* HEADER_H_ */
