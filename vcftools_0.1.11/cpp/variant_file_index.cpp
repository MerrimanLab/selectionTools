/*
 * variant_file_index.cpp
 *
 *  Created on: 3 Aug 2011
 *      Author: auton
 */

#include "variant_file.h"

bool variant_file::read_index_file(const string &index_filename)
{
	// Check index is newer than vcf file
	struct stat stat_idx, stat_vcf;
	stat(index_filename.c_str(), &stat_idx);
	stat(filename.c_str(), &stat_vcf);
	if (stat_vcf.st_mtime > stat_idx.st_mtime)
	{
		LOG.warning("Index file is older than variant file. Will regenerate.");
		return false;
	}

	LOG.printLOG("Reading Index file.\n");
	big_endian_machine = is_big_endian();
	gzFile in = gzopen(index_filename.c_str(), "rb");
	if (in == NULL)
		return false;

	char magic[7];
	idx_read(in, magic, 7, sizeof(char));
	if (strncmp(magic, "VCFIDX\1", 7) != 0)
	{	// Doesn't appear to be an index file
		gzclose(in);
		LOG.warning("Index file doesn't appear to be valid. Will (try to) overwrite.\n");
		return false;
	}

	uint32_t tmp;
	uint64_t tmp64;
	idx_read(in, &tmp, 1, sizeof(uint32_t));
	N_entries = tmp;
	idx_read(in, &tmp, 1, sizeof(uint32_t));
	N_indv = tmp;
	idx_read(in, &tmp, 1, sizeof(uint32_t));
	unsigned int l_meta = tmp;
	idx_read(in, &tmp, 1, sizeof(uint32_t));
	unsigned int l_indv = tmp;

	char *meta_buffer = new char [l_meta+1];
	char *indv_buffer = new char [l_indv+1];

	idx_read(in, meta_buffer, l_meta, 1);
	idx_read(in, indv_buffer, l_indv, 1);

	// Split the strings
	meta.resize(0);
	char * pch;
	pch = strtok(meta_buffer,"\n");
	while (pch != NULL)
	{
		meta.push_back(pch);
		pch = strtok(NULL, "\n");
	}

	indv.resize(0);
	pch = strtok (indv_buffer,"\n");
	while (pch != NULL)
	{
		indv.push_back(pch);
		pch = strtok (NULL, "\n");
	}

	delete [] indv_buffer;
	delete [] meta_buffer;

	entry_file_locations.resize(N_entries);
	for (unsigned int ui=0; ui<N_entries; ui++)
	{
		idx_read(in, &tmp64, 1, sizeof(uint64_t));
		entry_file_locations[ui] = tmp64;
	}

	gzclose(in);
	return true;
}

void variant_file::write_index_file(const string &index_filename)
{
	LOG.printLOG("Writing Index file.\n");
	big_endian_machine = is_big_endian();

	gzFile out = gzopen(index_filename.c_str(), "wb");
	if (out == NULL)
	{
		LOG.warning("Could not write index file.\n");
		return;
	}

	unsigned int l_meta = 0;
	for (unsigned int ui=0; ui<meta.size(); ui++)
		l_meta += meta[ui].size() + 1;

	unsigned int l_indv = 0;
	for (unsigned int ui=0; ui<indv.size(); ui++)
		l_indv += indv[ui].size() + 1;

	idx_write(out, (char*)"VCFIDX\1", 7, 1);

	uint32_t tmp;
	uint64_t tmp64;
	tmp = N_entries; idx_write(out, &tmp, 1, sizeof(uint32_t));
	tmp = N_indv; idx_write(out, &tmp, 1, sizeof(uint32_t));
	tmp = l_meta; idx_write(out, &tmp, 1, sizeof(uint32_t));
	tmp = l_indv; idx_write(out, &tmp, 1, sizeof(uint32_t));

	for (unsigned int ui=0; ui<meta.size(); ui++)
	{
		string str = meta[ui] + "\n";
		idx_write(out, (char*)(str.c_str()), str.size(), 1);
	}

	for (unsigned int ui=0; ui<indv.size(); ui++)
	{
		string str = indv[ui] + "\n";
		idx_write(out, (char*)(str.c_str()), str.size(), 1);
	}

	for (unsigned int ui=0; ui<entry_file_locations.size(); ui++)
	{
		tmp64 = entry_file_locations[ui];
		idx_write(out, &tmp64, 1, sizeof(uint64_t));
	}

	gzclose(out);
}

int variant_file::idx_read(gzFile &in, void *buffer, unsigned int len, size_t size)
{
	int ret = gzread(in, buffer, size*len);
	if ((big_endian_machine) && (size > 1))	// Note: don't bother swapping character arrays - index is defined as little endian.
	{
		unsigned int ui;
		for (ui=0; ui<len; ui++)
			ByteSwap((unsigned char *)buffer+(size*ui), size);
	}
	return ret;
}

void variant_file::idx_write(gzFile &out, void *buffer, unsigned int len, size_t size)
{
	if ((big_endian_machine) && (size > 1))
	{
		unsigned int ui;
		for (ui=0; ui<len; ui++)
			ByteSwap((unsigned char *)buffer+(size*ui), size);
	}
	gzwrite(out, buffer, size*len);
}
