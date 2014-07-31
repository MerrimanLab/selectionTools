/*
 * log.h
 *
 *  Created on: Nov 11, 2009
 *      Author: Adam Auton
 *      ($Revision: 91 $)
 */

#ifndef LOG_H_
#define LOG_H_

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <sstream>

using namespace std;

class output_log
{
public:
	output_log();
	~output_log() {};
	void open(const string &filename);
	void close();
	void printLOG(string s);
	void error(string err_msg, int error_code=0);
	void error(string err_msg, double value1, double value2, int error_code=0);
	void warning(string err_msg);
	void one_off_warning(string err_msg);
	void set_screen_output(bool do_screen_output);
	static string int2str(int n);
	static string longint2str(long int n);
	static string dbl2str(double n, int prc);
	static string dbl2str_fixed(double n, int prc);
private:
	bool output_to_screen;
	ofstream LOG;
};

#endif /* LOG_H_ */
