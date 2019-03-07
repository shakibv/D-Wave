/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains miscellaneous utility codes: conversion of numbers to
strings and vice versa, optional class, parsing of command-line
arguments, time function to measure time intervals.

---------------------------------------------------------------------

Copyright (C) 2012-2013 by Sergei Isakov <isakov@itp.phys.ethz.ch>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/

#ifndef __UTILS_H__
#define __UTILS_H__

#include <map>
#include <string>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <sys/time.h>
#include <time.h>

/*** strings ******************************************************************/

template <typename T>
inline std::string to_s(const T& val)
{
	std::ostringstream ss;
	ss << val;
	return ss.str();
}

template <typename C>
inline C to_val(const char* str);

template <>
inline int to_val(const char* str)
{
	return std::atol(str);
}

template <>
inline unsigned to_val(const char* str)
{
	return std::atol(str);
}

template <>
inline double to_val(const char* str)
{
	return std::atof(str);
}

/*** optional *****************************************************************/

template <typename T>
class opt {
private:
	struct bool_conv { int dummy; };
public:
	opt() : initialized(false) { }
	opt(const T& data) : initialized(true), data(data) {}

	T& operator*()
	{
		return data;
	}

	const T& operator*() const
	{
		return data;
	}

	T* operator->()
	{
		return &data;
	}

	const T* operator->() const
	{
		return &data;
	}

	// ugly conversion to bool
	operator int bool_conv::* () const
	{
		return initialized ? &bool_conv::dummy : 0;
	}
private:
	bool initialized;
	T data;
};

/*** command line arguments ***************************************************/

typedef std::map<std::string, std::string> amap_type;

inline amap_type parse_args(int argc, char *argv[])
{
	// note: negative numbers cannot be parsed

	amap_type args;

	bool have_key = false;
	std::string key;
	for (std::size_t i = 1; i < std::size_t(argc); ++i) {
		if (argv[i][0] == '-') {
			if (have_key)
				args[key] = "1";
			else
				have_key = true;
			key = argv[i] + 1;
		} else if (have_key) {
			args[key] = argv[i];
			have_key = false;
		}
	}

	if (have_key) args[key] = "1";

	return args;
}

inline opt<std::string> get_sarg(const amap_type& args, const std::string& o)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<std::string>();

	return opt<std::string>(it->second);
}

inline opt<std::string> get_sarg(const amap_type& args, const std::string& o, const std::string& def)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<std::string>(def);

	return opt<std::string>(it->second);
}

inline opt<int> get_iarg(const amap_type& args, const std::string& o)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<int>();

	return opt<int>(std::atoi(it->second.c_str()));
}

inline opt<int> get_iarg(const amap_type& args, const std::string& o, int def)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<int>(def);

	return opt<int>(std::atoi(it->second.c_str()));
}

inline opt<unsigned> get_uarg(const amap_type& args, const std::string& o)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<unsigned>();

	return opt<unsigned>(std::atoi(it->second.c_str()));
}

inline opt<unsigned> get_uarg(const amap_type& args, const std::string& o, unsigned def)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<unsigned>(def);

	return opt<unsigned>(std::atoi(it->second.c_str()));
}

inline opt<double> get_darg(const amap_type& args, const std::string& o)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<double>();

	return opt<double>(std::atof(it->second.c_str()));
}

inline opt<double> get_darg(const amap_type& args, const std::string& o, double def)
{
	amap_type::const_iterator it = args.find(o);
	if (it == args.end()) return opt<double>(def);

	return opt<double>(std::atof(it->second.c_str()));
}

/*** time *********************************************************************/

inline double get_time()
{
	struct timeval tv;
	if (gettimeofday(&tv, NULL) != 0)
		return 0.0;
	else
		return tv.tv_sec + 1e-6 * tv.tv_usec;
}

#endif

