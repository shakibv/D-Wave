/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains lattice-related functions (reading from file, etc).

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

#ifndef __LATTICE_H__
#define __LATTICE_H__

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>

#include "utils.h"

template <typename V, typename I>
class Lattice {
public:
	typedef V value_type;
	typedef I index_type;
private:
	struct Link {
		index_type s0;
		index_type s1;
		value_type cval;

		bool operator<(const Link& rhs) const
		{
			return std::fabs(cval) < std::fabs(rhs.cval);
		}
	};
public:
	Lattice(const std::string& lattice_file) : lattice_file(lattice_file)
	{
		std::ifstream fin;
		fin.open(lattice_file.c_str(), std::ios_base::in);
		if (!fin)
			throw std::runtime_error("cannot open file " + lattice_file + " to read lattice");

		maxs = 0;
		links.reserve(32768);

		bool first_line = true;
		while (1) {
			if (first_line) {
				std::string line;
				getline(fin, line, '\n');
				first_line = false;
				continue;
			}

			index_type s0, s1;
			value_type cval;
			fin >> s0 >> s1 >> cval;
			if (!fin) break;
			
			if (s0 < 0 || s1 < 0)
				throw std::runtime_error("negative spin index in file " + lattice_file);

			links.push_back({ index_type(s0), index_type(s1), cval });

			maxs = s0 > maxs ? s0 : maxs;
			maxs = s1 > maxs ? s1 : maxs;
		}

		fin.close();

		nsites = 0;
		std::vector<index_type> phys_sites(maxs + 1, index_type(-1));

		for (std::size_t i = 0; i < links.size(); ++i) {
			Link& link = links[i];

			if (phys_sites[link.s0] == index_type(-1))
				link.s0 = phys_sites[link.s0] = nsites++;
			else
				link.s0 = phys_sites[link.s0];

			if (phys_sites[link.s1] == index_type(-1))
				link.s1 = phys_sites[link.s1] = nsites++;
 			else
				link.s1 = phys_sites[link.s1];
		}

		// need this for higher ranges
		std::sort(links.begin(), links.end());
	}
	
	template <typename ST>
	void init_sites(std::vector<ST>& sites) const
	{
		sites.reserve(nsites);
		sites.resize(nsites);

		for (std::size_t i = 0; i < links.size(); ++i) {
			const Link& link = links[i];

			if (link.s0 == link.s1)
				sites[link.s0].hzv = link.cval;
			else {
				sites[link.s0].jzv.push_back(link.cval);
				sites[link.s0].neighbs.push_back(link.s1);
				++sites[link.s0].nneighbs;

				sites[link.s1].jzv.push_back(link.cval);
				sites[link.s1].neighbs.push_back(link.s0);
				++sites[link.s1].nneighbs;
			}
		}
	}

	template <typename ST>
	void init_sites(std::vector<ST>& sites, unsigned maxnb) const
	{
		sites.reserve(nsites);
		sites.resize(nsites);

		for (std::size_t i = 0; i < links.size(); ++i) {
			const Link& link = links[i];

			if (link.s0 == link.s1)
				sites[link.s0].hzv = link.cval;
			else {
				index_type& nneighb0 = sites[link.s0].nneighbs;
				if (nneighb0 >= maxnb)
					throw std::runtime_error("two many neighbors in lattice file "
						+ lattice_file + "; must be less or equal to " + to_s(maxnb));
				sites[link.s0].jzv[nneighb0] = link.cval;
				sites[link.s0].neighbs[nneighb0] = link.s1;
				++nneighb0;

				index_type& nneighb1 = sites[link.s1].nneighbs;
				if (nneighb1 >= maxnb)
					throw std::runtime_error("two many neighbors in lattice file "
						+ lattice_file + "; must be less or equal to " + to_s(maxnb));
				sites[link.s1].jzv[nneighb1] = link.cval;
				sites[link.s1].neighbs[nneighb1] = link.s0;
				++nneighb1;
			}
		}
	}
private:
	const std::string& lattice_file;

	std::size_t nsites;
	std::vector<Link> links;
	index_type maxs;
};

#endif
