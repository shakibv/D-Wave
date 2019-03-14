/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains a function that prints results.

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

#ifndef __OUTPUT_H__
#define __OUTPUT_H__

#include <map>
#include <string>
#include <vector>

template <typename value_type>
void print_results(const std::vector<value_type>& en,
	const std::string& latfile, unsigned rep0, unsigned nreps, bool lowest)
{
	std::map<value_type, std::size_t> map;
	for (std::size_t i = 0; i < en.size(); ++i) {
		/* typename std::map<value_type, std::size_t>::iterator it = map.find(en[i]); */

		auto it = map.begin();
		while (it != map.end()) {
			if (std::fabs(it->first - en[i]) < 1e-08)
				break;

			++it;
		}

		if (it == map.end())
			map[en[i]] = 1;
		else
			++it->second;
	}

	double scale = 1.0 / en.size();
	typename std::map<value_type, std::size_t>::const_iterator it = map.begin();
	for (; it != map.end(); ++it) {
		std::cout << std::setw(10) << it->first;
		std::cout << std::setw(10) << it->second;
		std::cout << std::setw(16) << double(it->second) * scale;
//		std::cout << std::setw(10) << rep0 << std::setw(10) << nreps;
		std::cout << "    " << latfile << "\n";

		if (lowest) break;
	}
}

#endif

