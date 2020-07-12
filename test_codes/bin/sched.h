/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains a function that reads from file or generates a schedule
to run simulated annealing.

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

#ifndef __SCHED_H__
#define __SCHED_H__

#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>

struct sched_entry {
	double beta;
};

inline std::vector<sched_entry> get_sched(const std::string sched_kind,
	unsigned nsweeps, double beta0, double beta1)
{
	std::vector<sched_entry> sched;

	if (sched_kind == "lin") {
		sched.resize(nsweeps);
		double bscale = nsweeps > 1 ? (beta1 - beta0) / (nsweeps - 1) : 0.0;
		for (std::size_t i = 0; i < nsweeps; ++i)
			sched[i].beta = beta0 + bscale * i;
	} else if (sched_kind == "exp") {
		sched.resize(nsweeps);
		sched[0].beta = beta0;
		double db = std::pow(beta1 / beta0, 1.0 / (nsweeps - 1));
		for (std::size_t i = 1; i < nsweeps; ++i)
			sched[i].beta = sched[i - 1].beta * db;
	} else {
		std::ifstream fin;
		fin.open(sched_kind.c_str(), std::ios_base::in);
		if (!fin)
			throw std::runtime_error("cannot open file " + sched_kind + " to read schedule");

		sched.reserve(10000);

		while (1) {
			double beta;
			fin >> beta;
			if (!fin) break;
			sched.push_back({beta});
		}

		fin.close();
	}

	return sched;
}

#endif

