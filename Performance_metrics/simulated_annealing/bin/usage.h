/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains a function that prints the usage of the programs.

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

#ifndef __USAGE_H__
#define __USAGE_H__

#include <string>
#include <iostream>
#include <stdexcept>

inline void usage(const std::string& msg, bool multi_threaded)
{
	std::cerr << "usage: " << "\n";
	std::cerr << "an.e -l lattice -s nsweeps -r nreps";
	std::cerr << " [-b0 beta0] [-b1 beta1] [-r0 rep0]";
	std::cerr << " [-v] [-sched sched_kind] [-t nthreads]\n";
	std::cerr << "where optional parameters are in square brackets\n";
	std::cerr << " -l lattice        --- lattice file\n";
	std::cerr << " -s nsweeps        --- number of sweeps\n";
	std::cerr << " -r nreps          --- number of repetitions\n";
	std::cerr << " -r0 rep0          --- start repetition; default value: 0\n";
	std::cerr << " -b0 beta0         --- initial inverse temperature; default value: 0.1\n";
	std::cerr << " -b1 beta1         --- final inverse temperature; default value: 3.0\n";
	std::cerr << " -sched sched_kind --- schedule kind: lin or exp or file name; default value: lin\n";
	std::cerr << " -v                --- verbose mode; prints some info including timing info\n";
    std::cerr << " -g                --- prints only the lowest energy solution\n";
	if (multi_threaded)
		std::cerr << " -t nthreads       --- number of threads\n";

	if (!msg.empty())
		throw std::runtime_error(msg);
}

#endif

