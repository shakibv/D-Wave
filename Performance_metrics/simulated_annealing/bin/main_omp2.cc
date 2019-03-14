/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Main function for multi-threaded codes using OPENMP.

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

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <iomanip>
#include <stdexcept>

#ifdef _OPENMP
#	include "omp.h"
#else
#	error "openmp is required"
#endif

#include "sched.h"
#include "usage.h"
#include "utils.h"
#include "output.h"

#ifndef ALGORITHM
#error "Please specify the algorithm"
#else
#include ALGORITHM
#endif

int main(int argc, char *argv[])
{
	try {
		double t0 = get_time();

		// command line arguments

		amap_type args = parse_args(argc, argv);

		opt<std::string> latfile = get_sarg(args, "l");
		if (!latfile) usage("lattice is not provided", false);
		opt<unsigned> nsweeps = get_uarg(args, "s");
		opt<unsigned> nreps = get_uarg(args, "r");
		if (!nreps) usage("nreps is not provided", false);
		opt<double> beta0 = get_darg(args, "b0", 0.1);
		opt<double> beta1 = get_darg(args, "b1", 3.0);
		opt<unsigned> rep0 = get_uarg(args, "r0", 0);
		opt<unsigned> verbose = get_uarg(args, "v", 0);
		opt<unsigned> lowest = get_uarg(args, "g", 0);
		opt<unsigned> nthreads = get_uarg(args, "t", omp_get_max_threads());
		opt<std::string> sched_kind = get_sarg(args, "sched", "lin");
		bool def_sched = *sched_kind == "lin" || *sched_kind == "exp";
		if (!nsweeps && def_sched)
			usage("nsweeps is not provided", false);

		typedef Algorithm<> alg_type;
		typedef alg_type::lattice_type lattice_type;

		// read lattice

		lattice_type lattice(*latfile);

		// schedule

		std::vector<sched_entry> sched = get_sched(*sched_kind, *nsweeps, *beta0, *beta1);
		*nsweeps = sched.size();

		// init annealing

		unsigned n = std::min(*nthreads, *nreps - *rep0);
		std::vector<alg_type> algs(n);

#ifdef OMP_VERSION_2
		#pragma omp parallel num_threads(n)
		{
			unsigned m = omp_get_thread_num();
			algs[m] = alg_type(lattice, sched);
		}
#else
		algs[0] = alg_type(lattice, sched);
#endif

		typedef alg_type::value_type value_type;
		std::vector<value_type> en(*nreps * alg_type::word_size, 0);

		if (*verbose) {
			if (def_sched)
				std::cout << "#" << *sched_kind << " schedule: nsweeps="
					<< *nsweeps << " b0=" << *beta0 << " b1=" << *beta1;
			else
				std::cout << "#schedule from file " << *sched_kind
					<< ": nsweeps=" << *nsweeps;
			std::cout << "; rep0=" << *rep0 << " nreps=" << *nreps << "\n";
			std::cout << "#" << algs[0].get_info() << "\n";
			std::cout << "#running " << algs.size() << " omp threads" << "\n";
		}

		double t1 = get_time();
		if (*verbose) std::cout << "#init done in " << t1 - t0 << " s\n";

		double t2 = get_time();

		#pragma omp parallel num_threads(n)
		{
			unsigned n = omp_get_num_threads();
			unsigned m = omp_get_thread_num();

#ifndef OMP_VERSION_2
			if (m > 0) algs[m] = algs[0];
#endif

			std::size_t r0 = *rep0 + *nreps * m / n;
			std::size_t r1 = *rep0 + *nreps * (m + 1) / n;
			std::size_t offs = *nreps * m / n * alg_type::word_size;

			// main loop

			for (std::size_t rep = r0; rep < r1; ++rep) {
				algs[m].reset_sites(rep);
				for (std::size_t sweep = 0; sweep < *nsweeps; ++sweep)
					algs[m].do_sweep(sweep);

				offs = algs[m].get_energies(en, offs);
			}
		}

		double t3 = get_time();
		if (*verbose) std::cout << "#work done in " << t3 - t2 << " s\n";

		double t4 = get_time();

		// print results

		print_results(en, *latfile, *rep0, *nreps, *lowest);

		double t5 = get_time();
		if (*verbose) std::cout << "#outp done in " << t5 - t4 << " s\n";
	} catch (std::exception& e) {
		std::cerr << "error: " << e.what() << std::endl;
	} catch (...) {
		std::cerr << "unknown error" << std::endl;
	}

	return 0;
}