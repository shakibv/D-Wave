/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of multi-spin simulated annealing algorithm for
Ising spin glasses with range-1 interactions without magnetic field
using approach one, see the paper for more details.

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

#ifndef __ALGORITHM_H__
#define __ALGORITHM_H__

#include <cmath>
#include <random>
#include <vector>
#include <string>

#include "bits.h"
#include "lattice.h"
#include "ms_config.h"
#include "utils.h"

template <typename T = uint64_t, std::size_t depth = 18>
class Algorithm {
public:
	typedef T word_type;
	typedef int value_type;
	typedef unsigned index_type;

	static const unsigned MAXNB = 6;

	typedef Lattice<value_type, index_type> lattice_type;

	struct site_type {
		word_type spin;
		word_type hzw;
		word_type jzw[MAXNB];
		value_type hzv;
		value_type jzv[MAXNB];
		index_type nneighbs;
		index_type neighbs[MAXNB];
	};

	struct sched_type {
		unsigned r6;
		unsigned r4;
		unsigned r2;
		unsigned r5;
		unsigned r3;
		unsigned r1;
	};

	static const std::size_t word_size = 8 * sizeof(word_type);
	static const std::size_t rand_size = std::size_t(1) << depth;

	typedef bitgen_lincon<word_type> bgen_type;
//	typedef bitgen_lagfib<word_type> bgen_type;

	Algorithm() {}

	template <typename SE>
	Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
	{
		lattice.init_sites(sites, MAXNB);

		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];
			
			if (!check_number_of_neighbors(site.nneighbs))
				throw std::runtime_error(to_s(site.nneighbs) +
					" neighbors is not defined in ms_config.h");

			std::size_t l = 0;
			for (; l < site.nneighbs; ++l)
				site.jzw[l] = site.jzv[l] == 1 ? word_type(-1) : 0;
		}

		sched.resize(sched0.size());
		for (std::size_t sweep = 0; sweep < sched0.size(); ++sweep) {
			double p = std::exp(-2 * sched0[sweep].beta);
			double p0 = p;
			sched[sweep].r1 = rand_size * p;
			p *= p0;
			sched[sweep].r2 = rand_size * p;
			p *= p0;
			sched[sweep].r3 = rand_size * p;
			p *= p0;
			sched[sweep].r4 = rand_size * p;
			p *= p0;
			sched[sweep].r5 = rand_size * p;
			p *= p0;
			sched[sweep].r6 = rand_size * p;
		}
	}

	void reset_sites(std::size_t rep)
	{
		rgen.seed(rep + 1);
		bgen.seed(rep + 1);

		// startconf
		for (std::size_t i = 0; i < sites.size(); ++i)
			sites[i].spin = random_word(rgen, sites[i].spin);
	}

	void do_sweep(std::size_t sweep)
	{
		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];
			switch (site.nneighbs) {
#	ifdef USE_1_NEIGHB
			case 1:
				update_site1(site, sched[sweep]);
				break;
#	endif
#	ifdef USE_2_NEIGHB
			case 2:
				update_site2(site, sched[sweep]);
				break;
#	endif
#	ifdef USE_3_NEIGHB
			case 3:
				update_site3(site, sched[sweep]);
				break;
#	endif
#	ifdef USE_4_NEIGHB
			case 4:
				update_site4(site, sched[sweep]);
				break;
#	endif
#	ifdef USE_5_NEIGHB
			case 5:
				update_site5(site, sched[sweep]);
				break;
#	endif
#	ifdef USE_6_NEIGHB
			case 6:
				update_site6(site, sched[sweep]);
				break;
#	endif
			}
		}
	}

	std::size_t get_energies(std::vector<value_type>& en, std::size_t offs) const
	{
		calc_energies(en, offs);
		return offs + word_size;
	}

	std::string get_info() const
	{
		return "algorithm: multi-spin, range-1 couplings, without fields";
	}
private:
	std::vector<site_type> sites;
	std::vector<sched_type> sched;

	std::mt19937 rgen;
	bgen_type bgen;

	void update_site1(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);

		if (r >= sched.r1) {
			word_type mask = l0;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void update_site2(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);

		word_type j0 = l0 ^ l1;
		word_type j1 = l0 & l1;

		if (r >= sched.r2) {
			word_type mask = j0 | j1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void update_site3(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);

		word_type j1 = l0 ^ l1;
		word_type j0 = j1 ^ l2;
		j1 = (l0 & l1) ^ (j1 & l2);

		if (r >= sched.r1) {
			word_type mask = j1;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = j1 | j0;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void update_site4(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		word_type l3 = site.jzw[3] ^ (spin ^ sites[site.neighbs[3]].spin);

		word_type j0 = l0 ^ l1;
		word_type j1 = l0 & l1;
		word_type j2 = l2 ^ l3;
		word_type j3 = l2 & l3;


		if (r >= sched.r2) {
			word_type mask = j1 | j3 | (j0 & j2);
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = j1 | j3 | j0 | j2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void update_site5(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		word_type l3 = site.jzw[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		word_type l4 = site.jzw[4] ^ (spin ^ sites[site.neighbs[4]].spin);

		word_type j1 = l0 ^ l1;
		word_type j0 = j1 ^ l2;
		j1 = (l0 & l1) ^ (j1 & l2);

		word_type j2 = l3 ^ l4;
		word_type j3 = l3 & l4;

		if (r >= sched.r1) {
			word_type mask = ((j1 | j3) & (j0 | j2)) | (j1 & j3);
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (j0 & j2) | j1 | j3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = j0 | j2 | j1 | j3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void update_site6(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		word_type l3 = site.jzw[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		word_type l4 = site.jzw[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		word_type l5 = site.jzw[5] ^ (spin ^ sites[site.neighbs[5]].spin);

		word_type j1 = l0 ^ l1;
		word_type j0 = j1 ^ l2;
		j1 = (l0 & l1) ^ (j1 & l2);

		word_type j3 = l3 ^ l4;
		word_type j2 = j3 ^ l5;
		j3 = (l3 & l4) ^ (j3 & l5);

		if (r >= sched.r2) {
			word_type mask = ((j1 | j3) & (j0 | j2)) | (j1 & j3);
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (j0 & j2) | j1 | j3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = j0 | j2 | j1 | j3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}

	void calc_energies(std::vector<value_type>& en, std::size_t offs) const
	{
		for (unsigned k = 0; k < word_size; ++k) {
			for (std::size_t i = 0; i < sites.size(); ++i) {
				const site_type& site = sites[i];

				int spin = 2 * int((site.spin >> k) & 1) - 1;

				value_type h = 0;
				for (std::size_t l = 0; l < site.nneighbs; ++l) {
					std::size_t j = site.neighbs[l];
					if (i > j) continue;

					int nspin = 2 * int((sites[j].spin >> k) & 1) - 1;
					h += site.jzv[l] * nspin;
				}

				en[offs + k] += h * spin;
			}
		}
	}
};

#endif
