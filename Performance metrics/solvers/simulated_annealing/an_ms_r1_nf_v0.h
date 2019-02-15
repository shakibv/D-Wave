/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of multi-spin simulated annealing algorithm for
Ising spin glasses with range-1 interactions with magnetic field
using approach two, see the paper for more details.

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
#include "utils.h"

#define OMP_VERSION_2

template <typename T = uint64_t>
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

	static const std::size_t word_size = 8 * sizeof(word_type);
	static const std::size_t lastbit = word_size - 1;

	typedef bitgen_lincon<word_type> bgen_type;
//	typedef bitgen_lagfib<word_type> bgen_type;

	Algorithm() {}

	template <typename SE>
	Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
	{
		lattice.init_sites(sites, MAXNB);

		maxnb = 0;

		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];

			if (site.nneighbs > maxnb) maxnb = site.nneighbs;

			std::size_t l = 0;
			for (; l < site.nneighbs; ++l)
				site.jzw[l] = site.jzv[l] == -1 ? word_type(-1) : 0;

			value_type cval = -1;
			for (; l < MAXNB; ++l) {
				site.jzv[l] = 0;
				site.jzw[l] = cval == -1 ? word_type(-1) : 0;
				site.neighbs[l] = i;

				cval = -cval;
			}
		}

		word_type p;
		word_type pch;
		sched.resize(sched0.size());
		for (std::size_t sweep 	= 0; sweep < sched0.size(); ++sweep) {
			sched[sweep].resize(2);
			double beta = sched0[sweep].beta;

			double p1 = std::exp(-2 * beta);
			double p0 = p1;

			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][1].r2[i] = p2mask(pch, i);

			p1 *= p0;
			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][0].r2[i] = p2mask(pch, i);

			p1 *= p0;
			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][1].r1[i] = p2mask(pch, i);

			p1 *= p0;
			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][0].r1[i] = p2mask(pch, i);

			p1 *= p0;
			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][1].r0[i] = p2mask(pch, i);

			p1 *= p0;
			p = beta != 0.0 ? word_type(word_type(-1) * p1) : word_type(-1);
			pch = p ^ (p >> 1);
			for (std::size_t i = 0; i < word_size; ++i)
				sched[sweep][0].r0[i] = p2mask(pch, i);
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

	void do_sweep(size_t sweep)
	{
		if (maxnb <= 4)
			for (std::size_t i = 0; i < sites.size(); ++i) {
				site_type& site = sites[i];
				update_site4(site, sched[sweep]);
			}
		else
			for (std::size_t i = 0; i < sites.size(); ++i) {
				site_type& site = sites[i];
				update_site6(site, sched[sweep]);
			}
	}

	std::size_t get_energies(std::vector<value_type>& en, std::size_t offs) const
	{
		calc_energies(en, offs);
		return offs + word_size;
	}

	std::string get_info() const
	{
		return "algorithm: multi-spin, range-1 couplings, without fields, vesion 0";
	}
private:
	std::vector<site_type> sites;

	std::mt19937 rgen;
	bgen_type bgen;

	unsigned maxnb;

	struct sched_type {
		word_type r0[word_size];
		word_type r1[word_size];
		word_type r2[word_size];
	};

	std::vector<std::vector<sched_type> > sched;

	word_type p2mask(word_type p, std::size_t i) const
	{
		return -((p >> i) & 1);
	}

	word_type flippable4(word_type mask, word_type mask1, word_type mask2, const sched_type& sched)
	{
		word_type s = mask;
		word_type r = (sched.r1[lastbit] & mask1) | (sched.r2[lastbit] & mask2);

		for (std::size_t i = lastbit - 1; s; --i) {
			s &= bgen();
			r ^= s & ((sched.r1[i] & mask1) | (sched.r2[i] & mask2));
		}

		return r;
	}

	word_type flippable6(word_type mask, word_type mask0, word_type mask1, word_type mask2, const sched_type& sched)
	{
		word_type s = mask;
		word_type r = (sched.r0[lastbit] & mask0) | (sched.r1[lastbit] & mask1) | (sched.r2[lastbit] & mask2);

		for (std::size_t i = lastbit - 1; s; --i) {
			s &= bgen();
			r ^= s & ((sched.r0[i] & mask0) | (sched.r1[i] & mask1) | (sched.r2[i] & mask2));
		}

		return r;
	}

	void update_site4(site_type& site, const std::vector<sched_type>& sched)
	{
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		word_type l3 = site.jzw[3] ^ (spin ^ sites[site.neighbs[3]].spin);

		word_type t0 = l0 & l1;
		l1 = l0 | l1;
		word_type t1 = l1 & l2;
		l2 = l1 | l2;
		word_type t2 = l2 & l3;
		l0 = t0 & t1;
		t1 = t0 | t1;
		l1 = t1 & t2;
		t0 = l0 & l1;
		l1 = l0 | l1;

		word_type mask1 = t0;
		word_type mask2 = l1 & ~mask1;

		word_type mask = mask1 | mask2;

		index_type k = site.nneighbs % 2;
		site.spin = spin ^ (~mask | flippable4(mask, mask1, mask2, sched[k]));
	}

	void update_site6(site_type& site, const std::vector<sched_type>& sched)
	{
		word_type spin = site.spin;

		word_type l0 = site.jzw[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		word_type l1 = site.jzw[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		word_type l2 = site.jzw[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		word_type l3 = site.jzw[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		word_type l4 = site.jzw[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		word_type l5 = site.jzw[5] ^ (spin ^ sites[site.neighbs[5]].spin);

		word_type t0 = l0 & l1;
		l1 = l0 | l1;
		word_type t1 = l1 & l2;
		l2 = l1 | l2;
		word_type t2 = l2 & l3;
		l3 = l2 | l3;
		word_type t3 = l3 & l4;
		l4 = l3 | l4;
		word_type t4 = l4 & l5;
		l0 = t0 & t1;
		t1 = t0 | t1;
		l1 = t1 & t2;
		t2 = t1 | t2;
		l2 = t2 & t3;
		t3 = t2 | t3;
		l3 = t3 & t4;
		t0 = l0 & l1;
		l1 = l0 | l1;
		t1 = l1 & l2;
		l2 = l1 | l2;
		t2 = l2 & l3;
		l0 = t0 & t1;
		t1 = t0 | t1;
		l1 = t1 & t2;
		t2 = t1 | t2;
		t0 = l0 & l1;
		l1 = l0 | l1;

		word_type mask0 = t0;
		word_type mask1 = l1 & ~mask0;
		word_type mask2 = t2 & ~mask0 & ~mask1;

		word_type mask = mask0 | mask1 | mask2;

		index_type k = site.nneighbs % 2;
		site.spin = spin ^ (~mask | flippable6(mask, mask0, mask1, mask2, sched[k]));
	}

	void calc_energies(std::vector<value_type>& en, std::size_t offs) const
	{
		for (unsigned k = 0; k < word_size; ++k)
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
};

#endif
