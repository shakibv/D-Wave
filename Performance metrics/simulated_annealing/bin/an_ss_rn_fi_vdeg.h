/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of single-spin simulated annealing algorithm for
Ising spin glasses with range-n interactions with magnetic field
and any number of neighbors.

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

template <typename T = uint64_t, std::size_t depth = 18>
class Algorithm {
public:
	typedef T word_type;
	typedef int value_type;
	typedef unsigned index_type;

	typedef Lattice<value_type, index_type> lattice_type;

	struct site_type {
		value_type spin;
		value_type hzv;
		std::vector<value_type> jzv;
		value_type de;
		index_type nneighbs;
		std::vector<index_type> neighbs;
	};

	struct sched_type {
		std::vector<unsigned> r;
	};

	static const std::size_t word_size = 1;
	static const std::size_t offs = 8 * sizeof(word_type) - depth;
	static const std::size_t rand_size = std::size_t(1) << depth;

	typedef bitgen_lincon<word_type> bgen_type;
//	typedef bitgen_lagfib<word_type> bgen_type;

	Algorithm() {}

	template <typename SE>
	Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
	{
		lattice.init_sites(sites);

		unsigned maxh = 0;
		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];

			unsigned mh = 0;
			for (std::size_t k = 0; k < site.nneighbs; ++k)
				mh += std::abs(site.jzv[k]);
			mh += std::abs(site.hzv);
			if (mh > maxh) maxh = mh;
		}

		sched.resize(sched0.size());
		for (std::size_t sweep = 0; sweep < sched0.size(); ++sweep) {
			double p0 = std::exp(-2 * sched0[sweep].beta);
			double p = 1.0;

			sched[sweep].r.resize(maxh + 1);
			for (unsigned k = 1; k <= maxh; ++k) {
				p *= p0;
				sched[sweep].r[k] = rand_size * p;
			}
		}
	}

	void reset_sites(std::size_t rep)
	{
		rgen.seed(rep + 1);
		bgen.seed(rep + 1);

		// startconf
		for (std::size_t i = 0; i < sites.size(); ++i)
			sites[i].spin = 2 * int((rgen() >> 29) & 1) - 1;

		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];
			value_type h = site.hzv;
			for (std::size_t k = 0; k < site.nneighbs; ++k)
				h += site.jzv[k] * sites[site.neighbs[k]].spin;
			site.de = -h * site.spin;
		}
	}

	void do_sweep(size_t sweep)
	{
		for (auto& site : sites)
			update_site(site, sched[sweep]);
	}

	std::size_t get_energies(std::vector<value_type>& en, std::size_t offs) const
	{
		en[offs] = calc_energy();
		return offs + 1;
	}

	std::string get_info() const
	{
		return "algorithm: single-spin, range-n couplings, with fields";
	}
private:
	std::vector<site_type> sites;
	std::vector<sched_type> sched;

	std::mt19937 rgen;
	bgen_type bgen;

	void update_site(site_type& site, const sched_type& sched)
	{
		if (site.de <= 0 || sched.r[site.de] > (bgen() >> offs)) {
			site.spin = -site.spin;
			site.de = -site.de;

			for (std::size_t k = 0; k < site.nneighbs; ++k) {
				site_type& neighbor = sites[site.neighbs[k]];
				neighbor.de -= 2 * site.jzv[k] * site.spin * neighbor.spin;
			}
		}
	}

	value_type calc_energy() const
	{
		value_type en = 0;

		for (std::size_t i = 0; i < sites.size(); ++i) {
			const site_type& site = sites[i];
			value_type h = site.hzv;
			for (std::size_t k = 0; k < site.nneighbs; ++k) {
				std::size_t l = site.neighbs[k];
				if (i < l) h += sites[l].spin * site.jzv[k];
			}

			en += h * site.spin;
		}

		return en;
	}
};

#endif
