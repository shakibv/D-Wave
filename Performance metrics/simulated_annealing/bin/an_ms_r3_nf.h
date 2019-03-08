/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of multi-spin simulated annealing algorithm for
Ising spin glasses with range-3 interactions without magnetic field
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
		word_type jzw0[MAXNB];
		word_type jzw1[MAXNB];
		value_type hzv;
		value_type jzv[MAXNB];
		index_type nneighbs;
		index_type neighbs[MAXNB];
		unsigned cs;
	};

	struct sched_type {
		unsigned r1;
		unsigned r2;
		unsigned r3;
		unsigned r4;
		unsigned r5;
		unsigned r6;
		unsigned r7;
		unsigned r8;
		unsigned r9;
		unsigned r10;
		unsigned r11;
		unsigned r12;
		unsigned r13;
		unsigned r14;
		unsigned r15;
		unsigned r16;
		unsigned r17;
		unsigned r18;
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

			site.cs = 0;
			for (std::size_t l = 0; l < site.nneighbs; ++l) {
				switch (site.jzv[l]) {
				case 1:
					site.jzw0[l] = word_type(-1);
					site.jzw1[l] = 0;
					break;
				case 2:
					site.jzw0[l] = 0;
					site.jzw1[l] = word_type(-1);
					break;
				case 3:
					site.jzw0[l] = word_type(-1);
					site.jzw1[l] = word_type(-1);
					break;
				}

				site.cs = std::abs(site.jzv[l]) + 4 * site.cs;
			}
			site.cs += 1000000 * site.nneighbs;
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
			p *= p0;
			sched[sweep].r7 = rand_size * p;
			p *= p0;
			sched[sweep].r8 = rand_size * p;
			p *= p0;
			sched[sweep].r9 = rand_size * p;
			p *= p0;
			sched[sweep].r10 = rand_size * p;
			p *= p0;
			sched[sweep].r11 = rand_size * p;
			p *= p0;
			sched[sweep].r12 = rand_size * p;
			p *= p0;
			sched[sweep].r13 = rand_size * p;
			p *= p0;
			sched[sweep].r14 = rand_size * p;
			p *= p0;
			sched[sweep].r15 = rand_size * p;
			p *= p0;
			sched[sweep].r16 = rand_size * p;
			p *= p0;
			sched[sweep].r17 = rand_size * p;
			p *= p0;
			sched[sweep].r18 = rand_size * p;
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
		// machine-generated code; do not edit
		for (std::size_t i = 0; i < sites.size(); ++i) {
			site_type& site = sites[i];
			switch (site.cs) {
		#ifdef USE_1_NEIGHB
			case 1000001:
				update_site1_1(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_1_NEIGHB
			case 1000002:
				update_site1_2(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_1_NEIGHB
			case 1000003:
				update_site1_3(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000005:
				update_site2_11(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000006:
				update_site2_12(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000010:
				update_site2_22(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000007:
				update_site2_13(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000011:
				update_site2_23(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_2_NEIGHB
			case 2000015:
				update_site2_33(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000021:
				update_site3_111(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000022:
				update_site3_112(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000026:
				update_site3_122(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000042:
				update_site3_222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000023:
				update_site3_113(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000027:
				update_site3_123(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000043:
				update_site3_223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000031:
				update_site3_133(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000047:
				update_site3_233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_3_NEIGHB
			case 3000063:
				update_site3_333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000085:
				update_site4_1111(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000086:
				update_site4_1112(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000090:
				update_site4_1122(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000106:
				update_site4_1222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000170:
				update_site4_2222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000087:
				update_site4_1113(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000091:
				update_site4_1123(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000107:
				update_site4_1223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000171:
				update_site4_2223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000095:
				update_site4_1133(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000111:
				update_site4_1233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000175:
				update_site4_2233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000127:
				update_site4_1333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000191:
				update_site4_2333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_4_NEIGHB
			case 4000255:
				update_site4_3333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000341:
				update_site5_11111(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000342:
				update_site5_11112(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000346:
				update_site5_11122(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000362:
				update_site5_11222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000426:
				update_site5_12222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000682:
				update_site5_22222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000343:
				update_site5_11113(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000347:
				update_site5_11123(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000363:
				update_site5_11223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000427:
				update_site5_12223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000683:
				update_site5_22223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000351:
				update_site5_11133(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000367:
				update_site5_11233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000431:
				update_site5_12233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000687:
				update_site5_22233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000383:
				update_site5_11333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000447:
				update_site5_12333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000703:
				update_site5_22333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000511:
				update_site5_13333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5000767:
				update_site5_23333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_5_NEIGHB
			case 5001023:
				update_site5_33333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001365:
				update_site6_111111(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001366:
				update_site6_111112(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001370:
				update_site6_111122(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001386:
				update_site6_111222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001450:
				update_site6_112222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001706:
				update_site6_122222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002730:
				update_site6_222222(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001367:
				update_site6_111113(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001371:
				update_site6_111123(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001387:
				update_site6_111223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001451:
				update_site6_112223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001707:
				update_site6_122223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002731:
				update_site6_222223(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001375:
				update_site6_111133(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001391:
				update_site6_111233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001455:
				update_site6_112233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001711:
				update_site6_122233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002735:
				update_site6_222233(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001407:
				update_site6_111333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001471:
				update_site6_112333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001727:
				update_site6_122333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002751:
				update_site6_222333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001535:
				update_site6_113333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6001791:
				update_site6_123333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002815:
				update_site6_223333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6002047:
				update_site6_133333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6003071:
				update_site6_233333(site, sched[sweep]);
				break;
		#endif
		#ifdef USE_6_NEIGHB
			case 6004095:
				update_site6_333333(site, sched[sweep]);
				break;
		#endif
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
		return "algorithm: multi-spin, range-3 couplings, without fields";
	}
private:
	std::vector<site_type> sites;
	std::vector<sched_type> sched;

	std::mt19937 rgen;
	bgen_type bgen;

	// machine-generated code; do not edit
	#ifdef USE_1_NEIGHB
	void update_site1_1(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;

		if (r >= sched.r1) {
			word_type mask = b0;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_1_NEIGHB
	void update_site1_2(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;

		if (r >= sched.r2) {
			word_type mask = b0 | b1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_1_NEIGHB
	void update_site1_3(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;

		if (r >= sched.r3) {
			word_type mask = b0 | b1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_11(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;

		if (r >= sched.r2) {
			word_type mask = b0 | b1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_12(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;

		if (r >= sched.r1) {
			word_type mask = b1;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b0 | b1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_22(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r4) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_13(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_23(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_2_NEIGHB
	void update_site2_33(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_111(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;

		if (r >= sched.r1) {
			word_type mask = b1;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b0 | b1;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_112(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_122(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_113(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_123(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_133(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_3_NEIGHB
	void update_site3_333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1111(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1112(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1122(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_2222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1113(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1123(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_2223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1133(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_2233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_1333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_2333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_4_NEIGHB
	void update_site4_3333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11111(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11112(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11122(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_12222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_22222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11113(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11123(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_12223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_22223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11133(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_12233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_22233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_11333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_12333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_22333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_13333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_23333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r14) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_5_NEIGHB
	void update_site5_33333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r3) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r15) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111111(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r2) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111112(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;

		if (r >= sched.r1) {
			word_type mask = b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b0 & b1) | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b1 | b2;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b0 | b1) | b2;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111122(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_112222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_122222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_222222(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111113(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111123(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_112223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_122223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_222223(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111133(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_112233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_122233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_222233(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r14) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_111333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_112333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_122333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r14) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_222333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r15) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_113333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r2) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r14) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_123333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;

		if (r >= sched.r1) {
			word_type mask = b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r3) {
			word_type mask = ((b0 & b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (b1 & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b0 | b1) & b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r9) {
			word_type mask = b2 | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = ((b0 & b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = (b1 | b2) | b3;
			site.spin = spin ^ mask;
		} else if (r >= sched.r15) {
			word_type mask = ((b0 | b1) | b2) | b3;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_223333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;

		if (r >= sched.r2) {
			word_type mask = b3 | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (((b0 & b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r6) {
			word_type mask = ((b1 & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (((b0 | b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b2 | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = (((b0 & b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r16) {
			word_type mask = (((b0 | b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_133333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		b0 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;

		if (r >= sched.r2) {
			word_type mask = b3 | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r4) {
			word_type mask = (((b0 & b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r8) {
			word_type mask = (((b0 | b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r10) {
			word_type mask = (b2 | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r14) {
			word_type mask = ((b1 | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r16) {
			word_type mask = (((b0 | b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_233333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0;

		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;

		if (r >= sched.r1) {
			word_type mask = (((b0 | b1) | b2) & b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r5) {
			word_type mask = (((b0 & b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r7) {
			word_type mask = ((b1 & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r11) {
			word_type mask = (b2 | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r13) {
			word_type mask = (((b0 & b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r17) {
			word_type mask = (((b0 | b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

	// machine-generated code; do not edit
	#ifdef USE_6_NEIGHB
	void update_site6_333333(site_type& site, const sched_type& sched)
	{
		unsigned r = bgen() >> (word_size - depth);
		word_type spin = site.spin;

		word_type c, s;

		word_type b0 = 0, b1 = 0, b2 = 0, b3 = 0, b4 = 0;

		c = site.jzw0[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		b1 = s;
		c = site.jzw1[0] ^ (spin ^ sites[site.neighbs[0]].spin);
		s = b1 ^ c;
		b1 = s;
		c = site.jzw0[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw1[1] ^ (spin ^ sites[site.neighbs[1]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		b2 = s;
		c = site.jzw0[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[2] ^ (spin ^ sites[site.neighbs[2]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[3] ^ (spin ^ sites[site.neighbs[3]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw1[4] ^ (spin ^ sites[site.neighbs[4]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		b3 = s;
		c = site.jzw0[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b0 ^ c;
		c = b0 & c;
		b0 = s;
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;
		c = site.jzw1[5] ^ (spin ^ sites[site.neighbs[5]].spin);
		s = b1 ^ c;
		c = b1 & c;
		b1 = s;
		s = b2 ^ c;
		c = b2 & c;
		b2 = s;
		s = b3 ^ c;
		c = b3 & c;
		b3 = s;
		s = b4 ^ c;
		b4 = s;

		if (r >= sched.r6) {
			word_type mask = (((b0 & b1) & b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r12) {
			word_type mask = (b2 | b3) | b4;
			site.spin = spin ^ mask;
		} else if (r >= sched.r18) {
			word_type mask = (((b0 | b1) | b2) | b3) | b4;
			site.spin = spin ^ mask;
		} else {
			word_type mask = word_type(-1);
			site.spin = spin ^ mask;
		}
	}
	#endif

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
				h += site.hzv;

				en[offs + k] += h * spin;
			}
	}
};

#endif
