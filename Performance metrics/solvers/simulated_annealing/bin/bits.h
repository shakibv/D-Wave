/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Contains lagged Fibonacci and linear congruential random number
generators.

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

#ifndef __BITS_H__
#define __BITS_H__

#include <random>

template <typename G, typename T>
inline T random_word(G& rgen, const T&)
{
	T word = T(rgen());
	for (unsigned i = 1; i < sizeof(T) / 4; ++i)
		word |= T(rgen()) << i * 32;

	return word;
}

template <typename T = uint64_t, std::size_t j = 418, std::size_t k = 1279>
class bitgen_lagfib {
public:
	typedef T word_type;

	bitgen_lagfib()
	{
		seed(1);
	}

	bitgen_lagfib(word_type seed_)
	{
		seed(seed_);
	}

	void seed(word_type seed)
	{
		std::mt19937 rgen(seed);

		for (unsigned i = 0; i < k; ++i)
			fibbuf[i] = random_word(rgen, fibbuf[i]);

		p = j;
		o = 0;
	}

	word_type operator()()
	{
		word_type r = fibbuf[p] += fibbuf[o];

		if (++p >= k) p = 0;
		if (++o >= k) o = 0;

		return r;
	}
private:
	std::size_t p;
	std::size_t o;
	word_type fibbuf[k];
};

template <typename T = uint64_t,
	T a = 6364136223846793005ULL, T c = 1442695040888963407ULL>
class bitgen_lincon {
public:
	typedef T word_type;

	bitgen_lincon()
	{
		r = 1;
	}

	bitgen_lincon(word_type seed)
	{
		r = seed;
	}

	void seed(word_type seed)
	{
		r = seed;
	}

	word_type operator()()
	{
		return r = r * a + c;
	}
private:
	word_type r;
	unsigned dummy[128];
};

#endif
