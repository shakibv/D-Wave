/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of single-spin simulated annealing algorithm for
Ising spin glasses with general interactions with magnetic field
and fixed number of neighbors.

---------------------------------------------------------------------

Copyright (C) 2013 by Ilia Zintchenko <zintchenko@itp.phys.ethz.ch>

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
#include <functional>

#include "lattice.h"
#include "ss_config.h"

#define OMP_VERSION_2

template<typename T = uint64_t>
  class Algorithm
  {
  public:
	
  typedef double value_type;
  typedef unsigned index_type;

  static const unsigned MAXNB = MAX_NUM_NEIGHBORS;
  static const std::size_t word_size = 1;

  struct site_type{
    int spin;
    value_type hzv;
    value_type jzv[MAXNB];
    value_type de;
    index_type nneighbs;
    index_type neighbs[MAXNB];
  };

  typedef Lattice<value_type, index_type> lattice_type;

  Algorithm() {}

  template <typename SE>
  Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
  : generator(41)
  {
    lattice.init_sites(sites, MAXNB);

    for(auto& site : sites)
      for(index_type k = site.nneighbs; k < MAXNB; ++k){
	site.jzv[k] = 0;
	site.neighbs[k] = 0;
      }

    rng_ = std::bind(std::uniform_real_distribution<double>(0, 1), std::ref(generator));

    bound_array.resize(sched0.size());

    auto ba = bound_array.begin();
    for(const auto& s : sched0){

      ba->resize(sites.size());
      for(auto& a : *ba)
        a = -std::log(rng_()) / (s.beta * 2);

      ++ba;
    }

  }

  void reset_sites(const std::size_t rep)
  {
    generator.seed(rep+1);

    for(auto& site : sites)
      site.spin = 2 * ((generator() >> 29) & 1) - 1;

    for(auto& site : sites){
      value_type tmp = site.hzv;
      for(index_type k = 0; k < MAXNB; ++k)
	tmp += site.jzv[k] * sites[site.neighbs[k]].spin;
      site.de = -tmp * site.spin;
    }
  }  

  void flip_spin(site_type& site)
  {
    site.spin = -site.spin;
    site.de = -site.de;

    for (index_type k = 0; k<MAXNB; ++k) {
      site_type& neighbor = sites[site.neighbs[k]];
      neighbor.de -= 2 * neighbor.spin * site.jzv[k] * site.spin;
    }    
  }

  void do_sweep(const std::size_t sweep)
  {
    const std::size_t l = generator() % sites.size();
    const auto& ba = bound_array[sweep];

    for(std::size_t i = 0; i<l; ++i)
      if(sites[i].de<  ba[i + sites.size() - l])
        flip_spin(sites[i]);

    for(std::size_t i = l; i<sites.size(); ++i)
      if(sites[i].de < ba[i - l])
        flip_spin(sites[i]);

  } 

  std::size_t get_energies(std::vector<value_type>& en, const std::size_t offs) const
  {
    value_type energy = 0;
    for(const auto& site : sites){
      value_type tmp = site.hzv;
      for(index_type k = 0; k < site.nneighbs; ++k)
	tmp += sites[site.neighbs[k]].spin * site.jzv[k] / 2;

      energy += tmp * site.spin;
    }

    en[offs] = energy;
    return offs+1;
  }
 
  std::string get_info() const {return "algorithm: single-spin generic";}

  private:

  std::vector<site_type> sites;
  std::vector<std::vector<double> > bound_array;

  std::mt19937 generator;
  std::function<double()> rng_; 

  };

#endif
