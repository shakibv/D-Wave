/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Implementation of single-spin simulated annealing algorithm for
Ising spin glasses on bipartite lattices with general interactions
with magnetic field and any number of neighbors.

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
#include <set>
#include <cassert>
#include <iterator>

#include "lattice.h"

#define OMP_VERSION_2

template<typename T = uint64_t>
  class Algorithm {

 public:

 typedef double value_type;
 typedef unsigned index_type;

 static const std::size_t word_size = 1;

 struct site_type{
   int spin;
   value_type hzv;
   std::vector<value_type> jzv;
   value_type de;
   index_type nneighbs;
   std::vector<index_type> neighbs;
 };

 typedef Lattice<value_type, index_type> lattice_type;

 Algorithm() {}

 template <typename SE>
 Algorithm(const lattice_type& lattice, const std::vector<SE>& sched0)
 : generator(41)
 {
   std::vector<site_type> sites0;
   lattice.init_sites(sites0);

   std::set<index_type> bin0, bin1;

   for(index_type s0 = 0; s0 < sites0.size(); ++s0)
     for(index_type k = 0; k < sites0[s0].nneighbs; ++k){

       const index_type s1 = sites0[s0].neighbs[k];

       if(bin0.find(s1) != bin0.end())
	 bin1.insert(s0);
       else if(bin0.find(s0) != bin0.end())
	 bin1.insert(s1);
       else if(bin1.find(s1) != bin1.end())
	 bin0.insert(s0);
       else if(bin1.find(s0) != bin1.end())
	 bin0.insert(s1);
       else{
	 bin0.insert(s0);
	 bin1.insert(s1);
       }

     }
    
   assert(bin0.size() + bin1.size() == sites0.size());

   if(bin1.size() < bin0.size())
     std::swap(bin0,bin1);

   const index_type ind0 = sites0.size();
   const index_type ind1 = sites0.size() + 1;

   max_edge = 0.0;
   for(const auto& site : sites0)
     max_edge += std::fabs(site.hzv);

   site_type sf0;
   for(const auto s1 : bin1){
     sf0.jzv.push_back(sites0[s1].hzv);
     sf0.neighbs.push_back(s1);
     ++sf0.nneighbs;

     sites0[s1].jzv.push_back(sites0[s1].hzv);
     sites0[s1].neighbs.push_back(ind0);
     ++sites0[s1].nneighbs;
   }

   site_type sf1;
   for(const auto s0 : bin0){
     sf1.jzv.push_back(sites0[s0].hzv);
     sf1.neighbs.push_back(s0);
     ++sf1.nneighbs;

     sites0[s0].jzv.push_back(sites0[s0].hzv);
     sites0[s0].neighbs.push_back(ind1);
     ++sites0[s0].nneighbs;
   }

   sf0.jzv.push_back(-max_edge);
   sf0.neighbs.push_back(ind1);
   ++sf0.nneighbs;

   sf1.jzv.push_back(-max_edge);
   sf1.neighbs.push_back(ind0);
   ++sf1.nneighbs;

   sites0.push_back(sf0);
   sites0.push_back(sf1);

   bin0.insert(ind0);
   bin1.insert(ind1);

   for(auto& site : sites0)
     site.hzv = 0;

   sites.resize(bin0.size());
   sums.resize(bin1.size());

   for(index_type s0 = 0; s0 < sites0.size(); ++s0)
     if(bin0.find(s0) != bin0.end()){

       const index_type ind = std::distance(bin0.begin(),bin0.find(s0));

       sites[ind].nneighbs = sites0[s0].nneighbs;
       for(index_type k = 0; k < sites0[s0].nneighbs; ++k){
	 sites[ind].neighbs.push_back(std::distance(bin1.begin(),bin1.find(sites0[s0].neighbs[k])));
	 sites[ind].jzv.push_back(sites0[s0].jzv[k]);
       }

     }

   rng_ = std::bind(std::uniform_real_distribution<double>(0, 1), std::ref(generator));

   bound_array.resize(sched0.size());

   auto ba = bound_array.begin();
   for(const auto& s : sched0){

     ba->resize(sites.size());
     for(auto& a : *ba)
       a = -std::log(rng_()) / s.beta;

     ba++;
   }

 }

 void reset_sites(const std::size_t rep)
 {
   generator.seed(rep+1);

   for(auto& site : sites)
     site.spin = 2 * ((generator() >> 29) & 1) - 1;

   std::fill(sums.begin(),sums.end(),0.0);
   for(const auto& site : sites)
     for(index_type k = 0; k < site.nneighbs; ++k)
       sums[site.neighbs[k]] += site.jzv[k] * site.spin;
 }

 value_type get_de(const site_type& site) const
 {
   value_type de = 0;
   for(index_type k = 0; k < site.nneighbs; ++k)
     de += std::fabs(sums[site.neighbs[k]]) - std::fabs(sums[site.neighbs[k]] - 2 * site.jzv[k] * site.spin);

   return de;
 }

 void flip_spin(site_type& site)
 {
   site.spin = -site.spin;
   for(index_type k = 0; k < site.nneighbs; ++k)
     sums[site.neighbs[k]] += 2 * site.jzv[k] * site.spin;
 }

 void do_sweep(const std::size_t sweep)
 {
   const std::size_t l = generator() % sites.size();
   const auto& ba = bound_array[sweep];

   for(std::size_t i = 0; i<l; ++i)
     if(get_de(sites[i]) < ba[i + sites.size() - l])
       flip_spin(sites[i]);

   for(std::size_t i = l; i<sites.size(); ++i)
     if(get_de(sites[i]) < ba[i - l])
       flip_spin(sites[i]);
 }

 std::size_t get_energies(std::vector<value_type>& en, const std::size_t offs) const
 {
   std::vector<value_type> sums0(sums.size(),0.0);

   for(const auto& site : sites)
     for(index_type k = 0; k < site.nneighbs; ++k)
       sums0[site.neighbs[k]] += site.jzv[k] * site.spin;

   value_type energy = 0.0;
   for(const auto& sum : sums0)
     energy -= std::fabs(sum);

   energy += max_edge;

   en[offs] = energy;
   return offs+1;
 }

 std::string get_info() const {return "algorithm: single-spin bipartite, variable degree";}

 private:

 std::vector<site_type> sites;
 std::vector<value_type> sums;
 std::vector<std::vector<double> > bound_array;

 std::mt19937 generator;
 std::function<double()> rng_; 

 value_type max_edge;
};

#endif
