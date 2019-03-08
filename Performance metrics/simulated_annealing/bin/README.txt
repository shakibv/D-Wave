=====================================================================
Simulated annealing codes 
v1.0
---------------------------------------------------------------------

Copyright (C) 2012-2013 by Sergei Isakov <isakov@itp.phys.ethz.ch>
                           Ilia Zintchenko <zintchenko@itp.phys.ethz.ch>

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

=====================================================================

---------------------------------------------------------------------
DESCRIPTION
---------------------------------------------------------------------

We present several efficient implementations of the simulated
annealing algorithm for Ising spin glasses on sparse graphs. In
particular, we provide a generic code for any choice of couplings, an
optimized code for bipartite graphs, and highly optimized
implementations using multi-spin coding for graphs with small maximum
degree and discrete couplings with a finite range. The latter codes
achieve up to 50 spin flips per nanosecond on modern Intel CPUs.

---------------------------------------------------------------------
REQUIREMENTS
---------------------------------------------------------------------

C++-11 conforming compiler
make

---------------------------------------------------------------------
INSTALLATION
---------------------------------------------------------------------

The codes can be built as follows: 

make <target>

where <target> specifies the code to build. To build a multi-threaded
version append <_omp> to target. Available targets are as follows

an_ms_r1_nf           Multi-spin code for range-1 interactions without magnetic field (approach one)

an_ms_r1_fi           Multi-spin code for range-1 interactions with magnetic field (approach one)

an_ms_r3_nf           Multi-spin code for range-3 interactions without magnetic field (approach one)

an_ms_r1_nf_v0        Multi-spin code for range-1 interactions without magnetic field (approach two)

an_ss_ge_fi           Single-spin code for general interactions with magnetic field (fixed number of neighbors)

an_ss_ge_fi_vdeg      Single-spin code for general interactions with magnetic field (any number of neighbors)

an_ss_ge_nf_bp        Single-spin code for general interactions on bipartite lattices without magnetic field (fixed number of neighbors)

an_ss_ge_nf_bp_vdeg   Single-spin code for general interactions on bipartite lattices without magnetic field (any number of neighbors)

an_ss_rn_fi           Single-spin code for range-n interactions with magnetic field (fixed number of neighbors)

an_ss_rn_fi_vdeg      Single-spin code for range-n interactions with magnetic field (any number of neighbors)


---------------------------------------------------------------------
USAGE
---------------------------------------------------------------------

Either of the previously described algorithms can be ran by using the
appropriate executable. Every algorithm follows the same command-line
interface which allows one to specify the lattice using -l, the
schedule -sched, the number of sweeps -s and the number of repetitions
-r. Because repetitions are independent of each other, the codes can
be trivially parallelized. The number of threads to run in parallel
can be specified by -t. In cases where one uses one of the
preprogrammed schedules, lin or exp, one can also specify the initial
inverse temperature beta0 and the final inverse temperature beta1. A
custom schedule can be loaded by using -sched followed by the name of
a text file. The full set of command-line arguments is as follows


-l [instance]         [instance] specifies instance file
-s [sweeps]           [sweeps] is number of sweeps                      
-r [reps]             [reps] is number of repetitions
-b0 [beta0]           [beta0] is initial inverse temperature. Default value: 0.1
-b1 [beta1]           [beta1] is final inverse temperature. Default value: 3.0
-r0 [rep0]            [rep0] is starting repetition. Default value: 0
-v                    if -v is set, timing and some other info is printed. Default value: not set
-g                    if -g is set, only the lowest energy solution is printed. Default value: not set
-sched [schedule]     [schedule] specifies a schedule. It can either be lin, exp or be a text file on the system which contains an inverse temperature on every line. Default value: lin
-t [threads]          [threads] is the number of threads to run in parallel. Default value: OMP NUM THREADS

The input lattice files are plain text files with following structure:
First line is the name of the lattice, and following N + M lines
contain N couplings and M local fields (not ordered). Each line
contains three values i, j and c. If i = j the line specifies a local
field on site i of size h_i = c. Otherwise, the line denotes a
coupling between spin i and j of value J_ij =c.

---------------------------------------------------------------------
SAMPLE INSTANCES AND RUNS
---------------------------------------------------------------------

The "example" sub-directory contains sample instances and sample input
and output data:

instance.txt          Ising spin glass instance on a 16x16 square lattice
instance126.txt       Ising spin glass instance on a 126-spin chimera graph
instance503.txt       Ising spin glass instance on a 503-spin chimera graph
sample_runs.txt       sample input and output data.