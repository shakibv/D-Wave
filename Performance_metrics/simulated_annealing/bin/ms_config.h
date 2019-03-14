/******************************************************************************

Simulated annealing codes 
v1.0

---------------------------------------------------------------------

Defines flags that specify the number of neighbors to use in
multi-spin codes. Contains a function that checks whether the
specified number of neighbors is implemented.

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

#ifndef __MS_CONFIG_H__
#define __MS_CONFIG_H__

#define USE_1_NEIGHB
#define USE_2_NEIGHB
#define USE_3_NEIGHB
#define USE_4_NEIGHB
#define USE_5_NEIGHB
#define USE_6_NEIGHB

inline bool check_number_of_neighbors(unsigned nneighbs)
{
	switch (nneighbs) {
	case 1:
#ifdef USE_1_NEIGHB
		return true;
#else
		return false;
#endif
	case 2:
#ifdef USE_2_NEIGHB
		return true;
#else
		return false;
#endif
	case 3:
#ifdef USE_3_NEIGHB
		return true;
#else
		return false;
#endif
	case 4:
#ifdef USE_4_NEIGHB
		return true;
#else
		return false;
#endif
	case 5:
#ifdef USE_5_NEIGHB
		return true;
#else
		return false;
#endif
	case 6:
#ifdef USE_6_NEIGHB
		return true;
#else
		return false;
#endif
	default:
		return false;
	}
}

#endif
