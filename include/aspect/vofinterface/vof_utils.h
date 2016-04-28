/*
 Copyright (C) 2016 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file doc/COPYING.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef __aspect__vofinterface_vof_utils_h
#define __aspect__vofinterface_vof_utils_h

#include <aspect/global.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    // Utils

    template<int dim>
    double vof_from_d (const Tensor<1, dim, double> normal,
                       const double d);

    template<int dim>
    double d_from_vof (const Tensor<1, dim, double> normal,
                       const double vol);

    template<int dim>
    double calc_vof_flux_edge (const unsigned int dir,
                               const double timeGrad,
                               const Tensor<1, dim, double> cell_normal,
                               const double d_face);
  }
}

#endif
