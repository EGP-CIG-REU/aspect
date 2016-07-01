/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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

#ifndef __aspect_vof_handler_h
#define __aspect_vof_handler_h

#include <aspect/simulator.h>

using namespace dealii;

namespace aspect
{

  /**
   * A member class that isolates the functions and variables that deal
   * with the volume of fluid implementation. If Volume of Fluid interface
   * tracking is not active, there is no instantiation of this class at
   * all.
   */
  template <int dim>
  class Simulator<dim>::VoFHandler
  {
    public:
      // Construtor
      VoFHandler(Simulator<dim> &sim, ParameterHandler &prm);

    private:
      // Parent simulator
      Simulator<dim> &sim;

      // Order for split update
      bool vof_dir_order_dsc;

      friend class Simulator<dim>;
      friend class SimulatorAccess<dim>;
  };

}

#endif
