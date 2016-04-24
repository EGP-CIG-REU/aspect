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

#ifndef __aspect__vofinterface_VOFEngineSig_h
#define __aspect__vofinterface_VOFEngineSig_h

#include <deal.II/base/parameter_handler.h>

// Aspect includes
#include <aspect/global.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    // Parameter handling
    void declare_parameters (const unsigned int dim,
                             ParameterHandler &prm);
    template<int dim>
    void parse_parameters (const Parameters<dim> Parameters,
                           ParameterHandler &prm);
    void parameter_connector ();

    template <int dim>
    void add_vof_var(std::vector<VariableDeclaration<dim>> &vars);

    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals);
  }
}

#endif
