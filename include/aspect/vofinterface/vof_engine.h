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

#ifndef __aspect__vofinterface_VoFEngineSig_h
#define __aspect__vofinterface_VoFEngineSig_h

#include <deal.II/base/parameter_handler.h>

// Aspect includes
#include <aspect/global.h>
#include <aspect/simulator_signals.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    template <int dim>
    class VoFEngine
    {
      public:

        static
        void add_vof_vars(const Parameters<dim> &parameters,
                          std::vector<VariableDeclaration<dim>> &vars);

      private:
    };

    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals);
  }
}

#endif
