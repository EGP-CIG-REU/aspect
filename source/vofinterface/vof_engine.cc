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

#include <aspect/vofinterface/vof_engine.h>

#include <deal.II/fe/fe_dgq.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    // Add required variables

    template <int dim>
    void VOFEngine<dim>::add_vof_vars(const Parameters<dim> &parameters,
                                      std::vector<VariableDeclaration<dim>> &vars)
    {
      if (!parameters.vof_tracking_enabled)
        return;

      vars.push_back(VariableDeclaration<dim>("vofs",
                                              std_cxx11::shared_ptr<FiniteElement<dim>>(
                                                new FE_DGQ<dim>(0)),
                                              1,
                                              1));

      vars.push_back(VariableDeclaration<dim>("vofsLS",
                                              std_cxx11::shared_ptr<FiniteElement<dim>>(
                                                new FE_DGQ<dim>(1)),
                                              1,
                                              1));
    }

    // ASPECT registration code

    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals)
    {
      signals.edit_finite_element_variables.connect(VOFEngine<dim>::add_vof_vars);
    }

    ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                      signal_connector<3>)
  }
}

