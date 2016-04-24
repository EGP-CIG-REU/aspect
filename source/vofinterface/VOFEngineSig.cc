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

#include <aspect/vofinterface/VOFEngineSig.h>

#include <deal.II/fe/fe_dgq.h>

namespace aspect
{
  namespace InterfaceTracker
  {
    using namespace dealii;

    // Parameters
    bool useVOFTracking;

    // Parameter Handling
    void declare_parameters (const unsigned int dim,
                             ParameterHandler &prm)
    {
      prm.declare_entry ("Use VoF tracking", "false",
                         Patterns::Bool (),
                         "When set to true, VoF fluid tracking will be used");
    }

    template <int dim>
    void parse_parameters (const Parameters<dim> Parameters,
                           ParameterHandler &prm)
    {
      useVOFTracking = prm.get_bool("Use VoF tracking");
    }

    void parameter_connector ()
    {
      SimulatorSignals<2>::declare_additional_parameters.connect(&declare_parameters);
      SimulatorSignals<3>::declare_additional_parameters.connect(&declare_parameters);

      SimulatorSignals<2>::parse_additional_parameters.connect(&parse_parameters<2>);
      SimulatorSignals<3>::parse_additional_parameters.connect(&parse_parameters<3>);
    }

    ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)

    template <int dim>
    void add_vof_var(std::vector<VariableDeclaration<dim>> &vars)
    {
      vars.push_back(VariableDeclaration<dim>("vofs",
                                              std_cxx11::shared_ptr<FiniteElement<dim>>(
                                                new FE_DGQ<dim>(0)),
                                              1,
                                              1));
    }

    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals)
    {
      // Return if VoF not intended to be active
      //
      if (!useVOFTracking)
        return;

      signals.edit_finite_element_variables.connect(&add_vof_var<dim>);
    }

    ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                      signal_connector<3>)
  }
}

