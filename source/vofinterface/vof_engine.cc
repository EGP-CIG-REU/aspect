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


    // TODO: Find better way to store parameters than as general variables
    bool vof_tracking_enabled;

    // Handle general parameters
    //
    void gen_declare_parameters (const unsigned int dim,
                                 ParameterHandler &prm)
    {
      prm.declare_entry ("Use VoF tracking", "false",
                         Patterns::Bool (),
                         "When set to true, VoF interface tracking will be used");
    }

    template <int dim>
    void gen_parse_parameters (const Parameters<dim> &parameters,
                               ParameterHandler &prm)
    {
      vof_tracking_enabled = prm.get_bool("Use VoF tracking");
      if (vof_tracking_enabled)
        Assert(dim==2,ExcMessage("VoF interface tracking not implemented for dim>2."));
    }

    // Add required variables

    template <int dim>
    void VOFEngine<dim>::add_vof_vars(std::vector<VariableDeclaration<dim>> &vars)
    {
      if (!vof_tracking_enabled)
        return;

      vars.push_back(VariableDeclaration<dim>("vofs",
                                              std_cxx11::shared_ptr<FiniteElement<dim>>(
                                                new FE_DGQ<dim>(0)),
                                              1,
                                              1));
    }

    // Engine constructor

    template <int dim>
    VOFEngine<dim>::VOFEngine()
    {
    }

    // ASPECT registration code

    void parameter_connector ()
    {
      SimulatorSignals<2>::declare_additional_parameters
      .connect(&gen_declare_parameters);
      SimulatorSignals<3>::declare_additional_parameters
      .connect(&gen_declare_parameters);

      SimulatorSignals<2>::parse_additional_parameters
      .connect(gen_parse_parameters<2>);
      SimulatorSignals<3>::parse_additional_parameters
      .connect(&gen_parse_parameters<3>);
    }

    ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)

    template <int dim>
    void signal_connector (SimulatorSignals<dim> &signals)
    {
      signals.edit_finite_element_variables.connect(VOFEngine<dim>::add_vof_vars);
    }

    ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>,
                                      signal_connector<3>)
  }
}

