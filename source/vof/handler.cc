
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

#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/vof/handler.h>
#include <aspect/vof/assembly.h>

#include <deal.II/fe/fe_dgq.h>

using namespace dealii;

namespace aspect
{

  template <int dim>
  Simulator<dim>::VoFHandler::VoFHandler (Simulator<dim> &simulator,
                                          ParameterHandler &prm)
    : sim (simulator),
      vof_initial_conditions (VoFInitialConditions::create_initial_conditions<dim>(prm))
  {
    parse_parameters (prm);

    sim.signals.edit_finite_element_variables.connect(std_cxx11::bind(&aspect::Simulator<dim>::VoFHandler::edit_finite_element_variables,
                                                                      std_cxx11::ref(*this),
                                                                      std_cxx11::_1));
  }

  template <int dim>
  void
  Simulator<dim>::VoFHandler::edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars)
  {
    vars.push_back(VariableDeclaration<dim>("vofs",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(0)),
                                            1,
                                            1));

    vars.push_back(VariableDeclaration<dim>("vofsN",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(0)),
                                            dim+1,
                                            1));

    vars.push_back(VariableDeclaration<dim>("vofsLS",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(1)),
                                            1,
                                            1));
  }

  template <int dim>
  void
  Simulator<dim>::VoFHandler::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("VoF config");
    {
      prm.declare_entry ("Small volume", "1e-6",
                         Patterns::Double (0, 1),
                         "Minimum significant volume. VOFs below this considered to be zero.");

      prm.declare_entry ("VoF solver tolerance", "1e-12",
                         Patterns::Double(0,1),
                         "The relative tolerance up to which the linear system for "
                         "the VoF system gets solved. See 'linear solver "
                         "tolerance' for more details.");

      prm.declare_entry ("VoF composition variable", "",
                         Patterns::Anything(),
                         "Name of compositional field to write VoF composition to.");
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void
  Simulator<dim>::VoFHandler::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("VoF config");
    {
      vof_epsilon = prm.get_double("Small volume");

      vof_solver_tolerance = prm.get_double("VoF solver tolerance");

      vof_composition_var = prm.get("VoF composition variable");

      if (vof_composition_var!="")
        {
          if (!sim.parameters.use_discontinuous_composition_discretization)
            {
              Assert(false, ExcMessage("VoF composition field not implemented for continuous composition."));
            }

          bool field_exists=false;

          for (unsigned int i=0; i<sim.parameters.n_compositional_fields; ++i)
            {
              field_exists = field_exists ||
                             (vof_composition_var==sim.parameters.names_of_compositional_fields[i]);
            }

          Assert(field_exists, ExcMessage("VoF composition field variable " +
                                          vof_composition_var +
                                          " does not exist."));
        }
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void
  Simulator<dim>::VoFHandler::initialize (ParameterHandler &prm)
  {
    // Do checks on required assumptions


    // Do initial conditions setup
    if (SimulatorAccess<dim> *sim_a = dynamic_cast<SimulatorAccess<dim>*>(vof_initial_conditions.get()))
      sim_a->initialize_simulator (sim);
    if (vof_initial_conditions.get())
      {
        vof_initial_conditions->parse_parameters (prm);
        vof_initial_conditions->initialize ();
      }

  }

  template <int dim>
  void Simulator<dim>::VoFHandler::do_vof_update ()
  {
    const unsigned int vof_block_idx = sim.introspection.variable("vofs").block_index;
    const unsigned int vofN_block_idx = sim.introspection.variable("vofsLS").block_index;

    // Reset current base to values at beginning of timestep

    // Due to dimensionally split formulation, use strang splitting
    // TODO: Reformulate for unsplit (may require flux limiter)
    bool update_from_old = true;
    for (unsigned int dir = 0; dir < dim; ++dir)
      {
        // Update base to intermediate solution
        if (!vof_dir_order_dsc)
          {
            assemble_vof_system(dir, update_from_old);
          }
        else
          {
            assemble_vof_system(dim-dir-1, update_from_old);
          }
        solve_vof_system ();
        // Copy current candidate normals.
        // primarily useful for exact linear translation
        sim.solution.block(vofN_block_idx) = sim.old_solution.block(vofN_block_idx);
        update_vof_normals (sim.solution);

        sim.current_linearization_point.block(vof_block_idx) = sim.solution.block(vof_block_idx);
        sim.current_linearization_point.block(vofN_block_idx) = sim.solution.block(vofN_block_idx);
        update_from_old = false;
      }
    // change dimension iteration order
    vof_dir_order_dsc = !vof_dir_order_dsc;
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>::VoFHandler;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
