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
  VoFField<dim>::VoFField(const FEVariable<dim> &fraction,
                          const FEVariable<dim> &reconstruction,
                          const FEVariable<dim> &level_set,
                          const std::string c_field_name)
    : fraction (fraction),
      reconstruction (reconstruction),
      level_set (level_set),
      c_field_name (c_field_name)
  {}

  template <int dim>
  VoFHandler<dim>::VoFHandler (Simulator<dim> &simulator,
                               ParameterHandler &prm)
    : sim (simulator),
      vof_initial_conditions (VoFInitialConditions::create_initial_conditions<dim>(prm))
  {
    this->initialize_simulator(sim);
    parse_parameters (prm);

    sim.signals.edit_finite_element_variables.connect(std_cxx11::bind(&aspect::VoFHandler<dim>::edit_finite_element_variables,
                                                                      std_cxx11::ref(*this),
                                                                      std_cxx11::_1));
    sim.signals.post_set_initial_state.connect(std_cxx11::bind(&aspect::VoFHandler<dim>::set_initial_vofs,
                                                               std_cxx11::ref(*this)));
  }

  template <int dim>
  void
  VoFHandler<dim>::edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars)
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
  VoFHandler<dim>::declare_parameters (ParameterHandler &prm)
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

      prm.declare_entry ("VoF field names", "",
                         Patterns::List(Patterns::Anything()),
                         "User-defined names for VoF fields.");

      // TODO: Replace with Map
      prm.declare_entry ("VoF composition variable", "",
                         Patterns::List(Patterns::Anything()),
                         "Name of compositional field to write VoF composition to.");
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void
  VoFHandler<dim>::parse_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("VoF config");
    {
      vof_epsilon = prm.get_double("Small volume");

      vof_solver_tolerance = prm.get_double("VoF solver tolerance");

      n_vof_fields = 1;

      vof_field_names = Utilities::split_string_list (prm.get("VoF field names"));
      AssertThrow((vof_field_names.size() == 0) ||
                  (vof_field_names.size() == n_vof_fields),
                  ExcMessage("The length of the list of names for the VoF fields "
                             "needs to either be empty or have length equal to the "
                             "number of compositional fields."));

      // check that names use only allowed characters, are not empty strings, and are unique
      for (unsigned int i=0; i<vof_field_names.size(); ++i)
        {
          Assert (vof_field_names[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                       "0123456789_") == std::string::npos,
                  ExcMessage("Invalid character in field " + vof_field_names[i] + ". "
                             "Names of VoF fields should consist of a "
                             "combination of letters, numbers and underscores."));
          Assert (vof_field_names[i].size() > 0,
                  ExcMessage("Invalid name of field " + vof_field_names[i] + ". "
                             "Names of VoF fields need to be non-empty."));
          for (unsigned int j=0; j<i; ++j)
            Assert (vof_field_names[i] != vof_field_names[j],
                    ExcMessage("Names of VoF fields have to be unique! " + vof_field_names[i] +
                               " is used more than once."));
        }

      // default names if not empty
      if (vof_field_names.size()==0)
        {
          for (unsigned int i=0; i<n_vof_fields; ++i)
            vof_field_names.push_back("F_" + Utilities::int_to_string(i+1));
        }

      vof_composition_vars = Utilities::split_string_list (prm.get("VoF composition variable"));
      AssertThrow((vof_composition_vars.size() == 0) ||
                  (vof_composition_vars.size() == n_vof_fields),
                  ExcMessage("The length of the list of names for the VoF fields "
                             "needs to either be empty or have length equal to the "
                             "number of compositional fields."));

      // check that names use only allowed characters, are not empty strings, and are unique
      for (unsigned int i=0; i<vof_composition_vars.size(); ++i)
        {
          Assert (vof_composition_vars[i].find_first_not_of("abcdefghijklmnopqrstuvwxyz"
                                                            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                                            "0123456789_") == std::string::npos,
                  ExcMessage("Invalid character in field " + vof_composition_vars[i] + ". "
                             "Names of VoF fields should consist of a "
                             "combination of letters, numbers and underscores."));

          for (unsigned int j=0; j<i; ++j)
            Assert (vof_composition_vars[i] != vof_composition_vars[j],
                    ExcMessage("Names of VoF fields have to be unique! " + vof_composition_vars[i] +
                               " is used more than once."));
        }

      if (vof_composition_vars.size()>0)
        {
          if (!sim.parameters.use_discontinuous_composition_discretization)
            {
              AssertThrow(false, ExcMessage("VoF composition field not implemented for continuous composition."));
            }

          for (unsigned int f=0; f<vof_composition_vars.size(); ++f)
            {
              bool field_exists=false;

              for (unsigned int i=0; i<sim.parameters.n_compositional_fields; ++i)
                {
                  field_exists = field_exists ||
                                 (vof_composition_vars[f]==sim.parameters.names_of_compositional_fields[i]);
                }

              Assert(field_exists, ExcMessage("VoF composition field variable " +
                                              vof_composition_vars[f] +
                                              " does not exist."));
            }
        }

      // no connection if not empty
      if (vof_composition_vars.size()==0)
        {
          for (unsigned int i=0; i<n_vof_fields; ++i)
            vof_composition_vars.push_back("");
        }
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void
  VoFHandler<dim>::initialize (ParameterHandler &prm)
  {
    // Do checks on required assumptions
    AssertThrow(dim==2,ExcMessage("Volume of Fluid Interface Tracking is currently only functional for dim=2."));

    AssertThrow(!sim.material_model->is_compressible(), ExcMessage("Volume of Fluid Interface Tracking currently assumes incompressiblity."));

    AssertThrow(!this->get_parameters().free_surface_enabled,
                ExcMessage("Volume of Fluid Interface Tracking is currently incompatible with the Free Surface implementation."));

    // Check for correct mapping

    if ( this->get_parameters().initial_adaptive_refinement > 0 ||
         this->get_parameters().adaptive_refinement_interval > 0 )
      {
        // AMR active so check refinement strategy includes 'vof boundary'
      }

    // Gather data

    data.push_back(VoFField<dim>(sim.introspection.variable("vofs"),
                                 sim.introspection.variable("vofsN"),
                                 sim.introspection.variable("vofsLS"),
                                 vof_composition_vars[0]));

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
  const VoFField<dim> &VoFHandler<dim>::get_field() const
  {
    return data[0];
  }

  template <int dim>
  void VoFHandler<dim>::do_vof_update ()
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
            assemble_vof_system(data[0], dir, update_from_old);
          }
        else
          {
            assemble_vof_system(data[0], dim-dir-1, update_from_old);
          }
        solve_vof_system (data[0]);
        // Copy current candidate normals.
        // primarily useful for exact linear translation
        sim.solution.block(vofN_block_idx) = sim.old_solution.block(vofN_block_idx);
        update_vof_normals (data[0], sim.solution);

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
  template class VoFHandler<dim>;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
