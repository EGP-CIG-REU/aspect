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
#include <aspect/simulator_access.h>
#include <aspect/vof_initial_conditions/interface.h>
#include <aspect/vof/assembly.h>

using namespace dealii;

namespace aspect
{
  template<int dim>
  struct VoFField
  {
    /**
     * Initialize th
     */
    VoFField(const FEVariable<dim> &fraction,
             const FEVariable<dim> &reconstruction,
             const FEVariable<dim> &level_set,
             const std::string c_field_name);

    const FEVariable<dim> &fraction, &reconstruction, &level_set;

    const std::string c_field_name;
  };

  /**
   * A member class that isolates the functions and variables that deal
   * with the volume of fluid implementation. If Volume of Fluid interface
   * tracking is not active, there is no instantiation of this class at
   * all.
   */
  template <int dim>
  class VoFHandler : public SimulatorAccess<dim>
  {
    public:
      // Construtor
      VoFHandler(Simulator<dim> &sim, ParameterHandler &prm);

      void edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars);

      // Parameter handling
      static
      void declare_parameters (ParameterHandler &prm);

      void parse_parameters (ParameterHandler &prm);

      // initialiation
      void initialize (ParameterHandler &prm);

      // Functions for initialization of state
      void set_initial_vofs ();
      void init_vof_compos ();
      void init_vof_ls ();

      // Do interface reconstruction
      void update_vof_normals (LinearAlgebra::BlockVector &solution);

      // Logic to handle dimensionally split update
      void do_vof_update ();

      // Assembly
      void assemble_vof_system (unsigned int dir,
                                bool update_from_old);
      void local_assemble_vof_system (const unsigned int calc_dir,
                                      bool update_from_old,
                                      const typename DoFHandler<dim>::active_cell_iterator &cell,
                                      internal::Assembly::Scratch::VoFSystem<dim> &scratch,
                                      internal::Assembly::CopyData::VoFSystem<dim> &data);
      void copy_local_to_global_vof_system (const internal::Assembly::CopyData::VoFSystem<dim> &data);
      // Solver
      void solve_vof_system ();


    private:
      // Parent simulator
      Simulator<dim> &sim;

      //Initial conditions
      const std_cxx11::unique_ptr<VoFInitialConditions::Interface<dim> >      vof_initial_conditions;

      // Store
      VoFField<dim> *data;

      // Minimal considered volume fraction
      double vof_epsilon;

      double vof_solver_tolerance;

      std::string vof_composition_var;

      // Order for split update
      bool vof_dir_order_dsc;

      friend class Simulator<dim>;
      friend class SimulatorAccess<dim>;
  };

}

#endif
