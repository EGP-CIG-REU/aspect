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
#include <aspect/utilities.h>
#include <aspect/vofinterface/vof_system.h>
#include <aspect/vofinterface/vof_utils.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        VoFSystem<dim>::VoFSystem (const FiniteElement<dim> &finite_element,
                                   const FiniteElement<dim> &vof_element,
                                   const Mapping<dim>       &mapping,
                                   const Quadrature<dim>    &quadrature,
                                   const Quadrature<dim-1>  &face_quadrature)
          :
          finite_element_values (mapping,
                                 finite_element, quadrature,
                                 update_values |
                                 update_gradients |
                                 update_JxW_values),
          face_finite_element_values (mapping,
                                      finite_element, face_quadrature,
                                      update_values |
                                      update_gradients |
                                      update_normal_vectors |
                                      update_JxW_values),
          neighbor_face_finite_element_values (mapping,
                                               finite_element, face_quadrature,
                                               update_values |
                                               update_gradients |
                                               update_normal_vectors |
                                               update_JxW_values),
          subface_finite_element_values (mapping,
                                         finite_element, face_quadrature,
                                         update_values |
                                         update_gradients |
                                         update_normal_vectors |
                                         update_JxW_values),
          local_dof_indices(finite_element.dofs_per_cell),
          phi_field (vof_element.dofs_per_cell, Utilities::signaling_nan<double>()),
          old_field_values (quadrature.size(), Utilities::signaling_nan<double>()),
          cell_i_n_values (quadrature.size(), Utilities::signaling_nan<Tensor<1, dim> > ()),
          cell_i_d_values (quadrature.size(), Utilities::signaling_nan<double> ()),
          face_current_velocity_values (face_quadrature.size(), Utilities::signaling_nan<Tensor<1, dim> >()),
          face_old_velocity_values (face_quadrature.size(), Utilities::signaling_nan<Tensor<1, dim> >()),
          face_old_old_velocity_values (face_quadrature.size(), Utilities::signaling_nan<Tensor<1, dim> >()),
          face_mesh_velocity_values (face_quadrature.size(), Utilities::signaling_nan<Tensor<1, dim> >())
        {}

        template <int dim>
        VoFSystem<dim>::VoFSystem (const VoFSystem &scratch)
          :
          finite_element_values (scratch.finite_element_values.get_mapping(),
                                 scratch.finite_element_values.get_fe(),
                                 scratch.finite_element_values.get_quadrature(),
                                 scratch.finite_element_values.get_update_flags()),
          face_finite_element_values (scratch.face_finite_element_values.get_mapping(),
                                      scratch.face_finite_element_values.get_fe(),
                                      scratch.face_finite_element_values.get_quadrature(),
                                      scratch.face_finite_element_values.get_update_flags()),
          neighbor_face_finite_element_values (scratch.neighbor_face_finite_element_values.get_mapping(),
                                               scratch.neighbor_face_finite_element_values.get_fe(),
                                               scratch.neighbor_face_finite_element_values.get_quadrature(),
                                               scratch.neighbor_face_finite_element_values.get_update_flags()),
          subface_finite_element_values (scratch.subface_finite_element_values.get_mapping(),
                                         scratch.subface_finite_element_values.get_fe(),
                                         scratch.subface_finite_element_values.get_quadrature(),
                                         scratch.subface_finite_element_values.get_update_flags()),
          local_dof_indices (scratch.finite_element_values.get_fe().dofs_per_cell),
          phi_field (scratch.phi_field),
          old_field_values (scratch.old_field_values),
          cell_i_n_values (scratch.cell_i_n_values),
          cell_i_d_values (scratch.cell_i_d_values),
          face_current_velocity_values (scratch.face_current_velocity_values),
          face_old_velocity_values (scratch.face_old_velocity_values),
          face_old_old_velocity_values (scratch.face_old_old_velocity_values),
          face_mesh_velocity_values (scratch.face_mesh_velocity_values)
        {}
      }

      namespace CopyData
      {
        /**
         * Constructor.
         * @param finite_element The element that describes the field for which we
         *    are trying to assemble a linear system. <b>Not</b> the global finite
         *    element.
         */
        template <int dim>
        VoFSystem<dim>::VoFSystem(const FiniteElement<dim> &finite_element)
          :
          local_matrix (finite_element.dofs_per_cell,
                        finite_element.dofs_per_cell),
          local_rhs (finite_element.dofs_per_cell),
          local_f_rhs (GeometryInfo<dim>::max_children_per_face *GeometryInfo<dim>::faces_per_cell,
                       Vector<double>(finite_element.dofs_per_cell)),
          assembled_rhs (GeometryInfo<dim>::max_children_per_face *GeometryInfo<dim>::faces_per_cell,
                         false),
          local_dof_indices (finite_element.dofs_per_cell),
          neighbor_dof_indices (GeometryInfo<dim>::max_children_per_face *GeometryInfo<dim>::faces_per_cell,
                                std::vector<types::global_dof_index>(finite_element.dofs_per_cell))
        {}

        template<int dim>
        VoFSystem<dim>::VoFSystem(const VoFSystem &data)
          :
          local_matrix (data.local_matrix),
          local_rhs (data.local_rhs),
          local_f_rhs (data.local_f_rhs),
          assembled_rhs (data.assembled_rhs),
          local_dof_indices (data.local_dof_indices),
          neighbor_dof_indices (data.neighbor_dof_indices)
        {}
      }
    }
  }

  template <int dim>
  void Simulator<dim>::do_vof_update ()
  {
    const unsigned int vof_block_idx = introspection.variable("vofs").block_index;
    const unsigned int vofN_block_idx = introspection.variable("vofsLS").block_index;

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
        solution.block(vofN_block_idx) = current_linearization_point.block(vofN_block_idx);
        update_vof_normals (solution);

        current_linearization_point.block(vof_block_idx) = solution(vof_block_idx);
        current_linearization_point.block(vofN_block_idx) = solution(vofN_block_idx);

        update_from_old = false;
      }
    // change dimension iteration order
    vof_dir_order_dsc = !vof_dir_order_dsc;
  }

  template <int dim>
  void Simulator<dim>::assemble_vof_system (unsigned int dir, bool update_from_old)
  {
    computing_timer.enter_section ("   Assemble VoF system");
    const unsigned int block_idx = introspection.variable("vofs").block_index;
    system_matrix.block(block_idx, block_idx) = 0;
    system_rhs = 0;

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    // const unsigned int vof_base_element = introspection.variable("vofs").base_index;
    const FiniteElement<dim> &vof_fe = (*introspection.variable("vofs").fe);

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     dof_handler.end()),
         std_cxx11::bind (&Simulator<dim>::
                          local_assemble_vof_system,
                          this,
                          dir,
                          update_from_old,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&Simulator<dim>::
                          copy_local_to_global_vof_system,
                          this,
                          std_cxx11::_1),
         // we have to assemble the term u.grad phi_i * phi_j, which is
         // of total polynomial degree
         //   stokes_deg - 1
         // (or similar for comp_deg). this suggests using a Gauss
         // quadrature formula of order
         //   stokes_deg/2
         // rounded up. do so. (note that x/2 rounded up
         // equals (x+1)/2 using integer division.)
         //
         // (note: we need to get at the advection element in
         // use for the scratch and copy objects below. the
         // base element for the compositional fields exists
         // only once, with multiplicity, so only query
         // introspection.block_indices.compositional_fields[0]
         // instead of subscripting with the correct compositional
         // field index.)
         internal::Assembly::Scratch::
         VoFSystem<dim> (finite_element,
                         vof_fe,
                         mapping,
                         QGauss<dim>((parameters.stokes_velocity_degree+1)/2),
                         QGauss<dim-1>((parameters.stokes_velocity_degree+1)/2)),
         internal::Assembly::CopyData::
         VoFSystem<dim> (vof_fe));

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);

    computing_timer.exit_section ();
  }

  template <int dim>
  void Simulator<dim>::local_assemble_vof_system (const unsigned int calc_dir,
                                                  bool update_from_old,
                                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                                  internal::Assembly::Scratch::VoFSystem<dim> &scratch,
                                                  internal::Assembly::CopyData::VoFSystem<dim> &data)
  {
    const bool old_velocity_avail = (timestep_number > 0);
    const bool old_old_velocity_avail = (timestep_number > 1);

    const unsigned int n_q_points    = scratch.finite_element_values.n_quadrature_points;
    const unsigned int n_f_q_points    = scratch.face_finite_element_values.n_quadrature_points;

    // also have the number of dofs that correspond just to the element for
    // the system we are currently trying to assemble
    const unsigned int vof_dofs_per_cell = data.local_dof_indices.size();

    Assert (vof_dofs_per_cell < scratch.finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
    Assert (vof_dofs_per_cell < scratch.face_finite_element_values.get_fe().dofs_per_cell, ExcInternalError());
    Assert (scratch.phi_field.size() == vof_dofs_per_cell, ExcInternalError());

    const FiniteElement<dim> &main_fe = scratch.finite_element_values.get_fe();

    const unsigned int vofN_component = introspection.variable("vofsN").first_component_index;
    const FEValuesExtractors::Vector vofN_n = FEValuesExtractors::Vector(vofN_component);
    const FEValuesExtractors::Scalar vofN_d = FEValuesExtractors::Scalar(vofN_component+dim);

    const unsigned int solution_component = introspection.variable("vofs").first_component_index;
    const FEValuesExtractors::Scalar solution_field = introspection.variable("vofs").extractor_scalar();
    const Quadrature<dim> &quadrature = scratch.finite_element_values.get_quadrature();

    //loop over all possible subfaces of the cell, and reset corresponding rhs
    for (unsigned int f = 0; f < GeometryInfo<dim>::max_children_per_face * GeometryInfo<dim>::faces_per_cell; ++f)
      {
        data.local_f_rhs[f] = 0;
        data.assembled_rhs[f] = false;
      }

    scratch.finite_element_values.reinit (cell);

    cell->get_dof_indices (scratch.local_dof_indices);
    for (unsigned int i=0; i<vof_dofs_per_cell; ++i)
      data.local_dof_indices[i] = scratch.local_dof_indices[main_fe.component_to_system_index(solution_component, i)];

    data.local_matrix = 0;
    data.local_rhs = 0;

    // Interface reconstruction data
    Tensor<1, dim, double> cell_i_normal;
    double cell_i_d = 0;

    if (update_from_old)
      {
        scratch.finite_element_values[solution_field].get_function_values (old_solution,
                                                                           scratch.old_field_values);

        scratch.finite_element_values[vofN_n].get_function_values (old_solution,
                                                                   scratch.cell_i_n_values);

        scratch.finite_element_values[vofN_d].get_function_values (old_solution,
                                                                   scratch.cell_i_d_values);
      }
    else
      {
        scratch.finite_element_values[solution_field].get_function_values (solution,
                                                                           scratch.old_field_values);

        scratch.finite_element_values[vofN_n].get_function_values (solution,
                                                                   scratch.cell_i_n_values);

        scratch.finite_element_values[vofN_d].get_function_values (solution,
                                                                   scratch.cell_i_d_values);
      }

    // interface values are constants, so can set from first value
    cell_i_normal = scratch.cell_i_n_values[0];
    cell_i_d = scratch.cell_i_d_values[0];

    const double cell_vol = cell->measure();
    // Obtain approximation to local interface
    for (unsigned int q = 0; q< n_q_points; ++q)
      {
        // Init FE field vals
        for (unsigned int k=0; k<vof_dofs_per_cell; ++k)
          scratch.phi_field[k] = scratch.finite_element_values[solution_field].value(main_fe.component_to_system_index(solution_component, k), q);

        // Init required local time
        for (unsigned int i = 0; i<vof_dofs_per_cell; ++i)
          {
            data.local_rhs[i] += scratch.old_field_values[q] *
                                 scratch.finite_element_values.JxW(q);
            for (unsigned int j=0; j<vof_dofs_per_cell; ++j)
              data.local_matrix (i, j) += scratch.phi_field[i] *
                                          scratch.phi_field[j] *
                                          scratch.finite_element_values.JxW(q);
          }
      }

    double face_flux;
    double dflux = 0;

    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      {
        unsigned int f_dim = f/2; // Obtain dimension
        bool f_dir_pos = (f%2==1);

        if (f_dim != calc_dir)
          continue;

        typename DoFHandler<dim>::face_iterator face = cell->face (f);

        scratch.face_finite_element_values.reinit (cell, f);

        scratch.face_finite_element_values[introspection.extractors.velocities]
        .get_function_values (current_linearization_point,
                              scratch.face_current_velocity_values);

        scratch.face_finite_element_values[introspection.extractors.velocities]
        .get_function_values (old_solution,
                              scratch.face_old_velocity_values);

        //scratch.face_finite_element_values[introspection.extractors.velocities]
        //.get_function_values (old_old_solution,
        //scratch.face_old_old_velocity_values);

        if (parameters.free_surface_enabled)
          scratch.face_finite_element_values[introspection.extractors.velocities]
          .get_function_values (free_surface->mesh_velocity,
                                scratch.face_mesh_velocity_values);

        face_flux = 0;
        double face_ls_d = 0;
        double face_ls_time_grad = 0;

        // Using VoF so need to accumulate flux through face
        for (unsigned int q=0; q<n_f_q_points; ++q)
          {

            Tensor<1,dim> current_u = scratch.face_current_velocity_values[q];

            //If old velocity available average to half timestep
            if (old_velocity_avail)
              current_u += 0.5*(scratch.face_old_velocity_values[q] -
                                scratch.face_current_velocity_values[q]);

            //Subtract off the mesh velocity for ALE corrections if necessary
            if (parameters.free_surface_enabled)
              current_u -= scratch.face_mesh_velocity_values[q];

            face_flux += time_step *
                         current_u *
                         scratch.face_finite_element_values.normal_vector(q) *
                         scratch.face_finite_element_values.JxW(q);

          }

        // Due to inability to reference this cell's values at the interface,
        // need to do explicit calculation
        if (f_dir_pos)
          {
            face_ls_d = cell_i_d - 0.5*cell_i_normal[f_dim];
            face_ls_time_grad = (face_flux/cell_vol)*cell_i_normal[f_dim];
          }
        else
          {
            face_ls_d = cell_i_d + 0.5*cell_i_normal[f_dim];
            face_ls_time_grad = -(face_flux/cell_vol)*cell_i_normal[f_dim];
          }

        // Calculate outward flux
        double flux_vof;
        if (face_flux < 0)
          {
            flux_vof = 0;
          }
        else
          {
            flux_vof = InterfaceTracker::calc_vof_flux_edge<dim> (f_dim,
                                                                  face_ls_time_grad,
                                                                  cell_i_normal,
                                                                  face_ls_d);
          }

        // Add fluxes to RHS
        for (unsigned int i=0; i<vof_dofs_per_cell; ++i)
          {
            data.local_rhs[i] -= flux_vof * face_flux;
          }

        if (face->at_boundary())
          {
            //TODO: Handle non-zero inflow VoF boundary conditions
          }
        else
          {
            Assert (cell->neighbor(f).state() == IteratorState::valid,
                    ExcInternalError());

            const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(f); //note: NOT active_cell_iterator, so this includes cells that are refined.
            if (!(face->has_children()))
              {
                Assert (cell->is_locally_owned(), ExcInternalError());
                //cell and neighbor are equal-sized, and cell has been chosen to assemble this face, so calculate from cell

                //how does the neighbor talk about this cell?
                const unsigned int neighbor2=cell->neighbor_of_neighbor(f);

                //set up neighbor values
                scratch.neighbor_face_finite_element_values.reinit (neighbor, neighbor2);

                std::vector<types::global_dof_index> neighbor_dof_indices (main_fe.dofs_per_cell);
                // get all dof indices on the neighbor, then extract those
                // that correspond to the solution_field we are interested in
                neighbor->get_dof_indices (neighbor_dof_indices);

                const unsigned int f_rhs_ind = f * GeometryInfo<dim>::max_children_per_face;

                for (unsigned int i=0; i<vof_dofs_per_cell; ++i)
                  data.neighbor_dof_indices[f_rhs_ind][i]
                    = neighbor_dof_indices[main_fe.component_to_system_index(solution_component, i)];

                data.assembled_rhs[f_rhs_ind] = true;

                // Add outward flux to neighbor
                for (unsigned int i=0; i<vof_dofs_per_cell; ++i)
                  {
                    data.local_f_rhs[f_rhs_ind][i] += flux_vof * face_flux;
                  }
              }
            else //face->has_children(), so always assemble from here.
              {
                // TODO: handle adaptive mesh
                Assert(false, ExcNotImplemented());
              }
          }
      }
  }

  template <int dim>
  void Simulator<dim>::copy_local_to_global_vof_system (const internal::Assembly::CopyData::VoFSystem<dim> &data)
  {
    // copy entries into the global matrix. note that these local contributions
    // only correspond to the advection dofs, as assembled above
    current_constraints.distribute_local_to_global (data.local_matrix,
                                                    data.local_rhs,
                                                    data.local_dof_indices,
                                                    system_matrix,
                                                    system_rhs);

    /* In the following, we copy DG contributions element by element. This
     * is allowed since there are no constraints imposed on discontinuous fields.
     */
    for (unsigned int f=0; f<GeometryInfo<dim>::max_children_per_face
         * GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (data.assembled_rhs[f])
          {
            for (unsigned int j=0; j<data.neighbor_dof_indices[f].size(); ++j)
              {
                system_rhs(data.neighbor_dof_indices[f][j]) += data.local_f_rhs[f][j];
              }
          }
      }
  }

  template <int dim>
  void Simulator<dim>::solve_vof_system ()
  {
    double vof_solver_tolerance = parameters.vof_solver_tolerance;
    unsigned int block_idx = introspection.variable("vofs").block_index;

    computing_timer.enter_section ("   Solve VoF system");
    pcout << "   Solving VoF system... " << std::flush;

    const double tolerance = std::max(1e-50,
                                      vof_solver_tolerance*system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

#ifdef ASPECT_USE_PETSC
    SolverCG<LinearAlgebra::Vector> solver(solver_control);
#else
    TrilinosWrappers::SolverCG solver(solver_control);
#endif

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a ConstraintMatrix
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      introspection.index_sets.system_partitioning,
      mpi_communicator);
    distributed_solution.block(block_idx) = current_linearization_point.block (block_idx);

    current_constraints.set_zero(distributed_solution);

    // solve the linear system:
    try
      {
#ifdef ASPECT_USE_PETSC
        solver.solve (system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      system_rhs.block(block_idx),
                      PreconditionIdentity());
#else
        solver.solve (system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      system_rhs.block(block_idx),
                      TrilinosWrappers::PreconditionIdentity());
#endif
      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          AssertThrow (false,
                       ExcMessage (std::string("The iterative advection solver "
                                               "did not converge. It reported the following error:\n\n")
                                   +
                                   exc.what()))
          else
            throw QuietException();
      }

    current_constraints.distribute (distributed_solution);
    solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    pcout << solver_control.last_step()
          << " iterations." << std::endl;

    // Do not add VoF solver iterations to statistics, duplicaiton due to
    // splitting messes with file format
    // statistics.add_value("Iterations for VoF solver",
    //                      solver_control.last_step());

    computing_timer.exit_section();
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::do_vof_update (); \
  template void Simulator<dim>::assemble_vof_system (unsigned int dir, \
                                                     bool update_from_old); \
  template void Simulator<dim>::local_assemble_vof_system (const unsigned int calc_dir, \
                                                           bool update_from_old, \
                                                           const typename DoFHandler<dim>::active_cell_iterator &cell, \
                                                           internal::Assembly::Scratch::VoFSystem<dim> &scratch, \
                                                           internal::Assembly::CopyData::VoFSystem<dim> &data); \
  template void Simulator<dim>::copy_local_to_global_vof_system (const internal::Assembly::CopyData::VoFSystem<dim> &data); \
  template void Simulator<dim>::solve_vof_system ();


  ASPECT_INSTANTIATE(INSTANTIATE)
}
