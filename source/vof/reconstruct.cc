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

#include <aspect/global.h>
#include <aspect/simulator.h>
#include <aspect/vof/handler.h>
#include <aspect/vof/utilities.h>

namespace aspect
{
  using namespace dealii;

  template <>
  void VoFHandler<2>::update_vof_normals (const VoFField<2> field,
                                          LinearAlgebra::BlockVector &solution)
  {
    const int dim = 2;

    LinearAlgebra::BlockVector initial_solution;

    sim.computing_timer.enter_section("   Reconstruct VoF interfaces");

    initial_solution.reinit(sim.system_rhs, false);

    // Boundary reference
    typename DoFHandler<dim>::active_cell_iterator endc =
      sim.dof_handler.end ();

    // Interface Reconstruction vars
    const unsigned int n_local = 9;

    Vector<double> local_vofs (n_local);
    std::vector<Point<dim>> resc_cell_centers (n_local);

    const unsigned int n_sums = 3;
    std::vector<double> strip_sums (dim * n_sums);

    const unsigned int n_normals = 6+1;
    std::vector<Tensor<1, dim, double>> normals (n_normals);
    std::vector<double> errs (n_normals);

    // Normal holding vars
    Point<dim> uReCen;
    Tensor<1, dim, double> normal;
    double d;

    for (unsigned int i=0; i<dim; ++i)
      uReCen[i] = 0.5;

    std::vector<types::global_dof_index> cell_dof_indicies (sim.finite_element.dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indicies (sim.finite_element.dofs_per_cell);

    const FEVariable<dim> &vof_var = data->fraction;
    const unsigned int vof_c_index = vof_var.first_component_index;
    const unsigned int vof_ind
      = sim.finite_element.component_to_system_index(vof_c_index, 0);

    const FEVariable<dim> &vofN_var = data->reconstruction;
    const unsigned int vofN_c_index = vofN_var.first_component_index;
    const unsigned int vofN_blockidx = vofN_var.block_index;

    const FEVariable<dim> &vofLS_var = data->level_set;
    const unsigned int vofLS_c_index = vofLS_var.first_component_index;
    const unsigned int n_vofLS_dofs = vofLS_var.fe->dofs_per_cell;
    const unsigned int vofLS_blockidx = vofLS_var.block_index;

    //
    bool use_vof_composition = data->c_field_name!="";
    const unsigned int c_var_index =
      ( (!use_vof_composition) ? numbers::invalid_unsigned_int
        : sim.introspection.compositional_index_for_name(data->c_field_name));
    Simulator<dim>::AdvectionField advf =
      ( (!use_vof_composition) ? Simulator<dim>::AdvectionField::temperature()
        : Simulator<dim>::AdvectionField::composition(c_var_index));
    const unsigned int base_element =
      ( (!use_vof_composition) ? numbers::invalid_unsigned_int
        : advf.base_element(sim.introspection));
    const std::vector<Point<dim> > support_points =
      ( (!use_vof_composition) ? std::vector<Point<dim>>(0)
        : sim.finite_element.base_element(base_element).get_unit_support_points());

    //Iterate over cells
    for (auto cell : sim.dof_handler.active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        double cell_vof;

        // Obtain data for this cell and neighbors
        cell->get_dof_indices (local_dof_indicies);
        cell_vof = solution(local_dof_indicies[vof_ind]);

        normal[0] = 0.0;
        normal[1] = 0.0;
        d = -1.0;

        if (cell_vof > 1.0 - vof_epsilon)
          {
            d = 1.0;
          }
        else if (cell_vof > vof_epsilon)
          {
            //Identify best normal
            // Get references to neighboring cells
            for (unsigned int i = 0; i < 3; ++i)
              {
                typename DoFHandler<dim>::active_cell_iterator cen;
                if (i == 0)
                  {
                    if (cell->neighbor(0)==endc||cell->neighbor_is_coarser(0)||cell->neighbor(0)->has_children())
                      cen = endc;
                    else
                      cen = cell->neighbor(0);
                  }
                if (i == 1)
                  cen = cell;
                if (i == 2)
                  {
                    if (cell->neighbor(1)==endc||cell->neighbor_is_coarser(1)||cell->neighbor(1)->has_children())
                      cen = endc;
                    else
                      cen = cell->neighbor(1);
                  }
                for (unsigned int j = 0; j < 3; ++j)
                  {
                    typename DoFHandler<dim>::active_cell_iterator curr;
                    if (cen == endc)
                      {
                        curr = endc;
                      }
                    else
                      {
                        if (j == 0)
                          {
                            if (cen->neighbor(2)==endc||cen->neighbor_is_coarser(2)||cen->neighbor(2)->has_children())
                              curr = endc;
                            else
                              curr = cen->neighbor(2);
                          }
                        if (j == 1)
                          curr = cen;
                        if (j == 2)
                          {
                            if (cen->neighbor(3)==endc||cen->neighbor_is_coarser(3)||cen->neighbor(3)->has_children())
                              curr = endc;
                            else
                              curr = cen->neighbor(3);
                          }
                      }
                    if (curr != endc && !curr->has_children())
                      {
                        curr->get_dof_indices (cell_dof_indicies);
                        resc_cell_centers[3 * j + i] = Point<dim> (-1.0 + i,
                                                                   -1.0 + j);
                      }
                    else
                      {
                        cell->get_dof_indices (cell_dof_indicies);
                        resc_cell_centers[3 * j + i] = Point<dim> (0.0,
                                                                   0.0);
                      }
                    local_vofs (3 * j + i) = sim.solution (cell_dof_indicies[vof_ind]);
                  }
              }
            // Gather cell strip sums
            for (unsigned int i = 0; i < dim * n_sums; ++i)
              strip_sums[i] = 0.0;

            for (unsigned int i = 0; i < 3; ++i)
              {
                for (unsigned int j = 0; j < 3; ++j)
                  {
                    strip_sums[3 * 0 + i] += local_vofs (3 * j + i);
                    strip_sums[3 * 1 + j] += local_vofs (3 * j + i);
                  }
              }

            // std::cout << "  Strip sums: " << std::endl;
            // for (unsigned int i = 0; i < dim * n_sums; ++i)
            //   std::cout << "    " << strip_sums[i] << std::endl;

            // Calculate normal vectors for the 6 candidates from the efficient
            // least squares approach
            for (unsigned int di = 0; di < dim; ++di)
              {
                unsigned int di2 = (di + 1) % dim;
                for (unsigned int i = 0; i < 3; ++i)
                  {
                    normals[3 * di + i][di] = 0.0;
                    normals[3 * di + i][di2] = 0.0;
                    if (i % 2 == 0)
                      {
                        //use low sum
                        normals[3 * di + i][di] += strip_sums[3 * di + 0];
                        normals[3 * di + i][di2] += 1.0;
                      }
                    else
                      {
                        //use high sum
                        normals[3 * di + i][di] += strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    if (i == 0)
                      {
                        //use low sum
                        normals[3 * di + i][di] -= strip_sums[3 * di + 1];
                        normals[3 * di + i][di2] += 0.0;
                      }
                    else
                      {
                        //use high sum
                        normals[3 * di + i][di] -= strip_sums[3 * di + 2];
                        normals[3 * di + i][di2] += 1.0;
                      }
                    if (strip_sums[3 * di2 + 2] > strip_sums[3 * di2 + 0])
                      normals[3 * di + i][di2] *= -1.0;
                  }
              }

            // Add time extrapolated local normal as candidate
            // this is not expected to be the best candidate in general, but
            // should result in exact reconstruction for linear interface
            // translations
            // Inclusion of this candidate will not reduce accuracy due to it
            // only being selected if it produces a better interface
            // approximation than the ELS candidates. Note that this will
            // render linear translation problems less dependent on the
            // interface reconstruction, so other tests will also be necessary.
            for (unsigned int i=0; i<dim; ++i)
              normals[6][i] = sim.solution(local_dof_indicies[sim.finite_element
                                                              .component_to_system_index(vofN_c_index+i, 0)]);

            // If candidate normal too small, remove from consideration
            if (normals[6]*normals[6]< vof_epsilon)
              {
                normals[6][0] = 0;
                normals[6][1] = 0;
              }

            unsigned int mn_ind = 0;
            for (unsigned int nind = 0; nind < n_normals; ++nind)
              {
                errs[nind] = 0.0;
                d = VolumeOfFluid::d_from_vof<dim> (normals[nind], cell_vof);
                for (unsigned int i = 0; i < n_local; ++i)
                  {
                    double dot = 0.0;
                    for (unsigned int di = 0; di < dim; ++di)
                      dot += normals[nind][di] * resc_cell_centers[i][di];
                    double val = local_vofs (i)
                                 - VolumeOfFluid::vof_from_d<dim> (normals[nind],
                                                                   d - dot);
                    errs[nind] += val * val;
                  }
                if (errs[mn_ind] >= errs[nind])
                  mn_ind = nind;
                // std::cout << "   " << normals[nind] << " e ";
                // std::cout  << errs[nind] << " " << mn_ind << std::endl;
              }

            normal = normals[mn_ind];
            d = VolumeOfFluid::d_from_vof<dim> (normal, cell_vof);
          }

        double n2 = (normal*normal);
        if (n2 > vof_epsilon)
          {
            normal = (normal / n2);
            d = VolumeOfFluid::d_from_vof<dim> (normal, cell_vof);
          }
        else
          {
            normal[0] = 0.0;
            normal[1] = 0.0;
          }

        for (unsigned int i=0; i<dim; ++i)
          initial_solution (local_dof_indicies[sim.finite_element
                                               .component_to_system_index(vofN_c_index+i, 0)]) = normal[i];

        initial_solution (local_dof_indicies[sim.finite_element
                                             .component_to_system_index(vofN_c_index+dim, 0)]) = d;

        for (unsigned int i=0; i<n_vofLS_dofs; ++i)
          {
            // Recenter unit cell on origin
            Tensor<1, dim, double> uSupp = vofLS_var.fe->unit_support_point(i)-uReCen;
            initial_solution (local_dof_indicies[sim.finite_element
                                                 .component_to_system_index(vofLS_c_index, i)])
              = d-uSupp*normal;
          }

        // If specified, write unit cell linear approximation to specified composition variable
        if (use_vof_composition)
          {
            Tensor<1, dim, double> nnormal;
            double nd=d;
            double normall1n = 0.0;
            for (unsigned int i=0; i<dim; ++i)
              {
                normall1n += numbers::NumberTraits<double>::abs(normal[i]);
                nnormal[i] = 0.0;
              }
            if (normall1n > vof_epsilon)
              {
                nnormal = normal / normall1n;
                nd = d / normall1n;
              }
            //Calculate correct factor to retain vol frac and [0,1] bound
            double fact = 0.5-abs(cell_vof-0.5);
            for (unsigned int i=0; i<sim.finite_element.base_element(base_element).dofs_per_cell; ++i)
              {
                const unsigned int system_local_dof
                  = sim.finite_element.component_to_system_index(advf.component_index(sim.introspection),
                                                                 /*dof index within component*/i);

                Tensor<1, dim, double> uSupp = support_points[i]-uReCen;

                const double value = cell_vof - fact*(uSupp*nnormal);

                initial_solution(local_dof_indicies[system_local_dof]) = value;
              }
          }
      }

    initial_solution.compress(VectorOperation::insert);

    sim.compute_current_constraints();
    sim.current_constraints.distribute(initial_solution);

    sim.solution.block(vofN_blockidx) = initial_solution.block(vofN_blockidx);
    sim.solution.block(vofLS_blockidx) = initial_solution.block(vofLS_blockidx);
    if (use_vof_composition)
      {
        const unsigned int blockidx = advf.block_index(sim.introspection);
        sim.solution.block(blockidx) = initial_solution.block(blockidx);
      }

    sim.computing_timer.exit_section();
  }


  template <>
  void VoFHandler<3>::update_vof_normals (const VoFField<3> field,
                                          LinearAlgebra::BlockVector &solution)
  {
    Assert(false, ExcNotImplemented());
  }
}
