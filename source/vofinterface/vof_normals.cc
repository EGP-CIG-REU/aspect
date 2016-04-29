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
#include <aspect/vofinterface/vof_utils.h>

namespace aspect
{
  using namespace dealii;

  template <>
  void Simulator<2>::update_vof_normals (LinearAlgebra::BlockVector &solution)
  {
    const int dim = 2;

    LinearAlgebra::BlockVector initial_solution;

    computing_timer.enter_section("   Reconstruct VoF interfaces");

    initial_solution.reinit(system_rhs, false);

    // Boundary reference
    typename DoFHandler<dim>::active_cell_iterator endc =
      dof_handler.end ();

    // Interface Reconstruction vars
    const unsigned int n_local = 9;

    Vector<double> local_vofs (n_local);
    std::vector<Point<dim>> resc_cell_centers (n_local);

    const unsigned int n_sums = 3;
    std::vector<double> strip_sums (dim * n_sums);

    const unsigned int n_normals = 6;
    std::vector<Tensor<1, dim, double>> normals (n_normals);
    std::vector<double> errs (n_normals);

    // Normal holding vars
    Point<dim> uReCen;
    Tensor<1, dim, double> normal;
    double d;

    for (unsigned int i=0; i<dim; ++i)
      uReCen[i] = 0.5;

    std::vector<types::global_dof_index> cell_dof_indicies (finite_element.dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indicies (finite_element.dofs_per_cell);

    const FEVariable<dim> &vof_var = introspection.variable("vofs");
    const unsigned int vof_c_index = vof_var.first_component_index;
    const unsigned int vof_ind
      = finite_element.component_to_system_index(vof_c_index, 0);

    const FEVariable<dim> &vofN_var = introspection.variable("vofsN");
    const unsigned int vofN_c_index = vofN_var.first_component_index;
    const unsigned int vofN_blockidx = vofN_var.block_index;

    const FEVariable<dim> &vofLS_var = introspection.variable("vofsLS");
    const unsigned int vofLS_c_index = vofLS_var.first_component_index;
    const unsigned int n_vofLS_dofs = vofLS_var.fe->dofs_per_cell;
    const unsigned int vofLS_blockidx = vofLS_var.block_index;


    //Iterate over cells
    for (auto cell : dof_handler.active_cell_iterators ())
      {
        double cell_vof;

        // Obtain data for this cell and neighbors
        cell->get_dof_indices (local_dof_indicies);
        cell_vof = solution(local_dof_indicies[vof_ind]);

        normal[0] = 0.0;
        normal[1] = 0.0;
        d = -0.5;

        if (cell_vof > 1.0 - parameters.voleps)
          {
            d = 0.5;
          }
        else if (cell_vof > parameters.voleps)
          {
            //Identify best normal
            // Get references to neighboring cells
            for (unsigned int i = 0; i < 3; ++i)
              {
                typename DoFHandler<dim>::active_cell_iterator cen;
                if (i == 0)
                  cen = cell->neighbor (0);
                if (i == 1)
                  cen = cell;
                if (i == 2)
                  cen = cell->neighbor (1);
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
                          curr = cen->neighbor (2);
                        if (j == 1)
                          curr = cen;
                        if (j == 2)
                          curr = cen->neighbor (3);
                      }
                    if (curr != endc)
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
                    local_vofs (3 * j + i) = solution (cell_dof_indicies[vof_ind]);
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

            unsigned int mn_ind = 0;
            for (unsigned int nind = 0; nind < n_normals; ++nind)
              {
                errs[nind] = 0.0;
                d = InterfaceTracker::d_from_vof<dim> (normals[nind], cell_vof);
                for (unsigned int i = 0; i < n_local; ++i)
                  {
                    double dot = 0.0;
                    for (unsigned int di = 0; di < dim; ++di)
                      dot += normals[nind][di] * resc_cell_centers[i][di];
                    double val = local_vofs (i)
                                 - InterfaceTracker::vof_from_d<dim> (normals[nind],
                                                                      d - dot);
                    errs[nind] += val * val;
                  }
                if (errs[mn_ind] >= errs[nind])
                  mn_ind = nind;
                // std::cout << "   " << normals[nind] << " e ";
                // std::cout  << errs[nind] << " " << mn_ind << std::endl;
              }

            normal = normals[mn_ind];
            d = InterfaceTracker::d_from_vof<dim> (normal, cell_vof);
          }

        double n2 = (normal*normal);
        if (n2 > parameters.voleps)
          {
            normal = (normal / n2);
            d = InterfaceTracker::d_from_vof<dim> (normal, cell_vof);
          }
        else
          {
            normal[0] = 0.0;
            normal[1] = 0.0;
          }

        for (unsigned int i=0; i<dim; ++i)
          initial_solution (local_dof_indicies[finite_element
                                               .component_to_system_index(vofN_c_index+i, 0)]) = normal[i];

        initial_solution (local_dof_indicies[finite_element
                                             .component_to_system_index(vofN_c_index+dim, 0)]) = d;

        for (unsigned int i=0; i<n_vofLS_dofs; ++i)
          {
            // Recenter unit cell on origin
            Tensor<1, dim, double> uSupp = vofLS_var.fe->unit_support_point(i)-uReCen;
            initial_solution (local_dof_indicies[finite_element
                                                 .component_to_system_index(vofLS_c_index, i)])
              = d-uSupp*normal;
          }
      }

    initial_solution.compress(VectorOperation::insert);

    compute_current_constraints();
    current_constraints.distribute(initial_solution);

    solution.block(vofN_blockidx) = initial_solution.block(vofN_blockidx);
    solution.block(vofLS_blockidx) = initial_solution.block(vofLS_blockidx);

    computing_timer.exit_section();
  }


  template <>
  void Simulator<3>::update_vof_normals (LinearAlgebra::BlockVector &solution)
  {
    Assert(false, ExcNotImplemented());
  }
}
