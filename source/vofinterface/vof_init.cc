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

// #include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/parsed_function.h>

namespace aspect
{

  template <int dim>
  void Simulator<dim>::set_initial_vofs ()
  {
    if (!parameters.vof_tracking_enabled)
      return;

    unsigned int n_i_samp = 3;

    LinearAlgebra::BlockVector initial_solution;

    initial_solution.reinit(system_rhs, false);

    // Main code

    Assert(dim==2, ExcMessage("VOF interface tracking not implemented for dim>2"));

    const QIterated<dim> quadrature (QMidpoint<1>(), n_i_samp);
    FEValues<dim, dim> fe_init (mapping, finite_element, quadrature,
                                update_JxW_values | update_quadrature_points);

    double h = 1.0/n_i_samp;

    std::vector<types::global_dof_index>
    local_dof_indicies (finite_element.dofs_per_cell);

    const FEVariable<dim> &vof_var = introspection.variable("vofs");
    const unsigned int base_element = vof_var.base_index;
    const unsigned int component_index = vof_var.first_component_index;
    const unsigned int blockidx = vof_var.block_index;
    const unsigned int vof_ind
      = finite_element.component_to_system_index(component_index, 0);
    std::cout << vof_ind << std::endl;

    std::cout << blockidx << std::endl;

    // if (sim.convert_output_to_years ())
    //   initFunc.set_time (sim.get_time ()/year_in_seconds);
    // else
    //   initFunc.set_time (sim.get_time());

    // Initialize state based on provided function
    for (auto cell : dof_handler.active_cell_iterators ())
      {
        if (!cell->is_locally_owned ())
          continue;

        // Calculate approximation for volume
        double cell_vol, cell_diam, d_func;
        cell->get_dof_indices (local_dof_indicies);


        cell_vol = cell->measure ();
        cell_diam = cell->diameter();
        //d_func = initFunc.value(cell->barycenter());
        fe_init.reinit (cell);

        double vof_val = 0.5;

        /*
         * Temporarily disable calculation
        if (d_func <=-0.5*cell_diam)
          {
            vof_val = 0.0;
          }
        else
          {
            if (d_func >= 0.5*cell_diam)
              {
                vof_val = 1.0;
              }
            else
              {

                for (unsigned int i = 0; i < fe_init.n_quadrature_points; ++i)
                  {
                    double d = 0.0;
                    Tensor<1, dim, double> grad;
                    Point<dim> xU = quadrature.point (i);
                    for (unsigned int di = 0; di < dim; ++di)
                      {
                        Point<dim> xH, xL;
                        xH = xU;
                        xL = xU;
                        xH[di] += 0.5*h;
                        xL[di] -= 0.5*h;
                        double dH = initFunc.value(cell->intermediate_point(xH));
                        double dL = initFunc.value(cell->intermediate_point(xL));
                        grad[di] = (dL-dH);
                        d += (0.5/dim)*(dH+dL);
                      }
                    double ptvof = vof_from_d<dim> (grad, d);
                    vof_val += ptvof * (fe_init.JxW (i) / cell_vol);
                  }
              }
          }
          */
        initial_solution (local_dof_indicies[vof_ind]) = vof_val;
      }

    initial_solution.compress(VectorOperation::insert);

    compute_current_constraints();
    current_constraints.distribute(initial_solution);

    solution.block(blockidx) = initial_solution.block(blockidx);
    old_solution.block(blockidx) = initial_solution.block(blockidx);
    old_old_solution.block(blockidx) = initial_solution.block(blockidx);
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::set_initial_vofs ();

  ASPECT_INSTANTIATE(INSTANTIATE)
}
