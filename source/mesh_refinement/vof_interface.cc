/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/mesh_refinement/vof_interface.h>
#include <aspect/vof/handler.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/geometry_info.h>

#include <deal.II/numerics/error_estimator.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

namespace aspect
{
  namespace MeshRefinement
  {

    template <int dim>
    void
    VoFInterface<dim>::tag_additional_cells() const
    {
        const QMidpoint<dim> qMidC;

        // Create a map from vertices to adjacent cells
        const std::vector<std::set<typename Triangulation<dim>::active_cell_iterator> >
                vertex_to_cells(GridTools::vertex_to_cell_map(this->get_triangulation()));

//        std::vector<std::set<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> > vert_cell_map =
//                vertex_to_cells(GridTools::vertex_to_cell_map(this->get_dof_handler().get_triangulation()));

        std::set<typename Triangulation<dim>::active_cell_iterator> marked_cells;
        FEValues<dim> fe_values (this->get_mapping(),
                                 this->get_fe(),
                                 qMidC,
                                 update_values |
                                 update_quadrature_points);
        for (unsigned int f=0; f<this->get_vof_handler().get_n_fields(); ++f) {

            const FEValuesExtractors::Scalar vof_field = this->get_vof_handler().get_field(f).fraction.extractor_scalar();
            std::vector<double> vof_q_values(qMidC.size());

            // Should be vof_epsilon, look into how to access that
            double voleps = vof_epsilon;

            typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active(),
                    endc = this->get_dof_handler().end();
            for (; cell != endc; ++cell) {
                // Skip if not local
                if (!cell->is_locally_owned())
                    continue;

                // Get cell vof
                double cell_vof;
                fe_values.reinit(cell);
                fe_values[vof_field].get_function_values(this->get_solution(),
                                                         vof_q_values);
                cell_vof = vof_q_values[0];

                // Handle overshoots
                if (cell_vof > 1.0)
                    cell_vof = 1.0;

                if (cell_vof < 0.0)
                    cell_vof = 0.0;

                // Check if at interface
                if (cell_vof > voleps && cell_vof < (1.0 - voleps))
                {
                    // Fractional volume
                    marked_cells.insert(cell);
                    cell->clear_coarsen_flag ();
                    cell->set_refine_flag ();
                }
            }
        }
        // Now mark for refinement all cells that are a neighbor of a cell that contains the interface

        typename std::set<typename parallel::distributed::Triangulation<dim>::active_cell_iterator>::const_iterator mcells = marked_cells.begin(),
                endmc = marked_cells.end();
        for (; mcells != endmc; mcells++)
        {
            typename parallel::distributed::Triangulation<dim>::active_cell_iterator mcell = *mcells;
            for(unsigned int vertex_index=0; vertex_index<GeometryInfo<dim>::vertices_per_cell;++vertex_index)
              {
                  std::set<typename Triangulation<dim>::active_cell_iterator> neighbor_cells = vertex_to_cells[mcell->vertex_index(vertex_index)];
                  typename std::set<typename Triangulation<dim>::active_cell_iterator>::const_iterator neighbor_cell = neighbor_cells.begin(),
                  end_neighbor_cell_index= neighbor_cells.end();
                  for(; neighbor_cell!=end_neighbor_cell_index;neighbor_cell++){
                      typename Triangulation<dim>::active_cell_iterator itr_tmp = *neighbor_cell;
                      itr_tmp->clear_coarsen_flag ();
                      itr_tmp->set_refine_flag ();
                  }
              }
        }

//        typename DoFHandler<dim>::active_cell_iterator
//                cell = this->get_dof_handler().begin_active(),
//                endc = this->get_dof_handler().end();
//        for (; cell != endc; ++cell) {
//            if (cell->is_locally_owned()){
//                if (marked_cells_and_neighbors.find(cell) != marked_cells_and_neighbors.end()){
//                    cell->clear_coarsen_flag();
//                    cell->set_refine_flag();
//                }
//            else
//            {
//                if (cell->refine_flag_set())
//                    cell->set_coarsen_flag();
//                }
//            }
//        }

//        std::set<typename Triangulation<dim>::active_cell_iterator> marked_cells_and_neighbors = marked_cells;
//        for (cell=marked_cells_and_neighbors);
//          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
//              marked_cells_and_neighbors.insert (vertex_to_cell_array[cell->vertex(v)])

    }

    template <int dim>
    void
    VoFInterface<dim>::
    declare_parameters (ParameterHandler &/*prm*/)
    {
    }

    template <int dim>
    void
    VoFInterface<dim>::parse_parameters (ParameterHandler &prm)
    {
      //TODO: Add check for vof active
      AssertThrow(this->get_parameters().vof_tracking_enabled,
                  ExcMessage("The 'vof boundary' mesh refinement strategy requires that the 'Use VoF tracking' parameter is enabled."));

      prm.enter_subsection ("VoF config");
      {
        vof_epsilon = prm.get_double("Small volume");
      }
      prm.leave_subsection ();
    }


  }
}

// explicit instantiations
namespace aspect
{
  namespace MeshRefinement
  {
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(VoFInterface,
                                              "vof interface",
                                              "A class that implements a mesh refinement criterion, which "
                                              "ensures a minimum level of refinement near the VoF interface boundary.")
  }
}
