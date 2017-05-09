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


#include <aspect/mesh_refinement/vof_boundary.h>
#include <aspect/vof/handler.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/fe/fe_values.h>

namespace aspect
{
  namespace MeshRefinement
  {
    template <int dim>
    void
    VoFBoundary<dim>::execute(Vector<float> &indicators) const
    {
      indicators = 0.0;

      const QMidpoint<dim> qMidC;
      const QMidpoint<dim-1> qMidF;

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               qMidC,
                               update_values |
                               update_quadrature_points);

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        qMidF,
                                        update_values |
                                        update_quadrature_points);

      FESubfaceValues<dim> fe_subface_values (this->get_mapping(),
                                              this->get_fe(),
                                              qMidF,
                                              update_values |
                                              update_quadrature_points);
      for (unsigned int f=0; f<this->get_vof_handler().get_n_fields(); ++f)
        {

          const FEValuesExtractors::Scalar vof_field = this->get_vof_handler().get_field(f).fraction.extractor_scalar();
          std::vector<double> vof_q_values(qMidC.size());

          // Should be vof_epsilon, look into how to accesss that
          double voleps = vof_epsilon;

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();
          unsigned int i=0;
          for (; cell!=endc; ++cell, ++i)
            {
              // Skip if not local
              if (!cell->is_locally_owned())
                continue;

              // Get cell vof
              double cell_vof;
              fe_values.reinit(cell);
              fe_values[vof_field].get_function_values (this->get_solution(),
                                                        vof_q_values);
              cell_vof = vof_q_values[0];

              // Handle overshoots
              if (cell_vof > 1.0)
                cell_vof = 1.0;

              if (cell_vof < 0.0)
                cell_vof = 0.0;

              // Check if at interface
              bool at_interface=false;
              if (cell_vof>voleps && cell_vof<(1.0-voleps))
                {
                  // Fractional volume
                  at_interface=true;
                }
              else
                {
                  // Check neighbors
                  std::vector<typename parallel::distributed::Triangulation<dim>::active_cell_iterator> neighbors;
                  GridTools::get_active_neighbors<parallel::distributed::Triangulation<dim> >(cell,neighbors);
                  for (unsigned int i=0; i<neighbors.size(); ++i)
                    {
                      fe_values.reinit(neighbors[i]);
                      fe_values.get_function_values (this->get_solution(),
                                                     vof_q_values);
                      double neighbor_vof = vof_q_values[0];
                      if (neighbor_vof>voleps && neighbor_vof<(1.0-voleps))
                        {
                          // Fractional volume
                          at_interface=true;
                        }
                      if (abs(cell_vof-neighbor_vof)>=voleps)
                        {
                          at_interface=true
                        }
                    }
                }

              if (at_interface)
                {
                  indicators(i) = 1.0;
                }
            }
        }
    }

    template <int dim>
    void
    VoFBoundary<dim>::tag_additional_cells() const
    {
      // Skip if do not have any vof data to use
      if (this->get_dof_handler().n_dofs()==0)
        return;

      // Currently do not need to do strong enforcement of refinement, will
      // consider at later point
    }

    template <int dim>
    void
    VoFBoundary<dim>::
    declare_parameters (ParameterHandler &prm)
    {
    }

    template <int dim>
    void
    VoFBoundary<dim>::parse_parameters (ParameterHandler &prm)
    {
      //TODO: Add check for vof active
      AssertThrow(this->get_parameters().vof_tracking_enabled,
                  ExcMessage("The 'vof boundary' mesh refinement strategy requires the 'Use VoF tracking' parameter be enabled."));

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
    ASPECT_REGISTER_MESH_REFINEMENT_CRITERION(VoFBoundary,
                                              "vof boundary",
                                              "A class that implements a mesh refinement criterion, which "
                                              "ensures a minimum level of refinement near the VoF interface boundary.")
  }
}
