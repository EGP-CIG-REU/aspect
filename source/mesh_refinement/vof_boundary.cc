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

      const FEValuesExtractors::Scalar vof_field = this->introspection().variable("vofs").extractor_scalar();
      std::vector<double> vof_q_values(qMidC.size());

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
          double voleps = this->get_parameters().voleps;
          if (cell_vof>voleps && cell_vof<(1.0-voleps))
            {
              // Fractional volume
              at_interface=true;
            }
          else
            {
              // Check neighbors
              for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
                {
                  if (cell->face(face_no)->at_boundary())
                    continue;

                  typename DoFHandler<dim>::face_iterator face = cell->face (face_no);

                  const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);

                  double face_vof, n_face_vof;

                  fe_face_values.reinit(cell, face_no);

                  fe_face_values[vof_field].get_function_values (this->get_solution(),
                                                                 vof_q_values);

                  face_vof = vof_q_values[0];

                  if (!(face->has_children()))
                    {
                      if (!cell->neighbor_is_coarser(face_no))
                        {

                          const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

                          fe_face_values.reinit(neighbor, neighbor2);

                          fe_face_values[vof_field].get_function_values (this->get_solution(),
                                                                         vof_q_values);

                          n_face_vof = vof_q_values[0];
                        }
                      else
                        {
                          std::pair<unsigned int, unsigned int> n2pair = cell->neighbor_of_coarser_neighbor(face_no);
                          fe_subface_values.reinit(neighbor, n2pair.first, n2pair.second);

                          fe_subface_values[vof_field].get_function_values (this->get_solution(),
                                                                            vof_q_values);

                          n_face_vof = vof_q_values[0];
                        }

                      if (n_face_vof > 1.0)
                        n_face_vof = 1.0;
                      if (n_face_vof < 0.0)
                        n_face_vof = 0.0;

                      if (abs(n_face_vof-cell_vof)>=voleps)
                        at_interface=true;
                    }
                  else
                    {
                      const unsigned int neighbor2 = cell->neighbor_face_no(face_no);

                      for (unsigned int subface_no=0; subface_no<face->number_of_children(); ++subface_no)
                        {
                          const typename DoFHandler<dim>::active_cell_iterator neighbor_child
                            = cell->neighbor_child_on_subface (face_no, subface_no);

                          fe_face_values.reinit (neighbor_child, neighbor2);

                          fe_face_values[vof_field].get_function_values (this->get_solution(),
                                                                         vof_q_values);

                          if (n_face_vof > 1.0)
                            n_face_vof = 1.0;
                          if (n_face_vof < 0.0)
                            n_face_vof = 0.0;

                          n_face_vof = vof_q_values[0];

                          if (abs(n_face_vof-cell_vof)>=voleps)
                            at_interface=true;
                        }
                    }
                }
            }

          if (at_interface)
            {
              indicators(i) = 1.0;
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
      prm.enter_subsection("Mesh refinement");
      {

        prm.enter_subsection("VoF boundary");
        {
          prm.declare_entry ("Boundary refinement", "0",
                             Patterns::Integer(0),
                             "Minimum refinement level near fluid interface");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    VoFBoundary<dim>::parse_parameters (ParameterHandler &prm)
    {
      //TODO: Add check for vof active
      prm.enter_subsection("Mesh refinement");
      {
        prm.enter_subsection("VoF boundary");
        {
          min_interface_level = prm.get_integer ("Boundary refinement");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
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
