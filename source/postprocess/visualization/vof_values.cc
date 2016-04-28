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

#include <aspect/postprocess/visualization/vof_values.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VoFValues<dim>::
      VoFValues ()
        :
        DataPostprocessor<dim> ()
      {}


      template <int dim>
      std::vector<std::string>
      VoFValues<dim>::
      get_names () const
      {
        return vof_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      VoFValues<dim>::
      get_data_component_interpretation () const
      {
        return interp;
      }


      template <int dim>
      UpdateFlags
      VoFValues<dim>::
      get_needed_update_flags () const
      {
        return update_values;
      }


      template <int dim>
      void
      VoFValues<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> > &,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError ());
        Assert (uh[0].size() == this->introspection().n_components, ExcInternalError ());

        const FiniteElement<dim> &finite_element = this->get_fe();

        const FEVariable<dim> &vof_var = this->introspection().variable("vofs");
        const unsigned int vof_ind = vof_var.first_component_index;
        const FEVariable<dim> &vofLS_var = this->introspection().variable("vofsLS");
        const unsigned int vofLS_ind = vofLS_var.first_component_index;

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            unsigned int out_ind = 0;
            computed_quantities[q][out_ind] = uh[q][vof_ind];
            ++out_ind;
            if (include_vofLS)
              {
                computed_quantities[q][out_ind] = uh[q][vofLS_ind];
                ++out_ind;
              }

            if (include_vofN)
              {
                Tensor<1, dim, double> normal = -duh[q][vofLS_ind];
                for (unsigned int i = 0; i<dim; ++i)
                  computed_quantities[q][out_ind] = normal[i];
              }
          }
      }


      template <int dim>
      void
      VoFValues<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("VoF values");
            {
              prm.declare_entry("Names of vofs", "VOF1",
                                Patterns::List (Patterns::Anything()),
                                "Names of vectors as they will appear in the output.");

              prm.declare_entry("Include internal reconstruction LS", "false",
                                Patterns::Bool (),
                                "Include the internal level set data use to save reconstructed interfaces");

              prm.declare_entry("Include normals", "false",
                                Patterns::Bool (),
                                "Include internal normal data in output (DEBUG)");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      VoFValues<dim>::parse_parameters (ParameterHandler &prm)
      {
        Assert (this->introspection ().variable_exists("vofs"),
                ExcMessage("VoF values not available"));

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("VoF values");
            {
              vof_names = Utilities::split_string_list(prm.get("Names of vofs"), ',');
              interp.push_back(DataComponentInterpretation::component_is_scalar);
              AssertThrow(vof_names.size() == 1,
                          ExcMessage("Only 1 VoF supported"));

              include_vofLS = prm.get_bool("Include internal reconstruction LS");
              if (include_vofLS)
                {
                  vof_names.push_back("vofsLS");
                  interp.push_back(DataComponentInterpretation::component_is_scalar);
                }

              include_vofN = prm.get_bool("Include normals");
              if (include_vofN)
                {
                  for (unsigned int i=0; i<dim; ++i)
                    {
                      vof_names.push_back("vofINormal");
                      interp.push_back(DataComponentInterpretation::component_is_part_of_vector);
                    }
                }
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VoFValues,
                                                  "vof values",
                                                  "A visualization output object that outputs the  vof data."
                                                  "Names are given in Postprocess/Visualization/VoF values")
    }
  }
}
