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
        return std::vector<DataComponentInterpretation::DataComponentInterpretation>
               (1, DataComponentInterpretation::component_is_scalar);
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
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> > &,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = uh.size();
        Assert (computed_quantities.size() == n_quadrature_points, ExcInternalError ());
        Assert (uh[0].size() == this->introspection().n_components, ExcInternalError ());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            computed_quantities[q][0] =
              uh[q][this->introspection().variable("vofs").first_component_index];
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
              AssertThrow(vof_names.size() == this->introspection().variable("vofs").n_components(),
                          ExcMessage("You must define names for all VoF components"));
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
