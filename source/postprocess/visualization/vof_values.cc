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
#include <aspect/vof/handler.h>


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

        unsigned int out_per_field=1;
        if (include_vofLS)
          out_per_field += 1;
        if (include_vofN)
          out_per_field += dim;

        for (unsigned int f=0; f<this->get_vof_handler().get_n_fields(); ++f)
          {
            VoFField<dim> field = this->get_vof_handler().get_field(f);

            const FEVariable<dim> &vof_var = field.fraction;
            const unsigned int vof_ind = vof_var.first_component_index;
            const FEVariable<dim> &vofLS_var = field.level_set;
            const unsigned int vofLS_ind = vofLS_var.first_component_index;

            for (unsigned int q=0; q<n_quadrature_points; ++q)
              {
                unsigned int out_ind = f*out_per_field;
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
                      {
                        computed_quantities[q][out_ind] = normal[i];
                        ++out_ind;
                      }
                  }
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
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("VoF values");
            {
              include_vofLS = prm.get_bool("Include internal reconstruction LS");
              include_vofN = prm.get_bool("Include normals");

              for (unsigned int f=0; f<this->get_vof_handler().get_n_fields(); ++f)
                {
                  std::string field_name = this->get_vof_handler().get_field_name(f);
                  vof_names.push_back("vof_"+field_name);
                  interp.push_back(DataComponentInterpretation::component_is_scalar);

                  if (include_vofLS)
                    {
                      vof_names.push_back("vofLS_"+field_name);
                      interp.push_back(DataComponentInterpretation::component_is_scalar);
                    }

                  if (include_vofN)
                    {
                      for (unsigned int i=0; i<dim; ++i)
                        {
                          vof_names.push_back("vofINormal_"+field_name);
                          interp.push_back(DataComponentInterpretation::component_is_part_of_vector);
                        }
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
