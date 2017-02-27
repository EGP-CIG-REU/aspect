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


#include <aspect/postprocess/visualization/gravity.h>



namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      Gravity<dim>::
      Gravity ()
        :
        DataPostprocessorVector<dim> ("gravity",
                                      update_q_points)
      {}



      template <int dim>
      void
      Gravity<dim>::
      compute_derived_quantities_vector (const std::vector<Vector<double> > &,
                                         const std::vector<std::vector<Tensor<1,dim> > > &,
                                         const std::vector<std::vector<Tensor<2,dim> > > &,
                                         const std::vector<Point<dim> > &,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const
      {
        const unsigned int n_quadrature_points = evaluation_points.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == dim,    ExcInternalError());

        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const Tensor<1,dim> g = this->get_gravity_model().gravity_vector (evaluation_points[q]);
            for (unsigned int k=0; k<dim; ++k)
              computed_quantities[q](k) = g[k];
          }
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(Gravity,
                                                  "gravity",
                                                  "A visualization output object that outputs the gravity vector.")
    }
  }
}
