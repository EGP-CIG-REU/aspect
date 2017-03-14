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

#include <aspect/particle/interpolator/bilinear_least_squares.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/lac/full_matrix.templates.h>

namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      template <int dim>
      std::vector<std::vector<double> >
      BilinearLeastSquares<dim>::properties_at_points(const std::multimap<types::LevelInd, Particle<dim> > &particles,
                                                      const std::vector<Point<dim> > &positions,
                                                      const ComponentMask &selected_properties,
                                                      const typename parallel::distributed::Triangulation<dim>::active_cell_iterator &cell) const
      {
        AssertThrow(dim == 2,
                    ExcMessage("Currently, the particle interpolator 'bilinear' is only supported for 2D models."));

        typename parallel::distributed::Triangulation<dim>::active_cell_iterator found_cell;

        Point<dim> approximated_cell_midpoint;

        if (cell == typename parallel::distributed::Triangulation<dim>::active_cell_iterator())
          {
            // We can not simply use one of the points as input for find_active_cell_around_point
            // because for vertices of mesh cells we might end up getting ghost_cells as return value
            // instead of the local active cell. So make sure we are well in the inside of a cell.
            Assert(positions.size() > 0,
                   ExcMessage("The particle property interpolator was not given any "
                              "positions to evaluate the particle cell_properties at."));

            approximated_cell_midpoint = std::accumulate (positions.begin(), positions.end(), Point<dim>())
                                         / static_cast<double> (positions.size());

            found_cell =
              (GridTools::find_active_cell_around_point<> (this->get_mapping(),
                                                           this->get_triangulation(),
                                                           approximated_cell_midpoint)).first;
          }
        else
          found_cell = cell;

        const types::LevelInd cell_index = std::make_pair<unsigned int, unsigned int> (found_cell->level(),found_cell->index());
        const std::pair<typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator,
              typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator> particle_range = particles.equal_range(cell_index);

        const unsigned int n_particles = std::distance(particle_range.first,particle_range.second);
        const unsigned int n_particle_properties = particles.begin()->second.get_properties().size();

        std::vector<std::vector<double> > cell_properties(positions.size(),
                                                          std::vector<double>(n_particle_properties, numbers::signaling_nan<double>()));

        unsigned int property_index = 0;
        for (unsigned int i=0; i < n_particle_properties; i++)
          if (selected_properties[i])
            property_index = i;

        AssertThrow(n_particles != 0,
                    ExcMessage("At least one cell contained no particles. The 'bilinear'"
                               "interpolation scheme does not support this case. "));


        const unsigned int matrix_dimension = 4;
        dealii::LAPACKFullMatrix<double> A(n_particles, matrix_dimension);
        Vector<double> r(n_particles);
        A = 0;
        r = 0;

        double max_value_for_particle_property = (particle_range.first)->second.get_properties()[property_index];
        double min_value_for_particle_property = (particle_range.first)->second.get_properties()[property_index];

        unsigned int index = 0;
        double cell_diameter = found_cell->diameter();
        for (typename std::multimap<types::LevelInd, Particle<dim> >::const_iterator particle = particle_range.first;
             particle != particle_range.second; ++particle, ++index)
          {
            const double particle_property_value = particle->second.get_properties()[property_index];
            r[index] = particle_property_value;

            const Point<dim> position = particle->second.get_location();
            A(index,0) = 1;
            A(index,1) = (position[0] - approximated_cell_midpoint[0])/cell_diameter;
            A(index,2) = (position[1] - approximated_cell_midpoint[1])/cell_diameter;
            A(index,3) = (position[0] - approximated_cell_midpoint[0]) * (position[1] - approximated_cell_midpoint[1])/std::pow(cell_diameter,2);

            if (max_value_for_particle_property < particle_property_value)
              max_value_for_particle_property = particle_property_value;
            if (min_value_for_particle_property > particle_property_value)
              min_value_for_particle_property = particle_property_value;
          }

        dealii::LAPACKFullMatrix<double> B(matrix_dimension, matrix_dimension);
        A.Tmmult(B, A, false);
        dealii::LAPACKFullMatrix<double> B_inverse(B);

        Vector<double> c_ATr(matrix_dimension);
        Vector<double> c(matrix_dimension);
        A.Tvmult(c_ATr,r);

        B_inverse.compute_svd();
        const double threshold = 1e-15;
        B_inverse.compute_inverse_svd(threshold);
        B_inverse.vmult(c, c_ATr);
        unsigned int index_positions = 0;
        for (typename std::vector<Point<dim> >::const_iterator itr = positions.begin(); itr != positions.end(); ++itr, ++index_positions)
          {
            Point<dim> support_point = *itr;
            double interpolated_value = c[0] +
                                        c[1]*(support_point[0] - approximated_cell_midpoint[0])/cell_diameter +
                                        c[2]*(support_point[1] - approximated_cell_midpoint[1])/cell_diameter +
                                        c[3]*(support_point[0] - approximated_cell_midpoint[0])*(support_point[1] - approximated_cell_midpoint[1])/std::pow(cell_diameter,2);

            if (interpolated_value > max_value_for_particle_property)
              interpolated_value = max_value_for_particle_property;
            else if (interpolated_value < min_value_for_particle_property)
              interpolated_value = min_value_for_particle_property;

            cell_properties[index_positions][property_index] = interpolated_value;
          }
        return cell_properties;
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Interpolator
    {
      ASPECT_REGISTER_PARTICLE_INTERPOLATOR(BilinearLeastSquares,
                                            "bilinear",
                                            "Interpolates particle properties onto a vector of points using a "
                                            "bilinear least squares method. Currently only 2D models are "
                                            "supported. ")
    }
  }
}
