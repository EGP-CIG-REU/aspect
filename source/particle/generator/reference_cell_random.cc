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
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#include <aspect/particle/generator/reference_cell_random.h>

#include <aspect/utilities.h>
#include <random>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <int dim>
      void
      ReferenceCellRandom<dim>::generate_particles(std::multimap<types::LevelInd, Particle<dim> > &particles)
      {
        types::particle_index n_particles_to_generate = this->get_triangulation().n_locally_owned_active_cells() * number_of_particles_per_cell;
        types::particle_index prefix_sum = 0;

        MPI_Scan(&n_particles_to_generate, &prefix_sum, 1, ASPECT_PARTICLE_INDEX_MPI_TYPE, MPI_SUM, this->get_mpi_communicator());

        types::particle_index particle_index = prefix_sum - n_particles_to_generate;

        typename Triangulation<dim>::active_cell_iterator
        cell = this->get_triangulation().begin_active(), endc = this->get_triangulation().end();

        int seed_for_random_number_generator = 0;

        for (; cell != endc; cell++)
          {
            if (cell->is_locally_owned())
              {
                seed_for_random_number_generator += 1;
                const std::vector<Point<dim> > particles_in_unit_cell = generate_random_particle_positions_in_unit_cell(seed_for_random_number_generator);

                for (typename std::vector<Point<dim> >::const_iterator itr_particles_in_unit_cell = particles_in_unit_cell.begin();
                     itr_particles_in_unit_cell != particles_in_unit_cell.end();
                     itr_particles_in_unit_cell++)
                  {
                    const Point<dim> position_real = this->get_mapping().transform_unit_to_real_cell(cell,
                                                     *itr_particles_in_unit_cell);
                    const Particle<dim> particle(position_real, *itr_particles_in_unit_cell, particle_index);
                    const types::LevelInd cell_index(cell->level(), cell->index());
                    particles.insert(std::make_pair(cell_index, particle));
                    ++particle_index;
                  }
              }
          }
      }


      template <int dim>
      std::vector<Point<dim> >
      ReferenceCellRandom<dim>::generate_random_particle_positions_in_unit_cell(int seed_for_random_number_generator)
      {
        std::vector<Point<dim> > particle_positions;
        //Initialize the random number generator seed
        //std::random_device rd;
        std::mt19937 generator(seed_for_random_number_generator);
        std::uniform_real_distribution<> dist(0,1);

        for (int i = 0; i < number_of_particles_per_cell; i ++)
          {
            if (dim == 2)
              {
                double x_coordinate = dist(generator);
                double y_coordinate = dist(generator);

                const Point<dim> position_unit = Point<dim>(x_coordinate,
                                                            y_coordinate);
                particle_positions.push_back(position_unit);
              }
            else if (dim == 3)
              {
                double x_coordinate = dist(generator);
                double y_coordinate = dist(generator);
                double z_coordinate = dist(generator);

                const Point<dim> position_unit = Point<dim>(x_coordinate,
                                                            y_coordinate,
                                                            z_coordinate);
                particle_positions.push_back(position_unit);
              }
            else
              ExcNotImplemented();
          }

        return particle_positions;
      }

      template <int dim>
      void
      ReferenceCellRandom<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell random");
              {
                prm.declare_entry ("Number of particles per cell", "16",
                                   Patterns::Integer(),
                                   "List of number of particles to randomly generate in each cell. ");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      ReferenceCellRandom<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("Reference cell random");
              {
                number_of_particles_per_cell = prm.get_integer("Number of particles per cell");
              }
              prm.leave_subsection();
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
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(ReferenceCellRandom,
                                         "reference cell random",
                                         "Generate a uniform distribution of particles per cell and spatial direction in "
                                         "the unit cell and transforms each of the particles back to real region in the model "
                                         "domain. Uniform here means the particles will be generated with an equal spacing in "
                                         "each spatial dimension")
    }
  }
}
