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


#include <aspect/postprocess/richardson_extrapolation.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>

#include <aspect/global.h>
#include <math.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/grid/grid_tools.h>

#include <aspect/simulator.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace Postprocess
  {

      template <int dim>
      void
      RichardsonExtrapolation<dim>::initialize() {
          output_file_name = this->get_output_directory()
                             + output_file_name + "_"
                             + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                             + ".dat";
          input_file_name = this->get_output_directory()
                            + input_file_name + "_"
                            + Utilities::int_to_string(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()))
                            + ".dat";
      }

      template <int dim>
      void
      RichardsonExtrapolation<dim>::read_in_data()
      {
          /**
         * Initializing private variables.
         */
          quadrature_points_input = new std::vector<Point<dim>>;
          temperature_input = new std::vector<double>;
          pressure_input = new std::vector<double>;
          velocity_input = new std::vector<Tensor<1,dim>>;
          compositional_fields_input = new std::vector<std::vector<double>>;

          for (unsigned int compositional_field_index=0; compositional_field_index < this->n_compositional_fields(); compositional_field_index++)
              compositional_fields_input->push_back(*(new std::vector<double>));

          /**
           * Looping and reading over data.
           */
          std::ifstream input(input_file_name);
          std::string line;

          while (std::getline(input, line))
          {
              Point<dim> q_point;
              Tensor<1,dim> velocity;
              double temperature, pressure;
              std::istringstream iss(line);

              iss >> q_point[0] >> q_point[1] >> velocity[0] >> velocity[1] >> pressure >> temperature;

              std::vector<std::vector<double>>::iterator itr_compositional_fields = compositional_fields_input->begin();
              for (; itr_compositional_fields != compositional_fields_input->end(); itr_compositional_fields++)
              {
                  double composition_value = 0;
                  iss >> composition_value;
                  itr_compositional_fields->push_back(composition_value);
              }

              quadrature_points_input->push_back(q_point);
              velocity_input->push_back(velocity);
              pressure_input->push_back(pressure);
              temperature_input->push_back(temperature);
          }

          input.close();

          /**
         * For debugging purposes, we now write out what;s been read in.
         */


          /**
          std::ofstream debugging;
          std::cout << input_file_name + ".g";
          debugging.open(input_file_name + ".g", std::ios_base::out);

          typename std::vector<double>::const_iterator itr_temperature = temperature_input->begin();
          typename std::vector<double>::const_iterator itr_pressure = pressure_input->begin();
          typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = velocity_input->begin();


          typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points_input->begin();

          for (int input_index = 0; itr_quadrature_points != quadrature_points_input->end(); input_index++, itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++)
          {
              debugging << (*itr_quadrature_points)[0] << "\t" << (*itr_quadrature_points)[1]
                        << "\t" << (*itr_velocity)[0] << "\t" << (*itr_velocity)[1]
                        << "\t" << *itr_pressure << "\t" << *itr_temperature;
              typename std::vector<std::vector<double>>::const_iterator itr_compositional_fields_input = compositional_fields_input->begin();
              for (; itr_compositional_fields_input != compositional_fields_input->end(); itr_compositional_fields_input++)
                debugging << "\t" << (*itr_compositional_fields_input)[input_index];
              debugging << std::endl;
          }
          debugging.close();
          **/


      }

      template <int dim>
      void
      RichardsonExtrapolation<dim>::write_out_data()
      {
          std::ofstream interpolated_data_stream;
          interpolated_data_stream.open(output_file_name, std::ios_base::out);
          interpolated_data_stream << std::setprecision(14);

          /**
          * Compute the Legendre gauss points at level 2 indirection.
          **/
          //QGaussLobatto<1> base_quadrature(this->get_stokes_velocity_degree() + 1);
          QGauss<1> base_quadrature(this->get_stokes_velocity_degree() + 1);
          QIterated<dim> quadrature_rule (base_quadrature, 2);

          FEValues<dim> fe_values(this->get_mapping(),
                                  this->get_fe(),
                                  quadrature_rule,
                                  update_values |
                                  update_quadrature_points |
                                  update_JxW_values);

          const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
          const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
          const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;

          // Declaring an iterator over all active cells on local mpi process
          typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
          for (; cell != this->get_dof_handler().end();
                 ++cell)
          {
              if (!(cell->is_locally_owned()))
                  continue;

              fe_values.reinit(cell);
              // Each vector is of the same length.
              const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();

              std::vector<double> interpolated_temperature(quadrature_points.size());
              std::vector<double> interpolated_pressure(quadrature_points.size());
              std::vector<Tensor<1,dim>> interpolated_velocity(quadrature_points.size());
              std::vector<std::vector<double>> interpolated_compositional_fields(this->n_compositional_fields(),
                                                                                 std::vector<double>(quadrature_points.size()));

              fe_values[extractor_pressure].get_function_values(this->get_solution(), interpolated_pressure);
              fe_values[extractor_temperature].get_function_values(this->get_solution(), interpolated_temperature);
              fe_values[extractor_velocity].get_function_values(this->get_solution(), interpolated_velocity);

              typename std::vector<std::vector<double>>::iterator itr_compositional_fields = interpolated_compositional_fields.begin();

              unsigned int index = 0;
              for(; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++) {
                      fe_values[this->introspection().extractors.compositional_fields[index]].get_function_values(
                              this->get_solution(),
                              *itr_compositional_fields);
                        index++;
                  }

              typename std::vector<double>::const_iterator itr_temperature = interpolated_temperature.begin();
              typename std::vector<double>::const_iterator itr_pressure = interpolated_pressure.begin();
              typename std::vector<Tensor<1,dim>>::const_iterator itr_velocity = interpolated_velocity.begin();

              typename std::vector<Point<dim>>::const_iterator itr_quadrature_points = quadrature_points.begin();

              unsigned int quadrature_point_index = 0;
              for (; itr_quadrature_points != quadrature_points.end();
                     itr_quadrature_points++, itr_velocity++, itr_pressure++, itr_temperature++, quadrature_point_index++)
              {
                  interpolated_data_stream << (*itr_quadrature_points)[0] << "\t" << (*itr_quadrature_points)[1]
                                           << "\t" << (*itr_velocity)[0] << "\t" << (*itr_velocity)[1]
                                           << "\t" << *itr_pressure << "\t" << *itr_temperature;
                  if (this->n_compositional_fields() != 0) {
                      unsigned int count = 0;
                      itr_compositional_fields = interpolated_compositional_fields.begin();
                      for(; itr_compositional_fields != interpolated_compositional_fields.end(); itr_compositional_fields++, count++)
                          interpolated_data_stream << "\t" << (*itr_compositional_fields)[quadrature_point_index] << "\t";
                  }
                  interpolated_data_stream << std::endl;
              }
          }

          interpolated_data_stream.close();
      }

    template <int dim>
    std::pair<std::string,std::string>
    RichardsonExtrapolation<dim>::execute (TableHandler &)
    {
      std::string info = "";

      if ( this->get_time() == end_time)
      {
        if (Utilities::fexists(input_file_name)) {
            /**
             * Initialize data structures with the solution data contained in ascii text file.
             */
            read_in_data();

            double velocity_l2_error = 0;
            double pressure_l2_error = 0;
            double temperature_l2_error = 0;
            double *compositional_field_l2_error = new double[this->n_compositional_fields()];
            for(unsigned int i = 0; i < this->n_compositional_fields(); i++)
                compositional_field_l2_error[i] = 0;

            const QGauss<dim> quadrature_formula(this->get_parameters().stokes_velocity_degree + 1);
//            const QGauss<dim> quadrature_formula_for_pressure(this->get_fe().base_element(this->introspection().base_elements.pressure).degree);

//            this->get_fe().base_element(this->introspection().base_elements.velocities).get_unit_support_points()
            FEValues<dim> fe_values(this->get_mapping(),
                                    this->get_fe(),
                                    quadrature_formula,
                                    update_values |
                                    update_quadrature_points |
                                    update_JxW_values);

          const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
          const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
            const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;

            // iterate over all active cells on local mpi process
            typename std::vector<double>::const_iterator itr_temperature_input = temperature_input->begin();
            typename std::vector<double>::const_iterator itr_pressure_input = pressure_input->begin();
            typename std::vector<Tensor<1, dim>>::const_iterator itr_velocity_input = velocity_input->begin();
            typename std::vector<std::vector<double>>::const_iterator itr_compositional_fields_input = compositional_fields_input->begin();

            typename std::vector<Point<dim>>::const_iterator itr_input_quadrature_points = quadrature_points_input->begin();

            unsigned int quadrature_point_index = 0;
            for (; itr_input_quadrature_points !=
                   quadrature_points_input->end(); quadrature_point_index++, itr_input_quadrature_points++, itr_velocity_input++, itr_pressure_input++, itr_temperature_input++) {
                  std::pair<const typename DoFHandler<dim>::active_cell_iterator,
                          Point<dim> > it = GridTools::find_active_cell_around_point(this->get_mapping(), this->get_dof_handler(), *itr_input_quadrature_points);
                  if (! it.first->is_locally_owned())
                      continue;

                fe_values.reinit(it.first);

                  std::vector<double> temperature(fe_values.n_quadrature_points);
                  std::vector<double> pressure(fe_values.n_quadrature_points);
                  std::vector<Tensor<1,dim>> velocity(fe_values.n_quadrature_points);
                  std::vector<std::vector<double>> compositional_fields (this->n_compositional_fields(), std::vector<double> (fe_values.n_quadrature_points));

                const std::vector<Point<dim>>  quadrature_points = fe_values.get_quadrature_points();
                  const std::vector<double>  jacobian_weights = fe_values.get_JxW_values();

                  fe_values[extractor_pressure].get_function_values(this->get_solution(), pressure);
                  fe_values[extractor_temperature].get_function_values(this->get_solution(), temperature);
                  fe_values[extractor_velocity].get_function_values(this->get_solution(), velocity);

                unsigned int index = 0;
                typename std::vector<std::vector<double>>::iterator itr_compositional_fields = compositional_fields.begin();

                for(; itr_compositional_fields != compositional_fields.end(); itr_compositional_fields++) {
                    fe_values[this->introspection().extractors.compositional_fields[index]].get_function_values(
                            this->get_solution(),
                            *itr_compositional_fields);
                    index++;
                }

                  index = 1000;

                const double tol_tmp = it.first->diameter() * 1e-2;
                  for (unsigned int i=0; i<fe_values.n_quadrature_points; i++)
                  {
                    Point<dim> tmp((*itr_input_quadrature_points) - quadrature_points[i]);
                      if ( std::sqrt(tmp.square()) <= tol_tmp) //&& *itr_weights== jacobian_weight_points[i])
                      {
                          index = i;
                          velocity_l2_error += (velocity[index] - (*itr_velocity_input))*(velocity[index] - (*itr_velocity_input))*(jacobian_weights[index]);
                          pressure_l2_error += (pressure[index] - (*itr_pressure_input))*(pressure[index] - (*itr_pressure_input))*(jacobian_weights[index]);
                          temperature_l2_error += (temperature[index] - (*itr_temperature_input))*(temperature[index] - (*itr_temperature_input))*(jacobian_weights[index]);

                          itr_compositional_fields = compositional_fields.begin();
                          itr_compositional_fields_input = compositional_fields_input->begin();
                          unsigned int compositional_field_index = 0;

                          for (; itr_compositional_fields != compositional_fields.end() && itr_compositional_fields_input != compositional_fields_input->end();
                                 itr_compositional_fields++, itr_compositional_fields_input++)
                          {
                              compositional_field_l2_error[compositional_field_index] += ((*itr_compositional_fields)[index] - (*itr_compositional_fields_input)[quadrature_point_index]) *
                                                                                         ((*itr_compositional_fields)[index] - (*itr_compositional_fields_input)[quadrature_point_index]) *
                                                                                         (jacobian_weights[index]);
                              compositional_field_index++;
                          }
                          break;
                      }
                      else if (i == (fe_values.n_quadrature_points-1))
                      {
                          std::cout << "Tolerance of: " << tol_tmp << std::endl;
                          std::cout << "x: "  << (*itr_input_quadrature_points)[0] << " y: " << (*itr_input_quadrature_points)[1] << std::endl;
                          for (unsigned int i=0; i<fe_values.n_quadrature_points; i++)
                          {
                              std::cout << "Index " << i << ": " << quadrature_points[i] << std::endl;
                          }
                          /**
                           * TODO: Make sure that AMR has not been set!!
                           */

                              Assert(false, ExcInternalError());
                      }
                  }
            }

          velocity_l2_error = std::sqrt(Utilities::MPI::sum(velocity_l2_error, this->get_mpi_communicator()));
          pressure_l2_error = std::sqrt(Utilities::MPI::sum(pressure_l2_error, this->get_mpi_communicator()));
          temperature_l2_error = std::sqrt(Utilities::MPI::sum(temperature_l2_error, this->get_mpi_communicator()));
          for(unsigned int compositional_field_index = 0; compositional_field_index < this->n_compositional_fields(); compositional_field_index++)
              compositional_field_l2_error[compositional_field_index] = std::sqrt(Utilities::MPI::sum(compositional_field_l2_error[compositional_field_index], this->get_mpi_communicator()));

          if(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0) {
              std::ofstream error_log(data_output_file_name, std::ios_base::app);
              info = std::to_string(velocity_l2_error);
              error_log << std::setprecision(14) << velocity_l2_error << " " << pressure_l2_error << " " << temperature_l2_error;
              for(unsigned int compositional_field_index = 0; compositional_field_index < this->n_compositional_fields(); compositional_field_index++)
                  error_log << " " << compositional_field_l2_error[compositional_field_index];
              error_log << std::endl;
              error_log.close();
          }
        }

        // Write out the current solution, interpolated at a higher resolved mesh.
        write_out_data();
      }

      return std::pair<std::string, std::string> ("Richardson Extrapolation: ",
                                                  info);
    }


    // Add function to compute the error between interpolated read in values with the current solution at current nodal points.

    template <int dim>
    void
    RichardsonExtrapolation<dim>::declare_parameters (ParameterHandler &prm)
    {
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Richardson extrapolation");
            {
                prm.declare_entry ("End time",
                        /* boost::lexical_cast<std::string>(std::numeric_limits<double>::max() /
                                             year_in_seconds) = */ "5.69e+300",
                                   Patterns::Double (),
                                   "The end time of the simulation. The default value is a number "
                                           "so that when converted from years to seconds it is approximately "
                                           "equal to the largest number representable in floating point "
                                           "arithmetic. For all practical purposes, this equals infinity. "
                                           "Units: Years if the "
                                           "'Use years in output instead of seconds' parameter is set; "
                                           "seconds otherwise.");
                prm.declare_entry("Input file name of interpolated data", "interpolated_data",
                                  Patterns::Anything (),
                                  "The file name containing the interpolated solution at the nodal values "
                                  "at the current grid resolution.");
                prm.declare_entry("Output file name of interpolated data", "interpolated_data",
                                  Patterns::Anything (),
                                  "The file name containing the interpolated solution at the nodal values "
                                  "at the current grid resolution * 4 cells.");
                prm.declare_entry("Output file name of errors", "Errors.dat",
                                  Patterns::Anything(),
                                  "A file name to append the computed little L2 error oo the interpolated "
                                  "solution to current solution at quadrature points.");
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
    }

    template <int dim>
    void
    RichardsonExtrapolation<dim>::parse_parameters (ParameterHandler &prm)
    {
        bool use_years = prm.get_bool ("Use years in output instead of seconds");
        prm.enter_subsection("Postprocess");
        {
            prm.enter_subsection("Richardson extrapolation");
            {
                // read end time from parameter file. if it is to be interpreted
                // in years rather than seconds, then do the conversion
                end_time = prm.get_double ("End time");
                if ( use_years == true)
                    end_time *= year_in_seconds;

                input_file_name = prm.get("Input file name of interpolated data");
                output_file_name = prm.get("Output file name of interpolated data");
                data_output_file_name = prm.get("Output file name of errors");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(RichardsonExtrapolation,
                                  "richardson extrapolation",
                                  "A postprocessor called Richardson Extrapolation.")
  }
}
