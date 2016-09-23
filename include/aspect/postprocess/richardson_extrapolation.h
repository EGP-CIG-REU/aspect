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


#ifndef __aspect__postprocess_richardson_extrapolation_h
#define __aspect__postprocess_richardson_extrapolation_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {

    /**
     * A postprocessor that evaluates the solution vector at individual
     * points.
     *
     * @ingroup Postprocessing
     */
    template <int dim>
    class RichardsonExtrapolation : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        virtual
        void initialize ();

        /**
         * Evaluate the solution and determine the values at the
         * selected points.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * Read in the solution data. Currently, we read in an ascii input file.
         * Ideally, we read from a binary file.
         */
        void
        read_in_data();

        /**
         * Write out the solution data. Again, this is written out into an ascii
         * file and instead, should be written out in binary format.
         */
        void
        write_out_data();

        void
        write_out_current_data();

        /**
         *  Function that computes the error between interpolated (read in) values
         *  with the current solution at current nodal points.
         */
        void
        compute_error ();
      private:
        /**
         * Run time parameters.
         */
        std::string input_file_name;
        std::string output_file_name;
        double end_time;

        /**
         * Data structure to store read in interpolated solution.
         */
        std::vector<Point<dim>> *quadrature_points_input;
        std::vector<double> *temperature_input;
        std::vector<double> *pressure_input;
        std::vector<Tensor<1,dim>> *velocity_input;
        std::vector<std::vector<double>> *compositional_fields_input;
    };
  }
}


#endif
