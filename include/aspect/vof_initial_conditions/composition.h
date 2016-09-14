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


#ifndef __aspect__vof_initial_conditions_composition_h
#define __aspect__vof_initial_conditions_composition_h

#include <aspect/vof_initial_conditions/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace VoFInitialConditions
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a specified composition field.
     *
     * @ingroup VoFInitialConditionsModels
     */
    template <int dim>
    class Composition : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */

        Function ();

        /**
         * Return number of sample points to use for initialization.
         */
        virtual
        unsigned int n_samples () const;

        /**
         * Return whether initialization is a signed distance level set function.
         */
        virtual
        typename VoFInitType::Kind init_type() const;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        virtual
        double initial_value (const Point<dim> &position) const;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        /**
         * Required information on initialization method.
         */

        unsigned int n_init_samples;

        /**
         * Name of composition field for consideration
         */
        unsigned int c_idx;
    };
  }
}


#endif
