/*
  Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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


#ifndef __aspect__vof_initial_conditions_interface_h
#define __aspect__vof_initial_conditions_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  /**
   * A namespace in which we define everything that has to do with defining
   * the initial conditions.
   *
   * @ingroup VoFInitialConditionsModels
   */
  namespace VoFInitialConditions
  {
    using namespace dealii;

    /**
     * A structure that contains enum values that identify type of
     * initialization. Interpolation approach has less resolution than desired,
     * and the nature of the data allows a more accurate approximation.
     */
    struct VoFInitType
    {
      enum Kind
      {
        /**
         * Initialization data is a composition value between 0 and 1 at all
         * points and should be accumulated by integration.
         */
        composition,
        /**
         * Initialization data is an interface defined by a signed distance
         * level set with positive value indicating fluid presence.
         */
        signed_distance_level_set
      };
    };

    /**
     * A base class for parameterizations of volume-of-fluid initial
     * conditions.
     *
     * @ingroup VoFInitialConditionsModels
     */
    template <int dim>
    class Interface
    {
      public:
        /**
         * Destructor. Made virtual to enforce that derived classes also have
         * virtual destructors.
         */
        virtual ~Interface();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        virtual
        void
        initialize ();

        /**
         * Return number of sample points to use for initialization.
         * Initialization is done by iterated midpoint quadrature integration
         * with a geometrically justified smoothing chosen to handle linear
         * interfaces well for the signed distance function initialization
         * data.
         */
        virtual
        unsigned int n_samples () const = 0;

        /**
         * Return which type of initialization data is being used.
         */
        virtual
        typename VoFInitType::Kind init_type() const = 0;

        /**
         * Return the initial value as a function of position.
         */
        virtual
        double initial_value (const Point<dim> &position, const unsigned int n_field) const = 0;


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

    };




    /**
     * Register a volume-of-fluid initial conditions model so that it can be
     * selected from the parameter file.
     *
     * @param name A string that identifies the volume-of-fluid initial
     * conditions model
     * @param description A text description of what this model does and that
     * will be listed in the documentation of the parameter file.
     * @param declare_parameters_function A pointer to a function that can be
     * used to declare the parameters that this initial conditions model wants
     * to read from input files.
     * @param factory_function A pointer to a function that can create an
     * object of this initial conditions model.
     *
     * @ingroup VoFInitialConditionsModels
     */
    template <int dim>
    void
    register_initial_conditions_model (const std::string &name,
                                       const std::string &description,
                                       void (*declare_parameters_function) (ParameterHandler &),
                                       Interface<dim> *(*factory_function) ());

    /**
     * A function that given the name of a model returns a pointer to an
     * object that describes it. Ownership of the pointer is transferred to
     * the caller.
     *
     * The model object returned is not yet initialized and has not read its
     * runtime parameters yet.
     *
     * @ingroup VoFInitialConditionsModels
     */
    template <int dim>
    Interface<dim> *
    create_initial_conditions (ParameterHandler &prm);


    /**
     * Declare the runtime parameters of the registered volume-of-fluid initial
     * conditions models.
     *
     * @ingroup VoFInitialConditionsModels
     */
    template <int dim>
    void
    declare_parameters (ParameterHandler &prm);



    /**
     * Given a class name, a name, and a description for the parameter file for
     * a volume-of-fluid initial conditions model, register it with the
     * functions that can declare their parameters and create these objects.
     *
     * @ingroup VoFInitialConditionsModels
     */
#define ASPECT_REGISTER_VOF_INITIAL_CONDITIONS(classname,name,description) \
  template class classname<2>; \
  template class classname<3>; \
  namespace ASPECT_REGISTER_VOF_INITIAL_CONDITIONS_ ## classname \
  { \
    aspect::internal::Plugins::RegisterHelper<aspect::VoFInitialConditions::Interface<2>,classname<2> > \
    dummy_ ## classname ## _2d (&aspect::VoFInitialConditions::register_initial_conditions_model<2>, \
                                name, description); \
    aspect::internal::Plugins::RegisterHelper<aspect::VoFInitialConditions::Interface<3>,classname<3> > \
    dummy_ ## classname ## _3d (&aspect::VoFInitialConditions::register_initial_conditions_model<3>, \
                                name, description); \
  }
  }
}


#endif
