/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/vof/handler.h>
#include <aspect/vof_initial_conditions/function.h>
#include <aspect/postprocess/interface.h>

namespace aspect
{
  namespace VoFInitialConditions
  {
    template <int dim>
    Function<dim>::Function ()
      :
      n_init_samples (3)
    {}

    template <int dim>
    unsigned int Function<dim>::n_samples () const
    {
      return n_init_samples;
    }

    template <int dim>
    typename VoFInitType::Kind Function<dim>::init_type () const
    {
      return function_init_type;
    }

    template <int dim>
    double
    Function<dim>::
    initial_value (const Point<dim> &position, const unsigned int n_field) const
    {
      return function->value(position, n_field);
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("VoF initial conditions");
      {
        prm.enter_subsection("Function");
        {
          prm.declare_entry("Signed distance init", "false",
                            Patterns::Bool (),
                            "When set to true, initialization will be assumed to be a"
                            "signed distance level set function.");

          prm.declare_entry ("Number initialization samples", "3",
                             Patterns::Integer (1),
                             "Number of sampled points per dimension when initializing from VOF");

          Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Function<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("VoF initial conditions");
      {
        prm.enter_subsection("Function");
        bool is_dist_init = prm.get_bool("Signed distance init");

        if (is_dist_init)
          function_init_type = VoFInitType::signed_distance_level_set;
        else
          function_init_type = VoFInitType::composition;

        n_init_samples = prm.get_integer ("Number initialization samples");

        try
          {
            function.reset(new Functions::ParsedFunction<dim>(this->get_vof_handler().get_n_fields()));
            function->parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'VoF initial conditions.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'\n"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
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
  namespace VoFInitialConditions
  {
    ASPECT_REGISTER_VOF_INITIAL_CONDITIONS(Function,
                                           "function",
                                           "Specify the composition in terms of an explicit formula. The format of these "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
