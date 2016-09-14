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


#include <aspect/vof_initial_conditions/composition.h>
#include <aspect/postprocess/interface.h>

namespace aspect
{
  namespace VoFInitialConditions
  {
    template <int dim>
    Function<dim>::Function ()
      :
      n_init_samples (3),
      field_name ("")
    {}

    template <int dim>
    unsigned int Function<dim>::n_samples () const
    {
      return n_init_samples;
    }

    template <int dim>
    typename VoFInitType::Kind Function<dim>::init_type () const
    {
      return VoFInitType::composition;
    }

    template <int dim>
    double
    Function<dim>::
    initial_value (const Point<dim> &position) const
    {
      const CompositionalInitialConditions::Interface c_init = this->get_compositional_initial_conditions();
      return c_init.initial_composition(position, c_idx);
    }

    template <int dim>
    void
    Function<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("VoF initial conditions");
      {
        prm.enter_subsection("Composition");
        {
          prm.declare_entry ("Number initialization samples", "3",
                             Patterns::Integer (1),
                             "Number of sampled points per dimension when initializing from VOF");

          prm.declare_entry ("Composition field", "",
                             Patterns::Anything(),
                             "Name of composition field to sample");
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
        {
          n_init_samples = prm.get_integer ("Number initialization samples");

          c_idx = this->introspection().compositional_index_for_name(prm.get("Composition field"));
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
    ASPECT_REGISTER_VOF_INITIAL_CONDITIONS(Composition,
                                           "composition",
                                           "Specify the VoF values by referencing data from the compositional initial conditions")
  }
}
