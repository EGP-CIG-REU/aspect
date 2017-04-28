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


#ifndef __aspect__postprocess_visualization_vof_values_h
#define __aspect__postprocess_visualization_vof_values_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived that implements a funciton that provides volume of
       * fluid interface tracking data for graphical output.
       */
      template <int dim>
      class VoFValues
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          VoFValues ();

          virtual
          std::vector<std::string>
          get_names () const;

          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          virtual
          UpdateFlags
          get_needed_update_flags () const;

          virtual
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const;

          static
          void
          declare_parameters (ParameterHandler &prm);

          virtual
          void
          parse_parameters (ParameterHandler &prm);

        private:
          std::vector<std::string> vof_names;
          std::vector<DataComponentInterpretation::DataComponentInterpretation> interp;
          bool include_vofLS;
          bool include_vofN;
      };
    }
  }
}

#endif
