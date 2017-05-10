/*
  Copyright (C) 2013 by the authors of the ASPECT code.

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


#ifndef __aspect_mesh_refinement_vof_boundary_h
#define __aspect_mesh_refinement_vof_boundary_h

#include <aspect/mesh_refinement/interface.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace MeshRefinement
  {

    /**
     * A class that implements a mesh refinement criterion that refines the
     * mesh near boundaries for the VoF interface tracking algorithm.
     *
     * @ingroup MeshRefinement
     */
    template <int dim>
    class VoFBoundary : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Mark large unrefined neighboring cells for refinement and prevent
         * coarsening
         */
        virtual
        void
        tag_additional_cells() const;


        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        int min_interface_level;
        double vof_epsilon;
    };
  }
}

#endif
