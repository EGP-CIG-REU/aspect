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

#ifndef __aspect__vofinterface_vof_system_h
#define __aspect__vofinterface_vof_system_h

#include <aspect/global.h>
#include <aspect/simulator_access.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  using namespace dealii;

  template <int dim>
  class Simulator;

  struct AdvectionField;

  namespace internal
  {
    namespace Assembly
    {
      namespace Scratch
      {
        template <int dim>
        struct VoFSystem
        {
          VoFSystem (const FiniteElement<dim> &finite_element,
                     const FiniteElement<dim> &vof_element,
                     const Mapping<dim>       &mapping,
                     const Quadrature<dim>    &quadrature,
                     const Quadrature<dim-1>  &face_quadrature);
          VoFSystem (const VoFSystem &data);

          FEValues<dim>          finite_element_values;
          FEFaceValues<dim>      face_finite_element_values;
          FEFaceValues<dim>      neighbor_face_finite_element_values;
          FESubfaceValues<dim>   subface_finite_element_values;

          std::vector<types::global_dof_index>   local_dof_indices;

          /**
           * Variables describing the values of the shape functions at the
           * quadrature points, as they are used in the advection assembly
           * function. note that the sizes of these arrays are equal to the
           * number of shape functions corresponding to the advected field (and
           * not of the overall field!), and that they are also correspondingly
           * indexed.
           */
          std::vector<double>         phi_field;
          std::vector<double>         face_phi_field;

          std::vector<double>         old_field_values;
          std::vector<Tensor<1,dim> > cell_i_n_values;
          std::vector<double>         cell_i_d_values;

          std::vector<Tensor<1,dim> > face_current_velocity_values;
          std::vector<Tensor<1,dim> > face_old_velocity_values;
          std::vector<Tensor<1,dim> > face_old_old_velocity_values;
          std::vector<Tensor<1,dim> > face_mesh_velocity_values;
        };
      }

      namespace CopyData
      {
        template <int dim>
        struct VoFSystem
        {
          /**
           * Constructor.
           * @param finite_element The element that describes the field for which we
           *    are trying to assemble a linear system. <b>Not</b> the global finite
           *    element.
           */
          VoFSystem(const FiniteElement<dim> &finite_element);
          VoFSystem(const VoFSystem &data);

          /**
           * Local contributions to the global matrix and right hand side
           * that correspond only to the variables listed in local_dof_indices
           */
          FullMatrix<double>          local_matrix;
          Vector<double>              local_rhs;
          /**
           * Local contributions to the global rhs from the face terms in the
           * discontinuous Galerkin interpretation of the VoF method.  The
           * vector is of length GeometryInfo<dim>::max_children_per_face *
           * GeometryInfo<dim>::faces_per_cell so as to hold one matrix for
           * each possible face or subface of the cell.
           **/
          std::vector<Vector<double>> local_f_rhs;

          /**
           * Denotes which face's rhs have actually been assembled in the DG
           * field assembly. Entries not used (for example, those corresponding
           * to non-existent subfaces; or faces being assembled by the
           * neighboring cell) are set to false.
           **/
          std::vector<bool>               assembled_rhs;

          /**
           * Indices of those degrees of freedom that actually correspond to
           * the vof field. Since this structure is used to represent just
           * contributions to the vof systems, there will be no contributions
           * to other parts of the system and consequently, we do not need to
           * list here indices that correspond to velocity or pressure degrees
           * (or, in fact any other variable outside the block we are currently
           * considering)
           */
          std::vector<types::global_dof_index>   local_dof_indices;
          /**
           * Indices of the degrees of freedom corresponding to the vof field
           * on all possible neighboring cells. This is used in the
           * discontinuous Galerkin interpretation of the VoF method. The outer
           * std::vector has length GeometryInfo<dim>::max_children_per_face *
           * GeometryInfo<dim>::faces_per_cell.
           **/
          std::vector<std::vector<types::global_dof_index> >   neighbor_dof_indices;
        };
      }
    }
  }
}


#endif
