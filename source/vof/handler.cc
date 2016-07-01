
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

#include <aspect/simulator.h>
#include <aspect/global.h>
#include <aspect/vof/handler.h>

#include <deal.II/fe/fe_dgq.h>

using namespace dealii;

namespace aspect
{

  template <int dim>
  Simulator<dim>::VoFHandler::VoFHandler (Simulator<dim> &simulator,
                                          ParameterHandler &prm)
    : sim (simulator)
  {
    //parse_parameters (prm);

    sim.signals.edit_finite_element_variables.connect(std_cxx11::bind(&aspect::Simulator<dim>::VoFHandler::edit_finite_element_variables,
                                                                      std_cxx11::ref(*this),
                                                                      std_cxx11::_1));
  }

  template <int dim>
  void
  Simulator<dim>::VoFHandler::edit_finite_element_variables (std::vector<VariableDeclaration<dim> > &vars)
  {
    vars.push_back(VariableDeclaration<dim>("vofs",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(0)),
                                            1,
                                            1));

    vars.push_back(VariableDeclaration<dim>("vofsN",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(0)),
                                            dim+1,
                                            1));

    vars.push_back(VariableDeclaration<dim>("vofsLS",
                                            std_cxx11::shared_ptr<FiniteElement<dim>>(
                                              new FE_DGQ<dim>(1)),
                                            1,
                                            1));
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template class Simulator<dim>::VoFHandler;

  ASPECT_INSTANTIATE(INSTANTIATE)
}
