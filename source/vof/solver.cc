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

#include <aspect/global.h>
#include <aspect/vof/handler.h>

#include <deal.II/lac/constraint_matrix.h>

#ifdef ASPECT_USE_PETSC
#include <deal.II/lac/solver_cg.h>
#else
#include <deal.II/lac/trilinos_solver.h>
#endif

#include <deal.II/lac/pointer_matrix.h>

#include <deal.II/fe/fe_values.h>

namespace aspect
{
  template <int dim>
  void Simulator<dim>::VoFHandler::solve_vof_system ()
  {
    unsigned int block_idx = sim.introspection.variable("vofs").block_index;

    sim.computing_timer.enter_section ("   Solve VoF system");
    sim.pcout << "   Solving VoF system... " << std::flush;

    const double tolerance = std::max(1e-50,
                                      vof_solver_tolerance*sim.system_rhs.block(block_idx).l2_norm());

    SolverControl solver_control (1000, tolerance);

#ifdef ASPECT_USE_PETSC
    SolverCG<LinearAlgebra::Vector> solver(solver_control);
    LinearAlgebra::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#else
    TrilinosWrappers::SolverCG solver(solver_control);
    TrilinosWrappers::PreconditionJacobi precondition;
    precondition.initialize(sim.system_matrix.block(block_idx, block_idx));
#endif

    // Create distributed vector (we need all blocks here even though we only
    // solve for the current block) because only have a ConstraintMatrix
    // for the whole system, current_linearization_point contains our initial guess.
    LinearAlgebra::BlockVector distributed_solution (
      sim.introspection.index_sets.system_partitioning,
      sim.mpi_communicator);
    distributed_solution.block(block_idx) = sim.current_linearization_point.block (block_idx);

    sim.current_constraints.set_zero(distributed_solution);

    // solve the linear system:
    try
      {
        solver.solve (sim.system_matrix.block(block_idx,block_idx),
                      distributed_solution.block(block_idx),
                      sim.system_rhs.block(block_idx),
                      precondition);
      }
    // if the solver fails, report the error from processor 0 with some additional
    // information about its location, and throw a quiet exception on all other
    // processors
    catch (const std::exception &exc)
      {
        if (Utilities::MPI::this_mpi_process(sim.mpi_communicator) == 0)
          AssertThrow (false,
                       ExcMessage (std::string("The iterative advection solver "
                                               "did not converge. It reported the following error:\n\n")
                                   +
                                   exc.what()))
          else
            throw QuietException();
      }

    sim.current_constraints.distribute (distributed_solution);
    sim.solution.block(block_idx) = distributed_solution.block(block_idx);

    // print number of iterations and also record it in the
    // statistics file
    sim.pcout << solver_control.last_step()
              << " iterations." << std::endl;

    // Do not add VoF solver iterations to statistics, duplicaiton due to
    // splitting messes with file format
    // statistics.add_value("Iterations for VoF solver",
    //                      solver_control.last_step());

    sim.computing_timer.exit_section();
  }
}

namespace aspect
{
#define INSTANTIATE(dim) \
  template void Simulator<dim>::VoFHandler::solve_vof_system ();


  ASPECT_INSTANTIATE(INSTANTIATE)
}
