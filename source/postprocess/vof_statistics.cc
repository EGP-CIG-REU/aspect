/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

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


#include <aspect/postprocess/vof_statistics.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VoFStatistics<dim>::execute (TableHandler &statistics)
    {
      const QGauss<dim> quadrature_formula (1);
      const unsigned int n_q_points = quadrature_formula.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<double> vof_values(n_q_points);
      FEValuesExtractors::Scalar vof = this->introspection().variable("vofs").extractor_scalar();

      double vof_vol_sum=0.0, vof_vol_corr=0.0;

      typename DoFHandler<dim>::active_cell_iterator
      cell = this->get_dof_handler().begin_active(),
      endc = this->get_dof_handler().end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[vof].get_function_values (this->get_solution(),
                                                vof_values);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                //Use Kahan sum for improved consistency
                double cell_term = vof_values[q] * fe_values.JxW(q);
                double c_nterm = cell_term - vof_vol_corr;
                double nsum = vof_vol_sum + c_nterm;
                vof_vol_corr = (nsum - vof_vol_sum) - c_nterm;
                vof_vol_sum = nsum;
              }
          }

      const double global_vof_sum
        = Utilities::MPI::sum (vof_vol_sum, this->get_mpi_communicator());

      std::ostringstream output;
      output.precision(3);

      statistics.add_value ("Global VoF vol", global_vof_sum);

      output << global_vof_sum
             << " m^3";

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Global VoF vol"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }
      }

      if (log_vol_delta)
        {
          if (!initial_vol_computed)
            {
              initial_vol_computed = true;
              init_vol = vof_vol_sum;
            }
          statistics.add_value ("Delta VoF vol",
                                vof_vol_sum-init_vol);

          output << " Delta "
                 << (global_vof_sum-init_vol)
                 << " m^3";

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const char *columns[] = { "Delta VoF vol"
                                    };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      return std::pair<std::string, std::string> ("Global VoF vol:",
                                                  output.str());
    }

    template <int dim>
    void
    VoFStatistics<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("VoF statistics");
        {
          prm.declare_entry ("Log volume change", "false",
                             Patterns::Bool (),
                             "Option to also log change in Global VoF volume in statistics file.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    VoFStatistics<dim>::parse_parameters (ParameterHandler &prm)
    {
      initial_vol_computed = false;
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("VoF statistics");
        {
          log_vol_delta = prm.get_bool ("Log volume change");
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
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VoFStatistics,
                                  "vof statistics",
                                  "A postprocessor that computes some statistics about the "
                                  "vof field.")
  }
}
