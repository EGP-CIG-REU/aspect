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
#include <aspect/vof/handler.h>
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
      AssertThrow(this->get_vof_handler().get_n_fields()!=numbers::invalid_unsigned_int,
                  ExcMessage("This postprocessor cannot be used without VoF fields"));
      const QGauss<dim> quadrature_formula (1);
      const unsigned int n_q_points = quadrature_formula.size();

      unsigned int n_vof_fields = this->get_vof_handler().get_n_fields();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);
      std::vector<double> vof_values(n_q_points);

      std::vector<double> local_vof_vol_sums(n_vof_fields);

      for (unsigned int f=0; f<n_vof_fields; ++f)
        {
          FEValuesExtractors::Scalar vof = this->get_vof_handler().get_field(f).fraction.extractor_scalar();
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
          local_vof_vol_sums[f] = vof_vol_sum;
        }

      std::vector<double> global_vof_vol_sums(n_vof_fields);

      Utilities::MPI::sum (local_vof_vol_sums, this->get_mpi_communicator(), global_vof_vol_sums);

      std::ostringstream output;
      output.precision(3);

      for (unsigned int f=0; f<n_vof_fields; ++f)
        {
          std::string col_name = "Global VoF volumes for " +
                                 this->get_vof_handler().get_field_name(f);
          statistics.add_value (col_name, global_vof_vol_sums[f]);

          output << global_vof_vol_sums[f];

          if (f+1 < n_vof_fields)
                 output << "/";

          // also make sure that the other columns filled by the this object
          // all show up with sufficient accuracy and in scientific notation
          {
            const std::string columns[] = { col_name
                                          };
            for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
              {
                statistics.set_precision (columns[i], 8);
                statistics.set_scientific (columns[i], true);
              }
          }
        }

      return std::pair<std::string, std::string> ("Global VoF volumes (m^3):",
                                                  output.str());
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
