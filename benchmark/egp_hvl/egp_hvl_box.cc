#include <aspect/material_model/simple.h>
#include <aspect/velocity_boundary_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>



namespace aspect
{
  /**
   * This is the "Sol Cx" benchmark defined in the following paper:
   * @code
   *  @Article{DMGT11,
   *    author =       {T. Duretz and D. A. May and T. V. Gerya and P. J. Tackley},
   *    title =        {Discretization errors and free surface stabilization in the
   *                  finite difference and marker-in-cell method for applied
   *                  geodynamics: {A} numerical study},
   *    journal =      {Geochemistry Geophysics Geosystems},
   *    year =         2011,
   *    volume =       12,
   *    pages =        {Q07004/1--26}}
   * @endcode
   *
   * The results are published in Kronbichler, Heister and Bangerth paper.
   */
  namespace EGPHVLBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
      // The following code has been taken from http://www.underworldproject.org/,
      // release 1.7.0. As mentioned in the Underworld Manual, this code has been
      // released under the GNU General Public License (GPL).

      void analytic_solution(
        double pos[],  /* Input parameters: density, viscosity A, viscosity B */
        double vel[], double *presssure)
      {
        /****************************************************************************************/
        /****************************************************************************************/
        /* Output */
        vel[0] = -1*cos(pos[1])*sin(pos[0]);
        vel[1] =  cos(pos[0])*sin(pos[1]);

        (*presssure) = cos(pos[0])*cos(pos[1]);
      }

      /**
       * The exact solution for the SolCx benchmark, given the value
       * of the jump in viscosity $\eta_B$.
       */
      template <int dim>
      class FunctionEGPHVL : public Function<dim>
      {
        public:
          FunctionEGPHVL ()
            :
            Function<dim>()
          {}

          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            AssertDimension(values.size(), 4);

            double pos[2]= {p(0),p(1)};
            double total_stress[3], strain_rate[3];

            // call the analytic function for the solution with a zero
            // background density
            AnalyticSolutions::analytic_solution
            (pos,
             &values[0], &values[2]);
          }
      };
    }

    template <int dim>
    class BenchmarkMaterialModel : public MaterialModel::InterfaceCompatibility<dim>
    {
      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;

        virtual double density (const double temperature,
                                const double pressure,
                                const std::vector<double> &compositional_fields,
                                const Point<dim> &position) const;

        virtual double compressibility (const double temperature,
                                        const double pressure,
                                        const std::vector<double> &compositional_fields,
                                        const Point<dim> &position) const;

        virtual double specific_heat (const double temperature,
                                      const double pressure,
                                      const std::vector<double> &compositional_fields,
                                      const Point<dim> &position) const;

        virtual double thermal_expansion_coefficient (const double      temperature,
                                                      const double      pressure,
                                                      const std::vector<double> &compositional_fields,
                                                      const Point<dim> &position) const;

        virtual double thermal_conductivity (const double temperature,
                                             const double pressure,
                                             const std::vector<double> &compositional_fields,
                                             const Point<dim> &position) const;
        /**
         * @}
         */

        /**
         * @name Qualitative properties one can ask a material model
         * @{
         */


        /**
         * Return whether the model is compressible or not.
         * Incompressibility does not necessarily imply that the density is
         * constant; rather, it may still depend on temperature or pressure.
         * In the current context, compressibility means whether we should
         * solve the contuity equation as $\nabla \cdot (\rho \mathbf u)=0$
         * (compressible Stokes) or as $\nabla \cdot \mathbf{u}=0$
         * (incompressible Stokes).
         */
        virtual bool is_compressible () const;
        /**
         * @}
         */


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm);


        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const;

        virtual double reference_density () const;

        virtual double reference_thermal_expansion_coefficient () const;

//TODO: should we make this a virtual function as well? where is it used?
        double reference_thermal_diffusivity () const;

        double reference_cp () const;
        /**
         * @}
         */

        /**
         * Returns the viscosity value on the right half of the domain,
         * typically 1 or 1e6
         */
        double get_eta_B() const;

        /**
         * Returns the background density of this model. See the
         * corresponding member variable of this class for more information.
         */
        double get_background_density() const;

      private:
        /**
         * Viscosity value on the right half of the domain, typically 1 or
         * 1e6
         */
        double eta_B;

        /**
         * A constant background density over which the density variations
         * are overlaid. This constant density has no effect on the dynamic
         * pressure and consequently on the flow field, but it contributes
         * to the total pressure via the adiabatic pressure. We use this
         * field to support our claim in the first ASPECT paper that the
         * accuracy of the solutions is guaranteed even if we don't subtract
         * the adiabatic pressure in our computations.
         */
        double background_density;
    };






    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    viscosity (const double,
               const double,
               const std::vector<double> &,       /*composition*/
               const SymmetricTensor<2,dim> &,
               const Point<dim> &p) const
    {
      // defined as given in the Duretz et al. paper
      return 1;
    }


    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    reference_viscosity () const
    {
      return 1;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    reference_density () const
    {
      return background_density;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    reference_cp () const
    {
      return 0;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    reference_thermal_diffusivity () const
    {
      return 0;
    }

    /**
      TODO: A debugging idea is that we write another identical benchmark setups where in one case density is overridden by interpolated value of density and in this case with the analytical solution.
    **/
    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &p) const
    {
      /** TODO: Rename background_density to reference_density
      **/
      // This function creates a closed circuit flow for the 2d box geometry model.
      //return 2.0 * p[1];
      // Use this function if running 2d spherical shell geometry model.
      double r = p.norm();
      double density = r*r;
      return 1;
    }


    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    BenchmarkMaterialModel<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }


    template <int dim>
    bool
    BenchmarkMaterialModel<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    BenchmarkMaterialModel<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("EGPHVL");
        {
          prm.declare_entry ("Viscosity jump", "1",
                             Patterns::Double (0),
                             "Viscosity in the right half of the domain.");
          prm.declare_entry ("Reference density", "1",
                             Patterns::Double (0),
                             "Density value upon which the variation of this testcase "
                             "is overlaid. Since this background density is constant "
                             "it does not affect the flow pattern but it adds to the "
                             "total pressure since it produces a nonzero adiabatic "
                             "pressure if set to a nonzero value.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    BenchmarkMaterialModel<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("EGPHVL");
        {
          eta_B = prm.get_double ("Viscosity jump");
          background_density = prm.get_double("Reference density");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = MaterialModel::NonlinearDependence::none;
      this->model_dependence.density = MaterialModel::NonlinearDependence::none;
      this->model_dependence.compressibility = MaterialModel::NonlinearDependence::none;
      this->model_dependence.specific_heat = MaterialModel::NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = MaterialModel::NonlinearDependence::none;
    }

    template <int dim>
    double
    BenchmarkMaterialModel<dim>::get_eta_B() const
    {
      return eta_B;
    }


    template <int dim>
    double
    BenchmarkMaterialModel<dim>::get_background_density() const
    {
      return background_density;
    }




    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      * The implementation of error evaluators that correspond to the
      * benchmarks defined in the paper Duretz et al. reference above.
      */
    template <int dim>
    class EGPHVLPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &statistics);
    };

    template <int dim>
    std::pair<std::string,std::string>
    EGPHVLPostprocessor<dim>::execute (TableHandler &statistics)
    {
      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      if (dynamic_cast<const BenchmarkMaterialModel<dim> *>(&this->get_material_model()) != NULL)
        {
          const BenchmarkMaterialModel<dim> *
          material_model
            = dynamic_cast<const BenchmarkMaterialModel<dim> *>(&this->get_material_model());

          ref_func.reset (new AnalyticSolutions::FunctionEGPHVL<dim>());
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor DuretzEtAl only works with the material model SolCx, SolKz, and Inclusion."));
        }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);
#ifdef DEBUG
     /* 
            FEValues<dim> fe_values(this->get_mapping(),
                                    this->get_fe(),
                                    quadrature_formula,
                                    update_values |
                                    update_quadrature_points |
                                    update_JxW_values);

          const FEValuesExtractors::Scalar &extractor_pressure = this->introspection().extractors.pressure;
          const FEValuesExtractors::Scalar &extractor_temperature = this->introspection().extractors.temperature;
            const FEValuesExtractors::Vector &extractor_velocity = this->introspection().extractors.velocities;

            unsigned int quadrature_point_index = 0;

      typename DoFHandler<dim>::active_cell_iterator cell = this->get_dof_handler().begin_active();
      for (; cell != this->get_dof_handler().end();
             ++cell) {
        if (!(cell->is_locally_owned()))
          continue;

        fe_values.reinit(cell);
        // Each vector is of the same length.
        const std::vector<Point<dim>> quadrature_points = fe_values.get_quadrature_points();

        std::vector<double> temperature(fe_values.n_quadrature_points);
        std::vector<double> pressure(fe_values.n_quadrature_points);
        std::vector<Tensor<1, dim>> velocity(fe_values.n_quadrature_points);
        std::vector<std::vector<double>> compositional_fields(this->n_compositional_fields(),
                                                              std::vector<double>(fe_values.n_quadrature_points));

        const std::vector<double> jacobian_weights = fe_values.get_JxW_values();

        fe_values[extractor_pressure].get_function_values(this->get_solution(), pressure);
        fe_values[extractor_temperature].get_function_values(this->get_solution(), temperature);
        fe_values[extractor_velocity].get_function_values(this->get_solution(), velocity);

        unsigned int quadrature_point_index = 0;
        for (;quadrature_point_index < fe_values.n_quadrature_points; quadrature_point_index++)
        {
          unsigned int my_rank = Utilities::MPI::this_mpi_process(this->get_mpi_communicator());
          std::ofstream solution_log("computed_solution_" + Utilities::int_to_string(my_rank) + ".log", std::ios_base::app);
          solution_log << std::setprecision(14) << quadrature_points[quadrature_point_index][0] << " "
                                                << quadrature_points[quadrature_point_index][1] << " "
                                                << velocity[quadrature_point_index][0] << " "
                                                << velocity[quadrature_point_index][1] << " "
                                                << temperature[quadrature_point_index] << " "
                                                << pressure[quadrature_point_index] << std::endl;
          solution_log.close();
        }
      }

      std::ofstream stream_solution("solution.dat", std::ios_base::app);
      this->get_solution().print(stream_solution, 14);
     */
#endif
      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          dim+2);
      ComponentSelectFunction<dim> comp_p(dim, dim+2);

      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_u,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_p,
                                         quadrature_formula,
                                         VectorTools::L1_norm,
                                         &comp_p);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_ul2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_u);
      VectorTools::integrate_difference (this->get_mapping(),this->get_dof_handler(),
                                         this->get_solution(),
                                         *ref_func,
                                         cellwise_errors_pl2,
                                         quadrature_formula,
                                         VectorTools::L2_norm,
                                         &comp_p);

      const double u_l1 = Utilities::MPI::sum(cellwise_errors_u.l1_norm(),this->get_mpi_communicator());
      const double p_l1 = Utilities::MPI::sum(cellwise_errors_p.l1_norm(),this->get_mpi_communicator());
      const double u_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_ul2.norm_sqr(),this->get_mpi_communicator()));
      const double p_l2 = std::sqrt(Utilities::MPI::sum(cellwise_errors_pl2.norm_sqr(),this->get_mpi_communicator()));

      std::ostringstream os;
      os << std::scientific << u_l1
         << ", " << p_l1
         << ", " << u_l2
         << ", " << p_l2;

      return std::make_pair("Errors u_L1, p_L1, u_L2, p_L2:", os.str());
    }

  }
}



// explicit instantiations
namespace aspect
{
  namespace EGPHVLBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(BenchmarkMaterialModel,
                                   "EGPHVLMaterial",
                                   "EGP and HVL benchmark material model.")


    ASPECT_REGISTER_POSTPROCESSOR(EGPHVLPostprocessor,
                                  "EGPHVLPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "derived analytical solution with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can also compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")


  }
}
