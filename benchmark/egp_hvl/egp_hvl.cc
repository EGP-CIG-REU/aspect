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
  namespace InclusionBenchmark
  {
    using namespace dealii;

    namespace AnalyticSolutions
    {
      // based on http://geodynamics.org/hg/cs/AMR/Discontinuous_Stokes with permission
      // The following code has been taken from http://www.underworldproject.org/,
      // release 1.7.0. As mentioned in the Underworld Manual, this code has been
      // released under the GNU General Public License (GPL).

      void _Velic_solCx(
        double pos[],
        double _eta_A, double _eta_B,   /* Input parameters: density, viscosity A, viscosity B */
        double _x_c, int _n,      /* Input parameters: viscosity jump location, wavenumber in x */
        double vel[], double *presssure,
        double total_stress[], double strain_rate[] )
      {
        /****************************************************************************************/
        /****************************************************************************************/
        /* Output */
        vel[0] = -pos[1];  
        vel[1] = pos[0];

        (*presssure) = 1;

        total_stress[0] = 0.0;
        total_stress[1] = 0.0;
        total_stress[2] = 0.0;

        strain_rate[0] = 0;
        strain_rate[1] = 0;
        strain_rate[2] = 0;
      }

      /**
       * The exact solution for the SolCx benchmark, given the value
       * of the jump in viscosity $\eta_B$.
       */
      template <int dim>
      class FunctionSolCx : public Function<dim>
      {
        public:
          FunctionSolCx (const double eta_B,
                         const double background_density)
            :
            Function<dim>(),
            eta_B_(eta_B),
            background_density (background_density)
          {}

          virtual void vector_value (const Point< dim >   &p,
                                     Vector< double >   &values) const
          {
            AssertDimension(values.size(), 4);

            double pos[2]= {p(0),p(1)};
            double total_stress[3], strain_rate[3];
            double eta_A=1.0;
            double eta_B=1.0;

            // call the analytic function for the solution with a zero
            // background density
            AnalyticSolutions::_Velic_solCx
            (pos,
             eta_A, eta_B,
             0.5, 1,
             &values[0], &values[2], total_stress, strain_rate );

            // then add the background pressure to the value we just got
            // values[2] += (0.5-p[1])*background_density;
            values[2] = 0.0;
          }
        private:
          double eta_B_, background_density;
      };
    }

    /**
      * A material model that describes the <i>SolCx</i> benchmark of the
      * paper cited in the documentation of the DuretzEtAl namespace.
      *
      * @note The SolCx benchmark only talks about the flow field, not about
      * a temperature field. All quantities related to the temperature are
      * therefore set to zero in the implementation of this class.
      *
      * @note The analytic solution of this benchmark is implemented in the
      * "SolCx error" postprocessor in
      * aspect::Postprocessor::DuretzEtAl::SolCx class and can be used to
      * assess the accuracy of the computed solution.
      *
      * @ingroup MaterialModels
      */
    template <int dim>
    class SolCxMaterial : public MaterialModel::InterfaceCompatibility<dim>
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
    SolCxMaterial<dim>::
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
    SolCxMaterial<dim>::
    reference_viscosity () const
    {
      return 1;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    reference_density () const
    {
      return background_density;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    reference_thermal_expansion_coefficient () const
    {
      return 0;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    specific_heat (const double,
                   const double,
                   const std::vector<double> &, /*composition*/
                   const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    reference_cp () const
    {
      return 0;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    thermal_conductivity (const double,
                          const double,
                          const std::vector<double> &, /*composition*/
                          const Point<dim> &) const
    {
      return 0;
    }

    template <int dim>
    double
    SolCxMaterial<dim>::
    reference_thermal_diffusivity () const
    {
      return 0;
    }

/**
  TODO: A debugging idea is that we write another identical benchmark setups where in one case density is overridden by interpolated value of density and in this case with the analytical solution.
**/
    template <int dim>
    double
    SolCxMaterial<dim>::
    density (const double,
             const double,
             const std::vector<double> &, /*composition*/
             const Point<dim> &p) const
    {
     /** TODO: Rename background_density to reference_density
     **/
      return background_density*(1-p[1])/2.0;
    }


    template <int dim>
    double
    SolCxMaterial<dim>::
    thermal_expansion_coefficient (const double temperature,
                                   const double,
                                   const std::vector<double> &, /*composition*/
                                   const Point<dim> &) const
    {
      return 0;
    }


    template <int dim>
    double
    SolCxMaterial<dim>::
    compressibility (const double,
                     const double,
                     const std::vector<double> &, /*composition*/
                     const Point<dim> &) const
    {
      return 0.0;
    }


    template <int dim>
    bool
    SolCxMaterial<dim>::
    is_compressible () const
    {
      return false;
    }

    template <int dim>
    void
    SolCxMaterial<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("SolCx");
        {
          prm.declare_entry ("Viscosity jump", "1",
                             Patterns::Double (0),
                             "Viscosity in the right half of the domain.");
          /** TODO: Rename Background density to reference density
          **/
          prm.declare_entry ("Background density", "1",
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
    SolCxMaterial<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("SolCx");
        {
          eta_B = prm.get_double ("Viscosity jump");
          background_density = prm.get_double("Background density");
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
    SolCxMaterial<dim>::get_eta_B() const
    {
      return eta_B;
    }


    template <int dim>
    double
    SolCxMaterial<dim>::get_background_density() const
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
    class SolCxPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
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
    SolCxPostprocessor<dim>::execute (TableHandler &statistics)
    {
      std_cxx1x::shared_ptr<Function<dim> > ref_func;
      if (dynamic_cast<const SolCxMaterial<dim> *>(&this->get_material_model()) != NULL)
        {
          const SolCxMaterial<dim> *
          material_model
            = dynamic_cast<const SolCxMaterial<dim> *>(&this->get_material_model());

          ref_func.reset (new AnalyticSolutions::FunctionSolCx<dim>(material_model->get_eta_B(),
                                                                    material_model->get_background_density()));
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Postprocessor DuretzEtAl only works with the material model SolCx, SolKz, and Inclusion."));
        }

      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.velocities).degree+2);

      Vector<float> cellwise_errors_u (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_p (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_ul2 (this->get_triangulation().n_active_cells());
      Vector<float> cellwise_errors_pl2 (this->get_triangulation().n_active_cells());

      ComponentSelectFunction<dim> comp_u(std::pair<unsigned int, unsigned int>(0,dim),
                                          this->get_fe().n_components());
      ComponentSelectFunction<dim> comp_p(dim, this->get_fe().n_components());

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
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SolCxMaterial,
                                   "SolCxMaterial",
                                   "A material model that corresponds to the 'SolCx' benchmark "
                                   "defined in Duretz et al., G-Cubed, 2011.")


    ASPECT_REGISTER_POSTPROCESSOR(SolCxPostprocessor,
                                  "SolCxPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "the Duretz et al., G-Cubed, 2011, paper with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")


  }
}
