#include <aspect/material_model/simple.h>
#include <aspect/boundary_velocity/interface.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/global.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <aspect/initial_composition/function.h>
#include <deal.II/base/parsed_function.h>


namespace aspect
{
  /**
   * This is the "simple_annulus" benchmark defined in the following paper:
   * @code
   *  @Article{thph17,
   *    author =       {C. Thieulot and E.G.P. Puckett and H. Lokavarapu},
   *    title =        {Stokes flow in an annulus: analytical solution and numerical benchmark},
   *    journal =      {xxx},
   *    year =         2017,
   *    volume =       xxx,
   *    number =       {x},
   *    publisher =    {xxx},
   *    pages =        {xxx--xxx}}
   * @endcode
   *
   */
  namespace InclusionBenchmark
  {
    using namespace dealii;


    namespace AnalyticSolutions
    {

      template <int dim>
      void analytic_solution(
        double pos[],
        double vel[],
        double *pressure,
        std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > pressure_function,
        std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > velocity_function)
      {
        /****************************************************************************************/
        /****************************************************************************************/
        /* Output */
        for (unsigned int i=0; i < dim; i++)
          vel[i] = velocity_function->value(Point<dim>(pos[0],pos[1]), i);

        //(*pressure) = 0;
        (*pressure) = pressure_function->value(Point<dim>(pos[0], pos[1]));
      }

      /**
       * The exact solution for the benchmark, given
       * density $\rho$.
       */
      template <int dim>
      class FunctionStreamline : public Function<dim>
      {
        public:
          FunctionStreamline (std_cxx11::shared_ptr<Functions::ParsedFunction<dim>> pressure,
                              std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > velocity)
            :
            Function<dim>(),
            pressure_function (pressure),
            velocity_function (velocity)
          {}

          virtual void vector_value (const Point< dim > &p,
                                     Vector< double >   &values) const
          {
            double pos[2]= {p(0),p(1)};

            AnalyticSolutions::analytic_solution<dim>
            (pos,
             &values[0], &values[2], pressure_function, velocity_function);
          }
        private:
          std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > pressure_function;
          std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > velocity_function;
      };
    }

    template <int dim>
    class SimpleAnnulusMaterialModel : public MaterialModel::Interface<dim>
    {
      private:
        std_cxx11::shared_ptr<Functions::ParsedFunction<dim> > density_function, pressure_function, velocity_function;

      public:
        /**
         * @name Physical parameters used in the basic equations
         * @{
         */
        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.densities[i] = density_function->value(in.position[i]);
              std::ostringstream os1;
              out.viscosities[i] = 1;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;

            }
        }

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
        virtual bool is_compressible () const
        {
          return false;
        }
        /**
         * @}
         */


        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Simple annulus");
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

              prm.enter_subsection("Density");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Pressure");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical velocity");
              {
                prm.enter_subsection("Function");
                {
                  Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }

        /**
         * Read the parameters this class declares from the parameter file.
         */
        virtual
        void
        parse_parameters (ParameterHandler &prm)
        {
          prm.enter_subsection("Material model");
          {
            prm.enter_subsection("Simple annulus");
            {
              eta_B = prm.get_double ("Viscosity jump");
              background_density = prm.get_double("Reference density");

              prm.enter_subsection("Density");
              {
                prm.enter_subsection("Function");
                {
                  density_function.reset(new Functions::ParsedFunction<dim>(1));
                  density_function->parse_parameters(prm);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Pressure");
              {
                prm.enter_subsection("Function");
                {
                  pressure_function.reset(new Functions::ParsedFunction<dim>(1));
                  pressure_function->parse_parameters(prm);
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();

              prm.enter_subsection("Analytical velocity");
              {
                prm.enter_subsection("Function");
                {
                  try
                    {
                      velocity_function.reset (new Functions::ParsedFunction<dim>(dim));
                      velocity_function->parse_parameters (prm);
                    }
                  catch (...)
                    {
                      std::cerr << "ERROR: FunctionParser failed to parse\n"
                                << "\t'Analytical velocity.Function'\n"
                                << "with expression\n"
                                << "\t'" << prm.get("Function expression") << "'\n"
                                << "More information about the cause of the parse error \n"
                                << "is shown below.\n";
                      throw;
                    }
                }
                prm.leave_subsection();
              }
              prm.leave_subsection();
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


        /**
         * @name Reference quantities
         * @{
         */
        virtual double reference_viscosity () const
        {
          return 1;
        }

        /**
         * Returns the viscosity value on the right half of the domain,
         * typically 1 or 1e6
         */
        double get_eta_B() const
        {
          return eta_B;
        }
        /**
         * Returns the background density of this model. See the
         * corresponding member variable of this class for more information.
         */

        std_cxx11::shared_ptr<Functions::ParsedFunction<dim>>  get_pressure() const
        {
          return pressure_function;
        }

        std_cxx11::shared_ptr<Functions::ParsedFunction<dim>>  get_velocity() const
        {
          return velocity_function;
        }

        double get_background_density() const
        {
          return background_density;
        }

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


    /**
      * A postprocessor that evaluates the accuracy of the solution.
      *
      */
    template <int dim>
    class SimpleAnnulusPostprocessor : public Postprocess::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Generate graphical output from the current solution.
         */
        virtual
        std::pair<std::string,std::string>
        execute (TableHandler &)
        {
          std_cxx1x::shared_ptr<Function<dim> > ref_func;
          if (dynamic_cast<const SimpleAnnulusMaterialModel<dim> *>(&this->get_material_model()) != NULL)
            {
              const SimpleAnnulusMaterialModel<dim> *
              material_model
                = dynamic_cast<const SimpleAnnulusMaterialModel<dim> *>(&this->get_material_model());

              ref_func.reset (new AnalyticSolutions::FunctionStreamline<dim>(material_model->get_pressure(), material_model->get_velocity()));
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
    };
  }
}



