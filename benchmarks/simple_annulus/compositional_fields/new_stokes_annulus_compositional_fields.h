#ifndef NEW_STOKES_ANNULUS_COMPOSITIONAL_FIELDS_H
#define NEW_STOMES_ANNULUS_COMPOSITIONAL_FIELDS_H

#include "../new_stokes_annulus.h"

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

    template <int dim>
    class SimpleAnnulusCompositionalMaterialModel : public SimpleAnnulusMaterialModel<dim>
    {
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
              out.densities[i] = in.composition[i][0];
              out.viscosities[i] = 1;
              out.compressibilities[i] = 0;
              out.specific_heat[i] = 0;
              out.thermal_expansion_coefficients[i] = 0;
              out.thermal_conductivities[i] = 0.0;

            }
        }
    };
 }
}
#endif


