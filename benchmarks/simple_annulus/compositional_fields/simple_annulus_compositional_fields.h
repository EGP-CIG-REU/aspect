#ifndef ASPECT_SIMPLE_ANNULUS_COMPOSITIONAL_FIELDS_H
#define ASPECT_SIMPLE_ANNULUS_COMPOSITIONAL_FIELDS_H

#include "../simple_annulus.h"



namespace aspect
{
  namespace InclusionBenchmark
  {
    using namespace dealii;

    template<int dim>
    class SimpleAnnulusCompositionalMaterialModel : public SimpleAnnulusMaterialModel<dim>
    {
      public:
        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const
        {
          for (unsigned int i=0; i < in.position.size(); ++i)
            {
              out.viscosities[i] = 1;
              out.densities[i] = in.composition[i][0];
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

