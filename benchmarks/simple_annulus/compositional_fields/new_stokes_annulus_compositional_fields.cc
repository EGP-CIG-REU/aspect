#include "new_stokes_annulus_compositional_fields.h"

// explicit instantiations
namespace aspect
{
  namespace InclusionBenchmark
  {
    ASPECT_REGISTER_MATERIAL_MODEL(SimpleAnnulusCompositionalMaterialModel,
                                   "SimpleAnnulusCompositionalMaterial",
                                   "Material model for the simple annulus benchmark.")


    ASPECT_REGISTER_POSTPROCESSOR(SimpleAnnulusPostprocessor,
                                  "SimpleAnnulusPostprocessor",
                                  "A postprocessor that compares the solution of the benchmarks from "
                                  "derived analytical solution with the one computed by ASPECT "
                                  "and reports the error. Specifically, it can also compute the errors for "
                                  "the SolCx, SolKz and inclusion benchmarks. The postprocessor inquires "
                                  "which material model is currently being used and adjusts "
                                  "which exact solution to use accordingly.")


  }
}
