#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, Hermes1DFunction<double>* lambda_al,
                                             std::string mat_cu, Hermes1DFunction<double>* lambda_cu,
                                             Hermes2DFunction<double>* src_term) : WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new DefaultMatrixFormDiffusion<double>(0, 0, mat_al, lambda_al));
  add_matrix_form(new DefaultMatrixFormDiffusion<double>(0, 0, mat_cu, lambda_cu));

  // Residual forms.
  add_vector_form(new DefaultVectorFormVol<double>(0, HERMES_ANY, src_term));
};
