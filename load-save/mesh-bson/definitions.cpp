#include "definitions.h"

/* Weak forms */

CustomWeakForm::CustomWeakForm() : WeakForm<double>(2)
{
  // Jacobian forms.
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(1, 1));

  // Residual forms.
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY, new Hermes2DFunction<double>(100.0)));
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(1, Hermes::HERMES_ANY, new Hermes2DFunction<double>(1.0)));
};

CustomWeakForm::CustomFormAdvection::CustomFormAdvection(int i, int j, std::string area) : MatrixFormVol<double>(i, j)
{
  this->set_area(area);
}

double CustomWeakForm::CustomFormAdvection::value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
  Func<double> *v, GeomVol<double> *e, Func<double> **ext) const
{
  double result = 0;
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (ext[0]->dx[i] * u->dx[i] + ext[0]->dy[i] * u->dy[i]) * v->val[i];
  }
  return result;
}

Ord CustomWeakForm::CustomFormAdvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
  GeomVol<Ord> *e, Func<Ord> **ext) const
{
  Ord result = Ord(0);
  for (int i = 0; i < n; i++) 
  {
    result += wt[i] * (ext[0]->dx[i] * u->dx[i] + ext[0]->dy[i] * u->dy[i]) * v->val[i];
  }
  return result;
}

MatrixFormVol<double>* CustomWeakForm::CustomFormAdvection::clone() const
{
  return new CustomFormAdvection(*this);
}