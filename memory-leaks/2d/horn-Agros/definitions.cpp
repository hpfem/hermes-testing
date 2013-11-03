#include "definitions.h"

CustomWeakFormAcoustics::CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                                                 double sound_speed, double omega) : WeakForm<complex>(1) 
{
  std::complex<double>  ii =  std::complex<double>(0.0, 1.0);

  // Jacobian.
  this->fns_1d.push_back(new Hermes1DFunction<complex>(1.0 / rho));
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion<complex>(0, 0, HERMES_ANY, this->fns_1d.back(), HERMES_SYM));
  this->fns_2d.push_back(new Hermes2DFunction<complex>(-sqr(omega) / rho / sqr(sound_speed)));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol<complex>(0, 0, HERMES_ANY, fns_2d.back(), HERMES_SYM));
  this->fns_2d.push_back(new Hermes2DFunction<complex>(-ii * omega / rho / sound_speed));
  add_matrix_form_surf(new WeakFormsH1::DefaultMatrixFormSurf<complex>(0, 0, bdy_newton, fns_2d.back()));

  // Residual.
  this->fns_1d.push_back(new Hermes1DFunction<complex>(1.0 / rho));
  add_vector_form(new WeakFormsH1::DefaultResidualDiffusion<complex>(0, HERMES_ANY, this->fns_1d.back()));
  this->fns_2d.push_back(new Hermes2DFunction<complex>(-sqr(omega) / rho / sqr(sound_speed)));
  add_vector_form(new WeakFormsH1::DefaultResidualVol<complex>(0, HERMES_ANY, fns_2d.back()));
  this->fns_2d.push_back(new Hermes2DFunction<complex>(-ii * omega / rho / sound_speed));
  add_vector_form_surf(new WeakFormsH1::DefaultResidualSurf<complex>(0, bdy_newton, fns_2d.back()));
}

CustomWeakFormAcoustics::~CustomWeakFormAcoustics()
{
  for (int i = 0; i < this->fns_1d.size(); i++)
    delete this->fns_1d[i];

  for (int i = 0; i < this->fns_2d.size(); i++)
    delete this->fns_2d[i];
}