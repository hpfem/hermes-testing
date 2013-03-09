#include "definitions.h"

/* Weak forms */

CustomWeakFormPoisson::CustomWeakFormPoisson(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                                             std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                                             Hermes::Hermes2DFunction<double>* src_term) : WeakForm<double>(1)
{
  // Jacobian forms.
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, mat_al, lambda_al));
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormDiffusion<double>(0, 0, mat_cu, lambda_cu));

  // Residual forms.
  add_vector_form(new WeakFormsH1::DefaultVectorFormVol<double>(0, Hermes::HERMES_ANY, src_term));
};

TestExactSolution1::TestExactSolution1(MeshSharedPtr mesh)
     : ExactSolutionScalar<double>(mesh) 
{
}

double TestExactSolution1::value(double x, double y) const 
{
  return std::sin(x);
}

void TestExactSolution1::derivatives(double x, double y, double& dx, double& dy) const 
{
  dx = -std::cos(x);
  dy = 0.0;
}

Ord TestExactSolution1::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

TestExactSolution1::~TestExactSolution1() 
{
}

MeshFunction<double>* TestExactSolution1::clone() const
{
	return new TestExactSolution1(this->mesh);
}

TestExactSolution2::TestExactSolution2(MeshSharedPtr mesh)
     : ExactSolutionScalar<double>(mesh) 
{
}

double TestExactSolution2::value(double x, double y) const 
{
  return 2*std::cos(y);
}

void TestExactSolution2::derivatives(double x, double y, double& dx, double& dy) const 
{
  dy = -2*std::sin(y);
  dx = 0.0;
}

Ord TestExactSolution2::ord(Ord x, Ord y) const 
{
  return Ord(10);
}

TestExactSolution2::~TestExactSolution2() 
{
}

MeshFunction<double>* TestExactSolution2::clone() const
{
	return new TestExactSolution2(this->mesh);
}