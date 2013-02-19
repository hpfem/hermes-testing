#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

/* Weak forms */

class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                        std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                        Hermes::Hermes2DFunction<double>* src_term);
};

class TestExactSolution1 : public ExactSolutionScalar<double>
{
public:
  TestExactSolution1(const Mesh* mesh);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~TestExactSolution1();
	
	virtual MeshFunction<double>* clone() const;
};

class TestExactSolution2 : public ExactSolutionScalar<double>
{
public:
  TestExactSolution2(const Mesh* mesh);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~TestExactSolution2();
	
	virtual MeshFunction<double>* clone() const;
};