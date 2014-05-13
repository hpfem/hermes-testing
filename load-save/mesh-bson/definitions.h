#include "hermes2d.h"

/* Weak forms */
using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm();

private:
    class CustomFormAdvection : public MatrixFormVol<double>
    {
    public:
      CustomFormAdvection(int i, int j, std::string area);

      virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
        Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

      virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *u, Func<Hermes::Ord> *v,
        GeomVol<Hermes::Ord> *e, Func<Ord> **ext) const;

      virtual MatrixFormVol<double>* clone() const;
    };
};
