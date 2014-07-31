#include "hermes2d.h"

#include "../../../testing-core/testing-core.h"

/* Weak forms */

class CustomWeakFormAcoustics : public WeakForm<::complex>
{
public:
  CustomWeakFormAcoustics(std::string bdy_newton, double rho,
                          double sound_speed, double omega);

  ~CustomWeakFormAcoustics();

  std::vector<Hermes1DFunction<::complex>*> fns_1d;
  std::vector<Hermes2DFunction<::complex>*> fns_2d;
};
