#include "hermes2d.h"

/* Weak forms */

class CustomWeakFormPoisson : public WeakForm<double>
{
public:
  CustomWeakFormPoisson(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                        std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                        Hermes::Hermes2DFunction<double>* src_term);
};
