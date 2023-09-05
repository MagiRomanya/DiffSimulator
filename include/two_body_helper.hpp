#ifndef TWO_BODY_HELPER_H_
#define TWO_BODY_HELPER_H_

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include <vector>

namespace two_body {
    void add_force(unsigned int p1, unsigned int p2, const Vec3& force, EnergyDerivatives* f);

    void add_force_derivatives(unsigned int p1, unsigned int p2, const Mat3& df_dx, const Mat3& df_dv, EnergyDerivatives* f);

    void add_force_parameter_derivatives(unsigned int p1, unsigned int p2, const Mat& df_dp, EnergyDerivatives* f, std::vector<unsigned int> parameters);
}

#endif // TWO_BODY_HELPER_H_
