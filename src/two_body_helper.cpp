#include "two_body_helper.hpp"

namespace two_body {
    void add_force(unsigned int p1, unsigned int p2, const Vec3& force, EnergyDerivatives* f) {
        f->force(p1)     += force.x(); // f0.x
        f->force(p1 + 1) += force.y(); // f0.y
        f->force(p1 + 2) += force.z(); // f0.z

        // 3º Newton law: Equal and oposite reaction
        f->force(p2)     += -force.x(); // f0.x
        f->force(p2 + 1) += -force.y(); // f0.y
        f->force(p2 + 2) += -force.z(); // f0.z
    }

    void add_force_derivatives(unsigned int p1, unsigned int p2, const Mat3& df_dx, const Mat3& df_dv, EnergyDerivatives* f) {
        for (int j=0; j < 3; j++){
            for (int k=0; k < 3; k++){
                // The df_dx derivative
                f->df_dx_triplets.push_back(Triplet(p1 + j, p1 + k, df_dx(j, k)));
                f->df_dx_triplets.push_back(Triplet(p1 + j, p2 + k, -df_dx(j, k)));
                f->df_dx_triplets.push_back(Triplet(p2 + j, p1 + k, -df_dx(j, k)));
                f->df_dx_triplets.push_back(Triplet(p2 + j, p2 + k, df_dx(j, k)));

                // The df_dv derivative
                f->df_dx_triplets.push_back(Triplet(p1 + j, p1 + k, df_dv(j, k)));
                f->df_dx_triplets.push_back(Triplet(p1 + j, p2 + k, -df_dv(j, k)));
                f->df_dx_triplets.push_back(Triplet(p2 + j, p1 + k, -df_dv(j, k)));
                f->df_dx_triplets.push_back(Triplet(p2 + j, p2 + k, df_dv(j, k)));
            }
        }
    }

    void add_force_parameter_derivatives(unsigned int p1, unsigned int p2, const Mat& df_dp, EnergyDerivatives* f, std::vector<unsigned int> parameters) {
        for (unsigned int i = 0; i < 3; i++) {
            for (unsigned int p = 0; p < parameters.size(); p++) {
                int parameter_index = parameters[p];
                f->df_dp(p1 + i, parameter_index) += df_dp(i, p);
                f->df_dp(p2 + i, parameter_index) += -df_dp(i, p);
            }
        }
    }
}
