#include "gravity.hpp"

void Gravity::calculate_energy_derivatives(const PhysicsState* state, EnergyDerivatives* f) {
    f->force(index)     += gravity_vector.x();
    f->force(index + 1) += gravity_vector.y();
    f->force(index + 2) += gravity_vector.z();
    return;
}
