#ifndef GRAVITY_H
#define GRAVITY_H

#include "physics_state.hpp"

class Gravity {
    public:
        Gravity(unsigned int index, Vec3 gravity)
                :  index(index), gravity_vector(gravity) { };

        void calculate_energy_derivatives(const PhysicsState* state, EnergyDerivatives* f) ;

    private:
        Vec3 gravity_vector;
        unsigned int index;
};


#endif // GRAVITY_H
