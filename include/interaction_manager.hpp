#ifndef INTERACTION_MANAGER_H_
#define INTERACTION_MANAGER_H_

#include <cstdio>
#include <iostream>
#include <vector>

#include "physics_state.hpp"
#include "spring.hpp"
#include "gravity.hpp"


class InteractionManager {
    public:
        void calculate_energy_derivatives(const PhysicsState* state, EnergyDerivatives* f) {
            // Springs
            for (unsigned int i=0; i < m_springs.size(); i++) {
                m_springs[i].calculate_energy_derivatives(state, f);
            }
            // Gravity
            for (unsigned int i=0; i < m_gravity.size(); i++) {
                m_gravity[i].calculate_energy_derivatives(state, f);
            }

        }

        std::vector<Spring> m_springs;
        std::vector<Gravity> m_gravity;
};

#endif // INTERACTION_MANAGER_H_
