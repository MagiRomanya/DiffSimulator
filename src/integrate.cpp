#include <iostream>
#include <vector>

#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "simulation_parameters.hpp"


void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vector* eq_vec, SparseMatrix* eq_mat);

void integrate(const SimulationParameters& sim_param, PhysicsState* state, const EnergyDerivatives& f) {
    const Scalar h = sim_param.TimeStep;
    const unsigned int nDoF = sim_param.q0.size();
    const unsigned int nParameters = sim_param.p.size();

    // Sparse Matrix creation
    // ----------------------------------------------------------------------------------
    SparseMatrix mass_matrix(nDoF, nDoF);
    mass_matrix.setFromTriplets(sim_param.mass.begin(), sim_param.mass.end());
    SparseMatrix df_dx(nDoF, nDoF), df_dv(nDoF, nDoF);
    df_dx.setFromTriplets(f.df_dx_triplets.begin(), f.df_dx_triplets.end());
    df_dv.setFromTriplets(f.df_dv_triplets.begin(), f.df_dv_triplets.end());
    // ----------------------------------------------------------------------------------

    // Construct the system of equations
    // ----------------------------------------------------------------------------------
    Vector equation_vector = h * (f.force + h * df_dx * state->q_dot);
    SparseMatrix equation_matrix = mass_matrix - h * df_dv - h * h * df_dx;
    handle_frozen_dof(sim_param.frozen_dof, &equation_vector, &equation_matrix);
    // ----------------------------------------------------------------------------------

    // Solving the system of equations
    // ----------------------------------------------------------------------------------
    // Gradient conjugate solving method class
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>> cg;
    cg.compute(equation_matrix);
    const Vector delta_q_dot = cg.solve(equation_vector);
    // ----------------------------------------------------------------------------------

    // Update the state with the result
    // ----------------------------------------------------------------------------------
    state->time += h;
    state->q_dot += delta_q_dot;
    state->q += state->q_dot * h;
}

struct FrozenDoFPredicate {
    FrozenDoFPredicate(const std::vector<unsigned int>& frozen_dof) : frozen_dof(frozen_dof) {}

    bool operator() (const Eigen::Index& row, const Eigen::Index& col, const double& value) const {
        // Keep elements in the diagonal and outside the dof column and row
        for (unsigned int i = 0; i < frozen_dof.size(); i++) {
            unsigned int dof = frozen_dof[i];
            if (row==col) return true;
            else if ((row==dof) or (col==dof)) return false;
        }
        return true;
    }
    private:
        const std::vector<unsigned int>& frozen_dof;
};

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vector* eq_vec, SparseMatrix* eq_mat) {
    // Eliminate non zeros from the rows and columns
    (*eq_mat).prune(FrozenDoFPredicate(frozen_dof));
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*eq_vec)[frozen_dof[i]] = 0.0;
    }
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, SparseMatrix* mat) {
    // Eliminate non zeros from the rows and columns
    (*mat).prune(FrozenDoFPredicate(frozen_dof));
}

void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vector* vec) {
    for (unsigned int i = 0; i < frozen_dof.size(); i++) {
        (*vec)[frozen_dof[i]] = 0.0;
    }
}
