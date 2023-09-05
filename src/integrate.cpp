#include "linear_algebra.hpp"
#include "physics_state.hpp"
#include "simulation_parameters.hpp"
#include <vector>


void handle_frozen_dof(const std::vector<unsigned int>& frozen_dof, Vector* eq_vec, SparseMatrix* eq_mat);

void integrate(const SimulationParameters& sim_param, PhysicsState* state, const EnergyDerivatives& f) {
    const Scalar h = sim_param.TimeStep;
    const unsigned int nDoF = sim_param.q0.size();
    const unsigned int nParameters = sim_param.p.size();

    SparseMatrix mass_matrix(nDoF, nDoF);
    mass_matrix.setFromTriplets(sim_param.mass.begin(), sim_param.mass.end());

    SparseMatrix df_dx(nDoF, nDoF), df_dv(nDoF, nDoF);
    df_dx.setFromTriplets(f.df_dx_triplets.begin(), f.df_dx_triplets.end());
    df_dv.setFromTriplets(f.df_dv_triplets.begin(), f.df_dv_triplets.end());

    Vector equation_vector = h * (f.force + h * df_dx * state->q_dot);
    SparseMatrix equation_matrix = mass_matrix - h * df_dv - h * h * df_dx;

    handle_frozen_dof(sim_param.frozen_dof, &equation_vector, &equation_matrix);

    // Gradient conjugate solving method class
    Eigen::ConjugateGradient<Eigen::SparseMatrix<Scalar>> cg;

    // Solving the system of equations
    cg.compute(equation_matrix);
    const Vector delta_q_dot = cg.solve(equation_vector);

    // Updating the state
    state->time += h;
    state->q_dot += delta_q_dot;
    state->q += state->q_dot;
}

struct FrozenDoFPredicate {
    FrozenDoFPredicate(const std::vector<unsigned int>& frozen_dof) : frozen_dof(frozen_dof) {}

    bool operator() (const Eigen::Index& row, const Eigen::Index& col, const double& value) const {
        // Keep elements in the diagonal and outside the dof column and row
        for (unsigned int i = 0; i < frozen_dof.size(); i++) {
            unsigned int dof = frozen_dof[i];
            if ((row==dof) or (col==dof)) return false;
            else if (row==col) return true;
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
        (*eq_vec)[i] = 0.0;
    }
}
