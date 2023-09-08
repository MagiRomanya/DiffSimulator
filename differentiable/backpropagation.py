from solve_system import solve_system
import meansquareloss


class Backpropagation:
    """
    Basic class which handles backpropagation.

    It stores information during forward
    and calculates the loss gradient in backward.
    """

    def __init__(self, mass: float, dx0dp, dv0dp, DeltaTime: float):
        """Initialize the backpropagation class."""
        # Important quantities
        self.equation_matrix_array = []
        self.dfdp_array = []
        self.dfdx_array = []

        self.mass = mass
        self.h = DeltaTime

        # Initial conditions parameter dependance
        self.dx0dp = dx0dp
        self.dv0dp = dv0dp

        # Loss information
        # self.loss = meansquareloss.MeanSquareLoss()
        self.loss = meansquareloss.MeanSquareLossLastState()
        self.g_array = []
        self.dgdx_array = []
        self.dgdv_array = []

        # Some non essential quantites
        self.x_array = []
        self.v_array = []
        self.f_array = []

    def step(self, x, v, x_t, v_t, equation_matrix, dfdp, dfdx, f=0):
        """Store the necessary information douring forward."""
        self.loss.update_containers(x, v, x_t, v_t)

        self.equation_matrix_array.append(equation_matrix)
        self.dfdp_array.append(dfdp)
        self.dfdx_array.append(dfdx)
        self.g_array.append(self.loss.evaluate())
        self.dgdx_array.append(self.loss.get_position_derivative())
        self.dgdv_array.append(self.loss.get_velocity_derivative())
        # Non essential quantities
        self.x_array.append(x)
        self.v_array.append(v)
        self.f_array.append(f)

    def get_dgdp(self):
        """Run the backpropagation algorithm to calculate the loss gradient."""
        n_states = len(self.dfdp_array)
        n_steps = n_states - 1
        h = self.h
        # self.dgdx and self.dgdv are the partial derivatives
        # dgdx and dgdv are the total derivatives
        dgdx = [a.astype(float) for a in self.dgdx_array]
        dgdv = [a.astype(float) for a in self.dgdv_array]
        S = [0] * n_states

        for i in range(n_steps - 1, -1, -1):
            A = self.equation_matrix_array[i+1]
            eq_vec = (h * dgdx[i+1] + dgdv[i+1])
            sys_solve = solve_system(A, eq_vec)

            dgdx[i] += sys_solve @ (h * self.dfdx_array[i+1]) + dgdx[i+1]
            dgdv[i] += sys_solve @ self.mass
            S[i+1] = sys_solve @ (h*self.dfdp_array[i+1])

        S[0] = dgdx[0] * self.dx0dp + dgdv[0] * self.dv0dp

        return sum(S)

    def get_g(self):
        """Return the loss function value of the recorded simulation."""
        return sum(self.g_array)
