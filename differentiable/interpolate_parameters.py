import numpy as np
from symulathon import Simulation


def __getPairwiseSpringNodeIndex(spring_node_indices, spring_bend_node_indices):
    """Return a list of all the pairs of node indices which form a spring.

    The ordering of this springs is the same as the one in the parameter vector
    (at least in theory).
    """
    param = spring_node_indices + spring_bend_node_indices
    result = []
    for p1, p2 in zip(param[0::2], param[1::2]):
        result.append((p1, p2))

    return result


sim = Simulation(1, 1)
N, M = sim.getGridDimensions()
spring_node_indices = sim.getSpringIndices()
bspring_node_indices = sim.getBendSpringIndices()
pairwise_node_indices = __getPairwiseSpringNodeIndex(spring_node_indices, bspring_node_indices)
nFlex = len(spring_node_indices) // 2
nBend = len(bspring_node_indices) // 2
nParam = nFlex + nBend


def __getSpringParameterIndex(pairwiseNodeSpringIndex, a: int, b: int):
    """Return the index of the spring which unites nodes a and b.

    pairwiseNodeStringIndex should be created by __getPairwieSpringIndex
    """
    for i, (p1, p2) in enumerate(pairwiseNodeSpringIndex):
        if (p1 == a and p2 == b) or (p1 == b and p2 == a):
            return i

    print("ERROR:GET_SPRING_PARAMETER_INDEX:",
          f"no matching spring with a={a} and b={b}")
    return False


def __find_corners(N: int, M: int):
    """Return 4 node indices of the 4 cloth corners."""
    bottomleft, bottomright, topleft, topright = (0, M*(N-1), M-1, (M-1)*N + N-1)
    return bottomleft, bottomright, topleft, topright


def __getCornerSpringNodes(N: int, M: int):
    """Return the 2 nodes that define the corner tension springs."""
    bl, br, tl, tr = __find_corners(N, M)
    bottomleft = (bl, bl+M)
    bottomright = (br, br-M)
    topleft = (tl, tl+M)
    topright = (tr, tr-M)
    return bottomleft, bottomright, topleft, topright


def __getCornerBendSpringNodes(N: int, M: int):
    """Return the 2 nodes that define the corner bend springs."""
    bl, br, tl, tr = __find_corners(N, M)
    bottomleft = (bl, bl+M+1)
    bottomright = (br-M, br+1)
    topleft = (tl-1, tl+M)
    topright = (tr-M-1, tr)
    return bottomleft, bottomright, topleft, topright


def __interpolate(u: float, v: float, boundary):
    bottomleft, bottomright, topleft, topright = boundary
    w11 = (1-u) * (1-v)
    w12 = (1-u) * v
    w21 = u * (1 - v)
    w22 = u*v
    value = bottomleft * w11 + topleft * w12 +\
        bottomright * w21 + topright * w22
    return value


def generate_parameters(k4: list, kbend4: list):
    """Generate a list of parameters interpolating from the corners."""
    corners = k4, kbend4
    sim = Simulation(1, 1)
    N, M = sim.getGridDimensions()
    mesh = sim.getMesh()
    pos = mesh.getPositions()
    vpos = np.array([pos[i:i+3] for i in range(0, len(pos), 3)])
    x, y, z = vpos.T
    maxX = np.max(x)
    maxZ = np.max(z)
    minX = np.min(x)
    minZ = np.min(z)
    deltaX = maxX - minX
    deltaZ = maxZ - minZ

    k_indices = sim.getSpringIndices()
    nFlex = len(k_indices) // 2
    parameters = [0] * len(pairwise_node_indices)
    for i, (a, b) in enumerate(pairwise_node_indices):
        index = __getSpringParameterIndex(pairwise_node_indices, a, b)
        u1 = x[a] / deltaX
        v1 = z[a] / deltaZ
        u2 = x[b] / deltaX
        v2 = z[b] / deltaZ
        u = (u1 + u2) / 2
        v = (v1 + v2) / 2
        parameters[index] = __interpolate(u, v, corners[i // nFlex])

    return parameters[:nFlex], parameters[nFlex:]


def get_parameter_indices() -> list:
    """Compute indices wrt parameter list of the corner springs."""
    param_indices = [0] * 8
    bottomleft, bottomright, topleft, topright = __getCornerSpringNodes(N, M)
    param_indices[0] = __getSpringParameterIndex(pairwise_node_indices, *bottomleft)
    param_indices[1] = __getSpringParameterIndex(pairwise_node_indices, *bottomright)
    param_indices[2] = __getSpringParameterIndex(pairwise_node_indices, *topleft)
    param_indices[3] = __getSpringParameterIndex(pairwise_node_indices, *topright)
    bottomleft, bottomright, topleft, topright = __getCornerBendSpringNodes(N, M)
    param_indices[4] = __getSpringParameterIndex(pairwise_node_indices, *bottomleft)
    param_indices[5] = __getSpringParameterIndex(pairwise_node_indices, *bottomright)
    param_indices[6] = __getSpringParameterIndex(pairwise_node_indices, *topleft)
    param_indices[7] = __getSpringParameterIndex(pairwise_node_indices, *topright)
    return param_indices
