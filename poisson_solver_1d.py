"""Possion Solver in 1D."""


import numpy as np
from scipy.signal import savgol_filter
from scipy.constants import elementary_charge, epsilon_0


class POISSON1D(object):
    """Create 1D mesh objet."""

    def __init__(self, width, nx):
        self.width = width  # domain width in x direction
        # number of nodes within the boundary/main body
        # total nodes will be (nx+2)
        self.nx = nx
        self.delx = width/(nx+1)  # delta_x
#       <------------------------- width -------------------------->
#       x_0=0.0 ---- x_1 ---- x_2 ---- ... ---- x_n ---- x_n+1=width
#       Metal ---------------- Vacuum ------------------------ Metal
        # init x coordinates
        self.x = np.linspace(0.0, self.width, self.nx+2)

    def calc_invA(self, idiag=0):
        """
        Construct 1d Possion matrix.

             2, -1,
            -1, 2, -1,
        A = [   -1,  2, ...         ]
                     ...
                                  -1
                              -1,  2
        compute invert A.
        """
        self.A = np.zeros((self.nx, self.nx), dtype=np.float)
        self.A[0, 0] = 2.0
        self.A[1, 0] = -1.0
        self.A[-2, -1] = -1.0
        self.A[-1, -1] = 2.0
        for i in range(1, self.nx-1):
            self.A[i-1, i] = -1.0
            self.A[i, i] = 2.0
            self.A[i+1, i] = -1.0
        if idiag == 1:
            print(self.A)
        self.invA = np.linself.Alg.inv(self.A)

    def calc_pot(self, charge, bc, ismooth=0):
        """
        Solve 1d Poisson's equation, return potential.

        Using finite diferences, inverse the matrix
        Poisson's equation in vacuum
        d2/dx2*phi(x) = -rho(x)/eps0
        discretize it with finite differences
        A*P(x)/dx2 = -rho(x)/eps0

             2, -1,
            -1, 2, -1,
        A = [   -1,  2, ...         ]
                     ...
                                  -1
                              -1,  2

        phi(x) = -invA*rho(x)/eps0*dx2 + ax + b
        a = (phi[-1] - phi[0])/width
        b = phi[0] + (invA*rho)[0]*/eps0*dx2
        $param charge: an array of total net charge density, size = nx
        $param bc=(bc_l,bc_r): a 2-tuple, boundary conditions
        """
        self.pot = np.zeros(self.nx+2)
        charge = charge*elementary_charge/epsilon_0
        self.pot[1:-1] = np.matmul(self.invA, charge)*self.delx*self.delx
        self.pot[0], self.pot[-1] = bc[0], bc[-1]
        if ismooth == 1:
            self.pot = savgol_filter(self.pot, 11, 3)


if __name__ == '__main__':
    width, nx = 100.0, 99
