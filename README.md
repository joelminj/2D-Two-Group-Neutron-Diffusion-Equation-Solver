# 2D-Two-Group-Neutron-Diffusion-Equation-Solver

This MATLAB script computes the steady-state neutron flux distribution in a 2D homogeneous rectangular reactor core using Two-Group Neutron Diffusion Theory. The continuous differential equations are discretized using the Finite Difference Method (FDM) and solved as a coupled system of linear algebraic equations.

<img width="1400" height="600" alt="rp_activity_1" src="https://github.com/user-attachments/assets/767a30de-7870-48c7-a12e-34fbbcad682f" />


### 1. Governing Physical Equations

The neutron energy spectrum is divided into two energy groups: Group 1 (Fast Neutrons) and Group 2 (Thermal Neutrons). We assume a fixed, independent source $S_{strength}$ that emits only fast neutrons. The steady-state coupled 2D diffusion equations are:

**Fast Group (Group 1):**


$$-D_1 \nabla^2 \phi_1(x,y) + \Sigma_{R1} \phi_1(x,y) = S_1(x,y)$$

*Where:*

* $D_1$ is the fast diffusion coefficient.
* $\Sigma_{R1} = \Sigma_{a1} + \Sigma_{s1\to2}$ is the macroscopic removal cross-section (representing fast neutrons lost to absorption and scattering down to the thermal group).
* $S_1(x,y)$ is the external fast neutron source.

**Thermal Group (Group 2):**


$$-D_2 \nabla^2 \phi_2(x,y) + \Sigma_{a2} \phi_2(x,y) = \Sigma_{s1\to2} \phi_1(x,y)$$

*Where:*

* $D_2$ is the thermal diffusion coefficient.
* $\Sigma_{a2}$ is the thermal absorption cross-section.
* The term $\Sigma_{s1\to2} \phi_1(x,y)$ acts as the volumetric source for the thermal group, representing fast neutrons that have thermalized.

### 2. Spatial Discretization (Finite Difference Method)

The continuous Laplacian operator $\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$ is approximated on a Cartesian grid using a central **5-point stencil**. For any internal grid node $(i,j)$ with spacing $\Delta x$ and $\Delta y$, the discretized Laplacian is:

$$\nabla^2 \phi_{i,j} \approx \frac{\phi_{i+1, j} - 2\phi_{i,j} + \phi_{i-1, j}}{\Delta x^2} + \frac{\phi_{i, j+1} - 2\phi_{i,j} + \phi_{i, j-1}}{\Delta y^2}$$

Applying this to a generic group diffusion equation ($-D\nabla^2\phi + \Sigma\phi = S$) yields:

$$-D \left[ \frac{\phi_{i+1, j} - 2\phi_{i,j} + \phi_{i-1, j}}{\Delta x^2} + \frac{\phi_{i, j+1} - 2\phi_{i,j} + \phi_{i, j-1}}{\Delta y^2} \right] + \Sigma \phi_{i,j} = S_{i,j}$$

### 3. Boundary Conditions

The code enforces Dirichlet boundary conditions where the flux vanishes at the extrapolated boundaries of the reactor:

$$\phi(0, y) = \phi(L_x, y) = \phi(x, 0) = \phi(x, L_y) = 0$$


To implement this efficiently, the system matrix is formulated *only* for the unknown interior nodes. The boundary values (zeros) are inherently enforced by their absence from the right-hand side of the interior node equations, and are appended purely for post-processing and visualization.

### 4. Matrix Formulation using Kronecker Products

To solve the 2D grid as a 1D vector system ($A\Phi = B$), the 2D indices $(i,j)$ are flattened into a 1D index $k$. To construct the sparse pentadiagonal matrix representing the 2D Laplacian without complex looping, the code utilizes the **Kronecker Tensor Product ($\otimes$)**:

$$L_{2D} = (I_y \otimes D_{xx}) + (D_{yy} \otimes I_x)$$

*Where:*

* $D_{xx}$ and $D_{yy}$ are 1D tridiagonal second-derivative matrices for the $x$ and $y$ directions.
* $I_x$ and $I_y$ are identity matrices of sizes $N_x \times N_x$ and $N_y \times N_y$.

This allows us to concisely define the system matrices for both energy groups:

$$A_1 = D_1 L_{2D} + \Sigma_{R1} I_{2D}$$

$$A_2 = D_2 L_{2D} + \Sigma_{a2} I_{2D}$$

### 5. Volumetric Source Definition

The physical source $S_{strength}$ is defined as a point/planar source at the dead center of the grid $(x_c, y_c)$. To ensure energy conservation within the finite difference volume, the source is divided by the volume of the central node:

$$S_1(x_c, y_c) = \frac{S_{strength}}{\Delta x \Delta y}$$

### 6. Sequential Solution Strategy

Because there is no up-scattering (thermal neutrons do not gain energy to become fast neutrons), the system is lower-triangular and can be solved sequentially without iteration:

1. **Solve Fast Group:**

$$\Phi_1 = A_1^{-1} S_1$$

2. **Calculate Thermal Source:**

$$S_2 = \Sigma_{s1\to2} \Phi_1$$

3. **Solve Thermal Group:**

$$\Phi_2 = A_2^{-1} S_2$$
