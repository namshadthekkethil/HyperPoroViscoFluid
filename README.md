# Unsteady Soft Tissue Mechanics Simulator (FEM-based Hyperelastic–Poroelastic–Viscoelastic Framework)

This repository provides a **C++/libMesh** implementation for simulating **unsteady large-deformation mechanics** in **soft biological tissues**, using the **Finite Element Method (FEM)**.  

The code supports **hyperelastic**, **poroelastic**, and **viscoelastic** material models and is designed for research in **soft tissue biomechanics**, **microvascular flow coupling**, and **cell–tissue interactions** in arbitrary 2D or 3D geometries.

---

## ✨ Features

- **FEM-based unsteady deformation solver** for time-dependent simulations.
- **Constitutive models**:
  - Hyperelastic (Neo-Hookean, Mooney–Rivlin, Ogden)
  - Poroelastic (solid-fluid coupling)
  - Viscoelastic (Maxwell, Kelvin–Voigt)
- **Arbitrary geometries** via libMesh mesh import.
- **Microvascular coupling**:using Darcy flow model.
- **Cell mechanics** 
- **Parallelized** using MPI through libMesh for large-scale FEM simulations.


---

## 🧰 Dependencies

| Dependency | Version | Purpose |
|-------------|----------|----------|
| **C++ Compiler** | C++11+ | Core implementation |
| **libMesh** | ≥ 1.3.0 | FEM framework |
| **PETSc** | ≥ 3.15 | Nonlinear solvers |
| **MPI** | Recommended | Parallel computation |

