## Program Overview

A 2D incompressible Navier-Stokes solver using the Immersed Boundary Method (IBM) with direct forcing. Rigid bodies are represented by Lagrangian markers on an Eulerian fluid grid, with no-slip conditions enforced at the boundary using a fractional-step method on a staggered grid. Arbitrary body shapes are supported with configurable translation, rotation, and oscillation. The solver is CPU-parallelised.

## Quick Start

**Requirements:** `gfortran` (with OpenMP), `ffmpeg`

Execute in terminal (macOS):

```bash
chmod +x run.sh && ./run.sh
```

## Solver Verification

https://github.com/user-attachments/assets/794f3b40-b65e-4d41-9666-d19e67a0cb5b

| Metric             | Current | Reference |
| ------------------ | :-----: | :-------: |
| $St$               |  0.195  |   0.195   |
| $C_{L,\text{avg}}$ |  0.000  |   0.000   |
| $C_{L,\text{min}}$ | -0.781  |  -0.600   |
| $C_{L,\text{max}}$ |  0.781  |   0.600   |
| $C_{L,\text{amp}}$ |  1.562  |   1.200   |
| $C_{D,\text{avg}}$ |  1.585  |   1.500   |
| $C_{D,\text{min}}$ |  1.400  |   1.300   |
| $C_{D,\text{max}}$ |  1.778  |   1.700   |
| $C_{D,\text{amp}}$ |  0.378  |   0.400   |

|          Parameter          | Description               |    Value    |
| :-------------------------: | :------------------------ | :---------: |
|          `nx, ny`           | Grid dimensions           |  700, 400   |
|          `xl, xr`           | Domain x-bounds           | -10.0, 25.0 |
|          `yb, yt`           | Domain y-bounds           | -10.0, 10.0 |
|            `nS`             | Lagrangian markers        |     64      |
|            `Lc`             | Characteristic length     |     1.0     |
|         `Xc0, Yc0`          | Initial body position     |  0.0, 0.0   |
|            `Th0`            | Initial body angle        |      0      |
| `Vx_mean, Vy_mean, Om_mean` | Mean velocities           |   0, 0, 0   |
|        `Ax, Ay, Ath`        | Oscillation amplitudes    |  0, 0.2, 0  |
|      `frx, fry, frth`       | Oscillation frequencies   | 0, 0.195, 0 |
|      `phx, phy, phth`       | Oscillation phases        |   0, 0, 0   |
|           `U_inf`           | Mean freestream velocity  |     1.0     |
|            `Re`             | Reynolds number           |     185     |
|         `wall_type`         | Boundary condition        |    slip     |
|            `dt`             | Time step                 |    0.002    |
|            `nt`             | Total time steps          |   100000    |
|          `n_field`          | Field save interval       |     100     |
|          `n_force`          | Force output interval     |     50      |
|            `tol`            | SOR convergence tolerance |   1.0e-5    |
|          `maxiter`          | Max SOR iterations        |   200000    |
|           `omega`           | SOR relaxation factor     |     1.9     |

## Governing Equations

<img width="920" height="958" alt="Image" src="https://github.com/user-attachments/assets/854a9466-f0ca-4e95-89e7-68fb019782f9" />

## Solution Algorithm

<img width="920" height="1025" alt="Image" src="https://github.com/user-attachments/assets/355d5319-02fc-431e-bdce-e3c240841eca" />

## Appendix

<img width="920" height="914" alt="Image" src="https://github.com/user-attachments/assets/22bbfd47-0b73-45af-ab69-c208c856b351" />

## References

- Chorin (1968). _Numerical solution of the Navier-Stokes equations_.
- Fadlun et al. (2000). _Combined immersed-boundary finite-difference methods for three-dimensional complex flow simulations_.
- Guilmineau & Queutey (2002). _A numerical simulation of vortex shedding from an oscillating circular cylinder_.
- Peskin (2002). _The immersed boundary method_.
- Uhlmann (2005). _An immersed boundary method with direct forcing for the simulation of particulate flows_.
- Taira & Colonius (2007). _The immersed boundary method: a projection approach_.
- Kempe & Fröhlich (2012). _An improved immersed boundary method with direct forcing for the simulation of particle laden flows_.
