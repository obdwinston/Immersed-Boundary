## Program Overview

A 2D incompressible Navier-Stokes solver using the Immersed Boundary Method (IBM) with direct forcing. Rigid bodies are represented by Lagrangian markers on an Eulerian fluid grid, with no-slip conditions enforced at the boundary using a fractional-step method on a staggered grid. Arbitrary body shapes are supported with configurable translation, rotation, and oscillation. The solver is CPU-parallelised.

https://github.com/user-attachments/assets/82c106be-e706-484a-aaa6-1fd105f438df

## Quick Start

**Requirements:** `gfortran` (with OpenMP), `ffmpeg`

Execute in terminal (macOS):

```bash
chmod +x run.sh && ./run.sh
```

## Solver Verification

https://github.com/user-attachments/assets/b32c15c9-b11c-461c-83b2-93ae837137fa

| Metric             | Current | Reference | Error (%) |
| ------------------ | :-----: | :-------: | :-------: |
| $St$               |  0.195  |   0.195   |    0.2    |
| $C_{L,\text{avg}}$ |  0.000  |   0.000   |    0.0    |
| $C_{L,\text{min}}$ | -0.653  |  -0.600   |    8.8    |
| $C_{L,\text{max}}$ |  0.653  |   0.600   |    8.9    |
| $C_{L,\text{amp}}$ |  1.306  |   1.200   |    8.9    |
| $C_{D,\text{avg}}$ |  1.649  |   1.500   |    9.9    |
| $C_{D,\text{min}}$ |  1.455  |   1.300   |   11.9    |
| $C_{D,\text{max}}$ |  1.843  |   1.700   |    8.4    |
| $C_{D,\text{amp}}$ |  0.389  |   0.400   |    2.9    |

|          Parameter          | Description               |    Value    |
| :-------------------------: | :------------------------ | :---------: |
|          `nx, ny`           | Grid dimensions           |  1400, 800  |
|          `xl, xr`           | Domain x-bounds           | -10.0, 25.0 |
|          `yb, yt`           | Domain y-bounds           | -10.0, 10.0 |
|            `nS`             | Lagrangian markers        |     128     |
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
