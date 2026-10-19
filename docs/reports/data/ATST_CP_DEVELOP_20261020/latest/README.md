# Latest develop build and CP smoke

- Source: `deepmodeling/abacus-develop` `origin/develop`
  `7be5df975` (`Precommit: Prettify source_relax (#8174)`).
- Build job: SAI `1864773`, `4V100`, 1 GPU, `COMPLETED`, 00:05:01.
- Built version: `ABACUS version v3.11.0-beta11`.
- Binary SHA256: `422636c79be9476958045610202fe4414db635b9b3ca8ca9113ba899231859c6`.
- Build flags: `ENABLE_LCAO=ON`, `ENABLE_ELPA=ON`, `USE_CUDA=ON`,
  `CMAKE_CUDA_ARCHITECTURES=70`, `ENABLE_LIBXC=OFF`.
- Modules: `cmake/3.31.6`, `cuda/12.9.1`,
  `openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto`,
  `elpa/2026.02.001-2603-gnu`, `fftw/3.3.10`.
- Cereal headers: LTS toolchain installation, installed locally as a header-only
  CMake package before the ABACUS build.

Both CP boundaries completed on the Pt/LCAO fixture:

| Case | `nelec` | `efermi` (eV) | `mu_calc` (eV) | residual (eV) | `energy` (eV) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `reference_fcp` | 215.966013526 | 5.034404363 | -4.999997234 | 2.77e-6 | -39571.255883550 |
| `compensated_gate` | 217.0 | 9.325994744 | 15.045970064 | -5.94e-7 | -39581.254287417 |

This is an execution/format check on the latest develop build, not an FCP v2
numerical parity or scientific-convergence claim.
