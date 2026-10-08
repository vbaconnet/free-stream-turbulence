# Free-Stream Turbulence (FST) Plugin for Neko

This plugin implements Free-Stream Turbulence (FST) generation for use as inflow boundary conditions in Neko. The implementation is based on the method described by Schlatter (2001) and extends the original work by Elektra Kluesberg, Prabal Negi, and Philipp Schlatter.

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Compilation](#compilation)
3. [Examples](#Examples)
4. [Mathematical Background](#mathematical-background)
5. [Configuration reference](#configuration-reference)
6. [Operating Modes](#operating-modes)
7. [File Formats](#file-formats)
8. [Troubleshooting](#troubleshooting)
9. [Best Practices](#best-practices)
10. [Limitations](#limitations)

---

## Quick Start

To enable FST in your simulation:

1. **Add FST parameters** to your case file under the `"case.FST"` JSON object
2. **Include the driver module** in your user file
3. **Apply user boundary conditions** on the inflow boundary

**Minimal case file configuration:**
```json
"FST": {
  "enabled": true,
  "t_start": 0.0,
  "t_ramp": 0.001,
  "Uinf": 1.0,
  "Tu": 3.7e-2,
  "L": 11.53e-3,
  "k_start": 63.77,
  "k_end": 1660.0,
  "n_shells": 80,
  "n_pts_per_shell": 40
}
```

<details>
<summary><b>Minimal user file (`user.f90`)</b></summary>

```fortran
module user
  use neko
  use fst_bc_driver
  implicit none
contains
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%initialize => initialize
    u%finalize => finalize
    u%dirichlet_conditions => user_bc
  end subroutine

  subroutine initialize(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, p
    type(coef_t), pointer :: coef
    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")
    p => neko_registry%get_field("p")
    coef => neko_user_access%case%fluid%c_Xh
    call fst_bc_driver_initialize(time%t, u, v, w, p, coef, neko_user_access%case%params)
  end subroutine

  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    if (trim(fields%items(1)%ptr%name) .eq. "u") then
      associate(u => fields%items(1)%ptr, v => fields%items(2)%ptr, w => fields%items(3)%ptr)
        coef => neko_user_access%case%fluid%c_Xh
        call fst_bc_driver_apply(u, v, w, bc, coef, time%t, time%tstep, 0.0_xp, .false.)
      end associate
    end if
  end subroutine

  subroutine finalize(time)
    type(time_state_t), intent(in) :: time
    call fst_bc_driver_finalize()
  end subroutine
end module user
```
</details>

---

## Compilation 

Compilation relies on the neko-generated script `makeneko`. Helper scripts 
are provided to compile on CPU and GPU backends: `makeneko_cpu` 
and `makeneko_gpu` respectively. 

Execute each file to see required arguments and options. 


Below are some quick copy-paste examples:

### GPU, CUDA backend
```bash
./makeneko_gpu user.f90 path/to/neko/install cuda path/to/fst
```

### GPU, HIP backend
```bash
./makeneko_gpu user.f90 path/to/neko/install hip path/to/fst
```

### CPU
```bash
./makeneko_cpu user.f90 path/to/neko/install path/to/fst
```

---

## Examples

Example configurations are provided in the `examples/` directory:

- **`examples/box/`** - Simple box domain with periodic y,z
  - `box.case` - Basic FST generation
  - `box.f90` - User file implementation

- **`examples/channel/`** - Channel flow with non-periodic y
  - `channel.case` - FST with custom fringe bounds
  - `channel.f90` - User file implementation

**To run an example:**
```bash
# Go to the desired example
cd examples/box

# Compile using the desired backend
./makeneko_cpu box.f90 path/to/neko ../../src

# Then run with pre-generated files
./neko box.case
```

---

## Mathematical Background

### Definitions

We define the following variables:

- $U_\infty$, the free-stream velocity.
- $Tu$, the turbulence intensity, defined as $ Tu = u'/U_\infty $
- $q$, the turbulent kinetic energy (TKE), defined as

$$ q = \frac{3}{2}u'^2=\frac{3}{2}(U_\infty Tu)^2 $$

- $L$, the turbulent length scale.
- $\mathbf{k}$, the total wavenumber $\mathbf{k}=(k_x, k_y, k_z)$, with $|\mathbf{k}|=k$.

- the Von Karman spectrum for homogenous, isotropic turbulence (HIT) given as

$$
E(k,L,q) = \frac{2}{3}\frac{1.606(kL)^4}{(1.350+(kL)^2)^{17/6}} Lq.
$$

### Turbulence Generation

The FST method generates synthetic turbulence using a sum of random Fourier modes:

$$\mathbf{u}_{FST}(\mathbf{x}, t) = \sum_{n=1}^{N} \hat{\mathbf{u}}_n A_n\sin(\mathbf{k}_n \cdot [\mathbf{x} - t\mathbf{U}_c]+ \phi_n)$$

Where:
- $\mathbf{U}_c=(U_\infty,0,0)$ the convective velocity
- $\hat{\mathbf{u}}_n$ are divergence-free unit vectors, i.e. $|\hat{\mathbf{u}}_n|=1$ and $\nabla \cdot \hat{\mathbf{u}}_n=0$
- $\phi_n$ are randomly generated phase shifts in $[0, 2\pi]$
- $A_n$ are the amplitudes of each mode, computed as

$$
  \frac{A_n^2}{2}=\frac{2}{M}E(k_n,L,q_0)\Delta k_n.
$$
  
  Here, $M$ is the total number of points per shell. For a given shell, 
  we distribute the amplitude uniformly over all the modes in said shell. The above equation
  can be retrieved based on the equality (definition)

$$
  \int_{0}^{\infty}{E(k,L,q)} = q = \frac{1}{2}\langle u_i^2\rangle
$$

  and inserting the equation for $\mathbf{u}_{FST}$.

- Note the use of $q_0$ instead of $q$ for the generation of amplitudes. This is because of the above equality

$$
  \int_{0}^{\infty}{E(k,L,q)dk} = q.
$$

  Because we discretize the spectrum in the interval $[k_{start}; k_{end}]$,
  if we integrate this discretized spectrum we only get an approximation of 
  the integral:

$$
  \int_{0}^{\infty}{E(k,L,q)dk} \simeq \int_{k_{start}}^{k_{end}}{E(k,L,q)dk} = \sum_{n=1}^{N} E(k_n)\Delta k_n
$$

  and consequently we do not recover the TKE by integrating the spectrum, i.e.

$$
  \sum_{n=1}^{N} E(k_n,L,q)\Delta k_n \ne q.
$$

  To fix this problem we search for a $q_0$ such that 

$$
  \sum_{n=1}^{N} E(k_n,L,q_0)\Delta k_n = q.
$$

  Thankfully, the Von Karman spectra is linear in $q$, therefore we can take 
  replace

$$
  E(k,L,q_0) = q_0 E(k,L,1),
$$

  and find the expression for $q_0$,

$$
 q_0 = \frac{q}{\sum_{n=1}^{N} E(k_n,L,1)\Delta k_n} = \frac{3/2(U_\infty Tu)^2}{\sum_{n=1}^{N} E(k_n,L,1)\Delta k_n}
$$

### Fringe Function

The spatial fringe ensures smooth blending at boundaries:

$$\lambda(y,z) = \lambda_y(y) \cdot \lambda_z(z),\quad \lambda_u = S\left(\frac{u - u_{start}}{\delta_{rise}}\right) - S\left(\frac{u - u_{end}}{\delta_{fall}} + 1\right)$$

With the smooth step function:

$$S(x) = \begin{cases} 0 & x \leq 0 \\ 1 & x \geq 1 \\ \frac{1}{1 + e^{1/(x-1) + 1/x}} & \text{otherwise} \end{cases}$$

Where $\delta_{rise} = \delta_{fall} = \alpha \cdot L_u$ and $L_u$ is the domain length.

### Temporal Ramp

Gradual introduction of turbulence:
$$\text{ramp}(t) = \begin{cases} 0 & t \leq t_{start} \\ \frac{t - t_{start}}{t_{ramp}} & t_{start} < t < t_{ramp} \\ 1 & t \geq t_{ramp} \end{cases}$$

---

## Configuration Reference

All FST parameters are specified in the case file under the `"FST"` JSON object (path: `case.FST.*`).

### Main Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `enabled` | logical | No | `true` | Enable/disable FST generation |
| `read_from_files` | logical | No | `false` | Read FST from pre-generated files instead of generating |
| `files_output_path` | string | No | `"./FST_output_files"` | Directory to write FST files |
| `read_files_path` | string | Yes* | - | Path to directory with pre-generated FST files |
| `Uinf` | real | Yes** | - | Free-stream velocity |
| `Tu` | real | Yes** | - | Turbulence intensity (e.g., 0.037 = 3.7%) |
| `L` | real | Yes** | - | Turbulent integral length scale |
| `k_start` | real | Yes** | - | Start of wavenumber range |
| `k_end` | real | Yes** | - | End of wavenumber range |
| `n_shells` | integer | Yes** | - | Number of spherical shells in wavenumber space |
| `n_pts_per_shell` | integer | Yes** | - | Maximum number of points per shell |
| `seed` | integer | No | `-143` | Random seed. **Must be negative** |
| `t_start` | real | No | `0.0` | Time to start applying FST |
| `t_ramp` | real | **Yes** | - | Duration of linear ramp. **Mandatory** |

*Required if `read_from_files=true`
**Required if `read_from_files=false`

### Fringe Parameters

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `periodic_x` | logical | No | `false` | **NOT SUPPORTED** - will error if true |
| `periodic_y` | logical | No | `false` | Periodicity in y-direction |
| `periodic_z` | logical | No | `false` | Periodicity in z-direction |
| `alpha` | real | No | `0.1` | Fringe width as fraction of domain length |
| `ystart` | real | No | `ymin` | Start of fringe region in y |
| `yend` | real | No | `ymax` | End of fringe region in y |
| `zstart` | real | No | `zmin` | Start of fringe region in z |
| `zend` | real | No | `zmax` | End of fringe region in z |

---

## Operating Modes

### Mode 1: Generate FST on the Fly (Default)

Set `read_from_files: false` (or omit it). The plugin generates all FST parameters at runtime.

**Required parameters:**
```json
"FST": {
  "enabled": true,
  "Uinf": 1.0,
  "Tu": 0.037,
  "L": 0.01153,
  "k_start": 63.77,
  "k_end": 1660.0,
  "n_shells": 80,
  "n_pts_per_shell": 40,
  "t_ramp": 0.001
}
```

**Generated files:** The plugin writes these to `files_output_path`:
- `fst.config` - Complete FST configuration
- `fst_spectrum.csv` - Spectral mode data
- `bb.txt` - Phase shifts
- `sphere.dat` - Shell information

### Mode 2: Read Pre-Generated FST Files

Set `read_from_files: true` and provide the path to existing FST files.

**Required parameters:**
```json
"FST": {
  "enabled": true,
  "read_from_files": true,
  "read_files_path": "./precomputed_fst",
  "t_ramp": 0.001
}
```

**Required files in `read_files_path`:**
- `fst_spectrum.csv` - **Mandatory**
- `bb.txt` - **Mandatory**
- `sphere.dat` - Required if `fst.config` is missing
- `fst.config` - Optional but recommended

**Note:** If the file `fst.config` does not exist, you must specify at least
the free-steam velocity `Uinf`. Be aware that estimated values that will be
logged will be wrong since there will be missing variables.

**Note:** Periodicity flags (`periodic_y`, `periodic_z`) must match those used during file generation, or an error will be raised.

---

## File Formats

### fst.config

This file is new to the Neko implementation of the FST code. Using the 
"on-the-fly" generation (operating mode 1) will print this file in your `files_output_path` directory. 

The file contains a summary of the parameters used for the FSt generation,
in a key-value pairs format:

```
Uinf 1.0
Tu 0.037
L 0.01153
k_start 63.77
k_end 1660.0
n_shells 80
n_max_pts_per_shell 40
n_eff_pts_per_shell 40
n_modes 6400
periodic_x F
periodic_y T
periodic_z T
seed -143
```

### fst_spectrum.csv

Comma-separated values with header:
```
shell,kx,ky,kz,amp,u_hat_pn1,u_hat_pn2,u_hat_pn3
1,12.34,5.67,8.90,0.123,0.456,-0.789,0.321
2,-3.45,6.78,-9.01,0.123,-0.123,0.456,-0.789
...
```

**Columns:**
1. `shell` - Shell number (1 to n_shells)
2. `kx` - Wavenumber in x-direction
3. `ky` - Wavenumber in y-direction
4. `kz` - Wavenumber in z-direction
5. `amp` - Amplitude for this shell
6-8. `u_hat_pn1-3` - Random divergence-free unit vector components

### bb.txt

This file (and its name) is more of a legacy thing from old implementations.
It prints the phase shifts and x-component of the random unit vectors, which
after being made divergence-free come out as `u_hat_pn1-3`.

**Note:** The second column will contain zeroes if running in "reading" mode,
(operating mode 2), simply because those random numbers cannot be accessed anymore.

Two columns per line (phase shift, unused):
```
-0.5235988 0.0
0.8817460 0.0
...
```

---


## Troubleshooting

### Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `Increase minimum total wave number!` | `k_start` too high to fit lowest wavenumber in periodic direction | Ensure `k_start` is slightly higher than $2\pi / L$, where L is the length of the periodic direction. |
| `t_start or t_ramp is invalid!` | Invalid time parameters | Ensure `0 <= t_start < t_ramp` |
| `Seed must be negative!` | seed >= 0 | Use negative seed (e.g., -143) |
| `Periodicity in x is not supported` | periodic_x=true | Set `periodic_x: false` |
| `Uinf must be provided` | Missing Uinf with read_from_files | Provide Uinf or fst.config |
| `Error opening fst_spectrum.csv` | Missing/invalid file | Check read_files_path |
| `fst_spectrum.csv should have 8 columns` | Malformed CSV | Regenerate FST files |
| `Unmatching periodicity` | Flags don't match files | Use same periodicity as generation |

---

## Best Practices

### Parameter Selection

| Parameter | Typical Range | Recommendation |
|-----------|---------------|----------------|
| n_shells | 60-120 | 80 |
| n_pts_per_shell | 20-40 | 40 |
| alpha | 0.05-0.2 | 0.1 recommended |

### Reproducibility

- Using the same seed will generate the exact same FST up to single precision 
digits. This is due to the random number generator which works with integers 
on 32 bits. 
- If you want to reproduce up to double precision digits, use the 
`read_from_files` mode so you can read exactly the same wavenumbers, etc.

---

## Limitations

### Current Limitations

1. **Periodicity in x-direction is NOT supported** - Setting `periodic_x: true` will cause an error
2. **Inlet must be in x-direction** - The (y,z) plane is assumed to be the inlet plane
3. **Single inflow boundary** only - Multiple inflows require separate FST objects
4. **Rotation of velocity components** limited to xy-plane via the `angleXY` parameter in `apply_BC`

---

## References

- Schlatter, P. (2001). *Spectral simulation of turbulent flow in a channel with inlet turbulence*. Masters' thesis, KTH Royal Institute of Technology.
- Original implementation by Elektra Kluesberg, Prabal Negi, and Philipp Schlatter.

---

*Last updated: September 2026*
