# FABM support in MOM6

## Status and scope

This branch provides an optional, generic FABM bridge for **interior** FABM
state variables.  FABM state variables declared in `fabm.yaml` are registered
as MOM6 tracers and therefore use MOM6 transport, diagnostics and restarts.
The bridge supplies temperature, practical salinity, reference density,
pressure, depth, layer thickness and photosynthetically active radiation (PAR).
It supports interior sources, surface and bottom fluxes on interior tracers,
and FABM vertical movement.  Surface-attached and bottom-attached FABM state
variables are deliberately rejected until their 2-D storage, diagnostics and
restart support are implemented.

FABM is a compile-time option.  A normal MOM6 executable can be built without
FABM.  Setting `USE_FABM_TRACERS = True` in such an executable produces a clear
fatal error rather than silently omitting biology.

## Reproducible source set

This document describes the source side of the bridge.  Coupled SIS2 builds
also require a compatible `MOM6-examples` checkout; see its `FABM_PORT.md`.
Record the three Git revisions for every experiment:

```text
MOM6:          git rev-parse HEAD
MOM6-examples: git -C ../MOM6-examples rev-parse HEAD
FABM:          git -C ../fabm rev-parse HEAD
```

## Ocean-only build

Install FABM with the same Fortran compiler, MPI implementation and real kind
as MOM6.  Let `FABM_ROOT` be the FABM installation prefix containing `include/`
and `lib/` or `lib64/`.

```bash
module purge
module load intel-oneapi-compilers/2021.4.0
module load intel-oneapi-mpi/2021.4.0
module load netcdf-fortran/4.5.3--intel-oneapi-mpi--2021.4.0--intel--2021.4.0

export FABM_ROOT=/absolute/path/to/fabm-install
cd MOM6
mkdir -p build_fabm
cd build_fabm
CPPFLAGS="-D_FABM_" \
FCFLAGS="-I${FABM_ROOT}/include" \
LDFLAGS="-L${FABM_ROOT}/lib64" \
LIBS="-lfabm" \
../ac/configure
make -j 8
```

Use `lib` instead of `lib64` if that is where the FABM installation places
`libfabm.a`.  Do not commit `build_fabm/`, object files or executables.

## Runtime configuration

Add the following to `MOM_override` or `MOM_input` for a FABM-enabled run:

```text
USE_FABM_TRACERS = True
FABM_CONFIG_FILE = "fabm.yaml"
FABM_LIGHT_ATTENUATION = 0.04
```

For a fresh start from a NetCDF file containing one variable per FABM interior
state (same variable names and units as FABM):

```text
FABM_INITIAL_CONDITIONS_FILE = "INPUT/fabm_npzd_initial.nc"
FABM_INITIAL_CONDITIONS_ON_GRID = True
```

The file is used only for a cold start.  Restart runs read MOM6 restart files,
not the initial-condition file.  If the source file is on another horizontal
grid, set `FABM_INITIAL_CONDITIONS_ON_GRID = False` and provide coordinates
supported by MOM6’s standard z-space tracer initializer.

## Current limitations

- No surface- or bottom-attached FABM states yet.
- Density is currently reference-density based; density-sensitive ecosystem
  configurations need an in-situ EOS connection.
- Scientific initial fields require documented remapping and unit/model
  conversions.  The synthetic NPZD file used in regression testing is not a
  scientific initialization.
