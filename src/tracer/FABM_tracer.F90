! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!> Minimal, column-by-column FABM interior-state bridge for MOM6.
!!
!! FABM owns the ecosystem formulation and state-variable metadata; MOM6 owns
!! prognostic tracer storage, transport, restarts and diagnostics. This module
!! links both models column by column. It supports FABM interior state variables
!! only; surface- and bottom-attached states need separate 2-D MOM6 storage.
! FABM is optional at compile time. Keep an empty module body in ordinary
! MOM6 builds so this source file need not be removed from a source list.
module FABM_tracer
#ifndef _FABM_
implicit none
#else
use fabm, only : fabm_create_model, type_fabm_model, fabm_standard_variables
use MOM_diag_mediator, only : diag_ctrl, register_diag_field, post_data
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_sponge, only : sponge_CS
use MOM_time_manager, only : time_type, get_date, get_time, set_date, operator(-)
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_initialization_from_Z, only : MOM_initialize_tracer_from_Z
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
use mpp_mod, only : mpp_pe, mpp_root_pe
implicit none ; private
#include <MOM_memory.h>
public register_FABM_tracer, initialize_FABM_tracer, FABM_tracer_set_forcing
public FABM_tracer_column_physics, FABM_tracer_surface_state, FABM_tracer_end
type, public :: FABM_tracer_CS ; private
  class(type_fabm_model), pointer :: model => NULL() !< FABM model instance.
  type(diag_ctrl), pointer :: diag => NULL() !< MOM diagnostic control structure.
  real, pointer :: tr(:,:,:,:) => NULL() !< MOM tracer storage [conc].
  real, pointer :: sources(:,:,:,:) => NULL() !< FABM interior sources [conc s-1].
  real, pointer :: surface_flux(:,:,:) => NULL() !< FABM surface fluxes [conc m s-1].
  real, pointer :: bottom_flux(:,:,:) => NULL() !< FABM bottom fluxes [conc m s-1].
  real, pointer :: bottom_state_sources(:,:,:) => NULL() !< FABM bottom-state sources [state s-1].
  real, pointer :: vertical_velocity(:,:,:,:) => NULL() !< FABM movement velocities [m s-1].
  ! A single local 3-D workspace is reused while diagnostics are posted.
  ! Do not retain one full grid per FABM diagnostic: full BFM configurations
  ! expose many diagnostics, making a four-dimensional cache prohibitively
  ! large on a distributed MOM6 grid.
  real, pointer :: diagnostics(:,:,:) => NULL() !< FABM interior diagnostic workspace [model units].
  real, pointer :: temperature(:,:,:) => NULL() !< Potential temperature [degC].
  real, pointer :: salinity(:,:,:) => NULL() !< Practical salinity [ppt].
  real, pointer :: density(:,:,:) => NULL() !< Reference density [kg m-3].
  real, pointer :: pressure(:,:,:) => NULL() !< Hydrostatic pressure [dbar].
  real, pointer :: depth(:,:,:) => NULL() !< Depth at layer centers [m].
  real, pointer :: thickness(:,:,:) => NULL() !< Layer thickness [m].
  real, pointer :: par(:,:,:) => NULL() !< Estimated PAR at layer centers [W m-2].
  real, pointer :: par_attenuation(:,:,:) => NULL() !< PAR attenuation coefficient [m-1].
  real, pointer :: surface_par(:,:) => NULL() !< Surface PAR [W m-2].
  real, pointer :: surface_shortwave(:,:) => NULL() !< Total downwelling shortwave [W m-2].
  real, pointer :: wind_speed(:,:) => NULL() !< 10-m wind speed [m s-1].
  real, pointer :: latitude(:,:) => NULL() !< Geographic latitude at tracer points [degrees_N].
  real, pointer :: longitude(:,:) => NULL() !< Geographic longitude at tracer points [degrees_E].
  real, pointer :: surface_air_pressure(:,:) => NULL() !< Air pressure supplied to FABM [Pa].
  real, pointer :: atmospheric_co2(:,:) => NULL() !< Atmospheric CO2 mole fraction [1e-6].
  real, pointer :: ice_area_fraction(:,:) => NULL() !< Sea-ice area fraction supplied to FABM [1].
  real, pointer :: day_of_year => NULL() !< Decimal days since 1 January [days].
  real, pointer :: mask(:,:,:) => NULL() !< FABM wet-cell mask [nondim].
  integer, pointer :: bottom_index(:,:) => NULL() !< Bottom layer index [nondim].
  character(len=200) :: config_file = 'fabm.yaml' !< FABM configuration.
  character(len=200) :: initial_conditions_file = '' !< Optional cold-start 3-D FABM state file.
  logical :: initial_conditions_ongrid = .false. !< If true, IC file is already on the MOM6 horizontal grid.
  real :: light_attenuation !< PAR attenuation [m-1].
  real :: fallback_surface_par !< Explicit PAR fallback [native heat-flux units].
  real :: fallback_wind_speed !< Wind fallback used only if MOM6 has no u10_sqr [m s-1].
  real :: fallback_surface_air_pressure !< Constant air-pressure fallback [Pa].
  real :: fallback_atmospheric_co2 !< Constant atmospheric CO2 fallback [1e-6].
  real :: fallback_ice_area_fraction !< Constant sea-ice area-fraction fallback [1].
  real :: minimum_active_thickness !< Smallest layer treated as active by FABM [H ~> m or kg m-2].
  logical :: contracted_layer_bug = .true. !< If true, expose contracted MOM6 layers to FABM.
  integer :: maximum_vertical_substeps = 0 !< Maximum sinking substeps; 0 disables the check [nondim].
  logical :: profile_first_call = .false. !< Print a CPU-time breakdown for one FABM update.
  integer :: ntr = 0 !< Number of FABM interior states.
  integer :: nidiag = 0 !< Number of FABM interior diagnostic variables.
  integer, allocatable :: id_source(:) !< MOM diagnostic IDs for FABM state sources.
  integer, allocatable :: id_surface_flux(:) !< MOM diagnostic IDs for FABM surface fluxes.
  integer, allocatable :: id_bottom_flux(:) !< MOM diagnostic IDs for FABM bottom fluxes.
  integer, allocatable :: id_vertical_velocity(:) !< MOM diagnostic IDs for FABM movement velocities.
  integer, allocatable :: id_diagnostic(:) !< MOM diagnostic IDs for FABM diagnostics.
end type FABM_tracer_CS
contains
!> Create FABM and register every YAML-declared interior state as a MOM6 tracer.
!! MOM6 consequently owns advection, restart and standard tracer diagnostics.
logical function register_FABM_tracer(G, GV, US, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  type(unit_scale_type), intent(in) :: US
  type(param_file_type), intent(in) :: param_file
  type(FABM_tracer_CS), pointer :: CS
  type(tracer_registry_type), pointer :: tr_Reg
  type(MOM_restart_CS), target, intent(inout) :: restart_CS
# include "version_variable.h"
  character(len=40) :: mdl = 'FABM_tracer'
  character(len=48) :: flux_units
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: enable_bugs
  integer :: n
  if (associated(CS)) call MOM_error(FATAL, 'register_FABM_tracer: associated control structure.')
  allocate(CS)
  call log_version(param_file, mdl, version, '')
  call get_param(param_file, mdl, 'FABM_CONFIG_FILE', CS%config_file, &
                 'YAML configuration file for FABM.', default='fabm.yaml')
  call get_param(param_file, mdl, 'FABM_INITIAL_CONDITIONS_FILE', CS%initial_conditions_file, &
                 'Optional NetCDF file used to initialize all FABM interior states on a cold start. '//&
                 'Variables must have the same names and units as the FABM interior state variables. '//&
                 'A restart always takes precedence.', default='')
  call get_param(param_file, mdl, 'FABM_INITIAL_CONDITIONS_ON_GRID', CS%initial_conditions_ongrid, &
                 'If true, FABM_INITIAL_CONDITIONS_FILE is already on the MOM6 horizontal tracer grid. '//&
                 'If false, MOM6 laterally regrids the file.', default=.false.)
  call get_param(param_file, mdl, 'FABM_LIGHT_ATTENUATION', CS%light_attenuation, &
                 'Exponential attenuation used to estimate in-water PAR.', units='m-1', &
                 default=0.04, scale=US%m_to_L**(-1))
  call get_param(param_file, mdl, 'FABM_SURFACE_PAR', CS%fallback_surface_par, &
                 'Constant surface PAR used only when MOM6 supplies no shortwave forcing. '//&
                 'A negative value requires shortwave forcing.', units='W m-2', &
                 default=-1.0, scale=US%W_m2_to_QRZ_T)
  call get_param(param_file, mdl, 'FABM_WIND_SPEED', CS%fallback_wind_speed, &
                 'Constant 10-m wind speed used only when MOM6 has no u10_sqr forcing.', &
                 units='m s-1', default=0.0)
  call get_param(param_file, mdl, 'FABM_SURFACE_AIR_PRESSURE', CS%fallback_surface_air_pressure, &
                 'Constant surface air pressure supplied to FABM. MOM6 does not retain this '//&
                 'meteorological field in its generic ocean forcing type.', units='Pa', default=101325.0)
  call get_param(param_file, mdl, 'FABM_ATMOSPHERIC_CO2', CS%fallback_atmospheric_co2, &
                 'Constant atmospheric CO2 mole fraction supplied when MOM6 has no atm_co2 forcing.', &
                 units='1e-6', default=420.0)
  call get_param(param_file, mdl, 'FABM_ICE_AREA_FRACTION', CS%fallback_ice_area_fraction, &
                 'Constant sea-ice area fraction supplied when no SIS2 ice field is linked.', &
                 units='1', default=0.0)
  call get_param(param_file, mdl, 'FABM_PROFILE_FIRST_CALL', CS%profile_first_call, &
                 'If true, print a CPU-time breakdown for the first FABM tracer update on each PE.', &
                 default=.false.)
  call get_param(param_file, mdl, 'ENABLE_BUGS_BY_DEFAULT', enable_bugs, &
                 default=.true., do_not_log=.true.) ! Logged centrally by MOM.F90.
  call get_param(param_file, mdl, 'FABM_CONTRACTED_LAYER_BUG', CS%contracted_layer_bug, &
                 'If true, recover a bug that exposes MOM6 layers contracted to the minimum thickness '//&
                 'to FABM and its vertical-movement calculation.', default=enable_bugs)
  call get_param(param_file, mdl, 'FABM_MINIMUM_ACTIVE_THICKNESS', CS%minimum_active_thickness, &
                 'Layers at or below this thickness are excluded from FABM when '//&
                 'FABM_CONTRACTED_LAYER_BUG=False.', units='m', default=1.0e-1, scale=GV%m_to_H, &
                 do_not_log=CS%contracted_layer_bug)
  call get_param(param_file, mdl, 'FABM_MAXIMUM_VERTICAL_SUBSTEPS', CS%maximum_vertical_substeps, &
                 'Maximum number of FABM vertical-movement substeps before terminating. '//&
                 'A value of 0 disables this safety check.', units='nondim', default=0)
  if (CS%minimum_active_thickness < 0.0) call MOM_error(FATAL, &
    'FABM_tracer.F90, register_FABM_tracer: FABM_MINIMUM_ACTIVE_THICKNESS must be nonnegative.')
  if (CS%maximum_vertical_substeps < 0) call MOM_error(FATAL, &
    'FABM_tracer.F90, register_FABM_tracer: FABM_MAXIMUM_VERTICAL_SUBSTEPS must be nonnegative.')
  ! Defer FABM initialization so only the root MPI rank writes FABM's optional
  ! coupling/task debug logs. Multiple ranks replacing the same files can fail
  ! on a shared filesystem; this does not change model configuration or science.
  CS%model => fabm_create_model(path=trim(CS%config_file), initialize=.false.)
  if (mpp_pe() /= mpp_root_pe()) CS%model%log = .false.
  call CS%model%initialize()
  ! No NPZD names are hard-coded: another FABM interior-state model can use this bridge.
  CS%ntr = size(CS%model%interior_state_variables)
  if (CS%ntr == 0) call MOM_error(FATAL, 'register_FABM_tracer: no FABM interior states.')
  allocate(CS%tr(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%ntr), source=0.0)
  if (GV%Boussinesq) then ; flux_units = 'conc m3 s-1'
  else ; flux_units = 'conc kg s-1' ; endif
  do n=1,CS%ntr
    tr_ptr => CS%tr(:,:,:,n)
    call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, &
      name=trim(CS%model%interior_state_variables(n)%name), &
      longname=trim(CS%model%interior_state_variables(n)%long_name), &
      units=trim(CS%model%interior_state_variables(n)%units), &
      registry_diags=.true., flux_units=flux_units, restart_CS=restart_CS)
  enddo
  register_FABM_tracer = .true.
end function register_FABM_tracer
!> Allocate FABM work arrays, link MOM6 environmental fields and initialize FABM.
!! FABM receives local 1:nx,1:ny array sections after set_domain; MOM6 retains
!! its local-domain bounds, so the links below are views with explicit sections.
subroutine initialize_FABM_tracer(restart, day, G, GV, US, h, param_file, diag, OBC, CS, sponge_CSp, tv)
  logical, intent(in) :: restart
  type(time_type), target, intent(in) :: day
  type(ocean_grid_type), intent(inout) :: G
  type(verticalGrid_type), intent(in) :: GV
  type(unit_scale_type), intent(in) :: US
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  type(param_file_type), intent(in) :: param_file
  type(diag_ctrl), target, intent(in) :: diag
  type(ocean_OBC_type), pointer :: OBC
  type(FABM_tracer_CS), pointer :: CS
  type(sponge_CS), pointer :: sponge_CSp
  type(thermo_var_ptrs), intent(in) :: tv
  integer :: i, j, k, n
  real, pointer :: tr_ptr(:,:,:) => NULL()
  character(len=96) :: field_name
  if (.not.associated(CS)) return
  allocate(CS%sources(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%ntr), source=0.0)
  allocate(CS%surface_flux(G%isd:G%ied,G%jsd:G%jed,CS%ntr), source=0.0)
  allocate(CS%bottom_flux(G%isd:G%ied,G%jsd:G%jed,CS%ntr), source=0.0)
  ! Zero-length placeholder: attached bottom states are not implemented, but FABM
  ! bottom fluxes on interior tracers can still be queried through this interface.
  allocate(CS%bottom_state_sources(G%isd:G%ied,G%jsd:G%jed,0), source=0.0)
  allocate(CS%vertical_velocity(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%ntr), source=0.0)
  if (size(CS%model%surface_state_variables) > 0 .or. size(CS%model%bottom_state_variables) > 0) &
    call MOM_error(FATAL, "FABM bridge does not yet support surface- or bottom-attached state variables.")
  CS%diag => diag
  allocate(CS%id_source(CS%ntr), source=0)
  allocate(CS%id_surface_flux(CS%ntr), source=0)
  allocate(CS%id_bottom_flux(CS%ntr), source=0)
  allocate(CS%id_vertical_velocity(CS%ntr), source=0)
  do n=1,CS%ntr
    field_name = "fabm_"//trim(CS%model%interior_state_variables(n)%name)//"_source"
    CS%id_source(n) = register_diag_field("ocean_model", trim(field_name), diag%axesTL, day, &
      "FABM source tendency: "//trim(CS%model%interior_state_variables(n)%long_name), &
      trim(CS%model%interior_state_variables(n)%units)//" s-1")
    field_name = "fabm_"//trim(CS%model%interior_state_variables(n)%name)//"_surface_flux"
    CS%id_surface_flux(n) = register_diag_field("ocean_model", trim(field_name), diag%axesT1, day, &
      "FABM surface flux: "//trim(CS%model%interior_state_variables(n)%long_name), &
      trim(CS%model%interior_state_variables(n)%units)//" m s-1")
    field_name = "fabm_"//trim(CS%model%interior_state_variables(n)%name)//"_bottom_flux"
    CS%id_bottom_flux(n) = register_diag_field("ocean_model", trim(field_name), diag%axesT1, day, &
      "FABM bottom flux: "//trim(CS%model%interior_state_variables(n)%long_name), &
      trim(CS%model%interior_state_variables(n)%units)//" m s-1")
    field_name = "fabm_"//trim(CS%model%interior_state_variables(n)%name)//"_vertical_velocity"
    CS%id_vertical_velocity(n) = register_diag_field("ocean_model", trim(field_name), diag%axesTL, day, &
      "FABM vertical velocity (positive upward): "//trim(CS%model%interior_state_variables(n)%long_name), "m s-1")
  enddo
  CS%nidiag = size(CS%model%interior_diagnostic_variables)
  allocate(CS%id_diagnostic(CS%nidiag), source=0)
  do n=1,CS%nidiag
    field_name = "fabm_"//trim(CS%model%interior_diagnostic_variables(n)%name)
    CS%id_diagnostic(n) = register_diag_field("ocean_model", trim(field_name), diag%axesTL, day, &
      "FABM diagnostic: "//trim(CS%model%interior_diagnostic_variables(n)%long_name), &
      trim(CS%model%interior_diagnostic_variables(n)%units))
    CS%model%interior_diagnostic_variables(n)%save = CS%id_diagnostic(n) > 0
  enddo
  allocate(CS%temperature(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%salinity(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%density(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=US%R_to_kg_m3 * GV%Rho0)
  allocate(CS%pressure(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%depth(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%thickness(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%par(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  ! Expose the same background attenuation used by the bridge to FABM light models.
  allocate(CS%par_attenuation(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=CS%light_attenuation)
  allocate(CS%surface_par(G%isd:G%ied,G%jsd:G%jed), source=0.0)
  ! Horizontal FABM dependencies are allocated independently of a particular
  ! ecosystem model. Models that do not request them simply ignore these links.
  allocate(CS%surface_shortwave(G%isd:G%ied,G%jsd:G%jed), source=0.0)
  allocate(CS%wind_speed(G%isd:G%ied,G%jsd:G%jed), source=CS%fallback_wind_speed)
  allocate(CS%latitude(G%isd:G%ied,G%jsd:G%jed), source=0.0)
  allocate(CS%longitude(G%isd:G%ied,G%jsd:G%jed), source=0.0)
  allocate(CS%surface_air_pressure(G%isd:G%ied,G%jsd:G%jed), source=CS%fallback_surface_air_pressure)
  allocate(CS%atmospheric_co2(G%isd:G%ied,G%jsd:G%jed), source=CS%fallback_atmospheric_co2)
  ! The ocean tracer interface has no SIS2 ice-state argument. Keep an
  ! explicit constant fallback so FABM light models can run in ocean-only
  ! configurations; coupled ice-ocean support should replace this link with
  ! the instantaneous SIS2 ice concentration.
  allocate(CS%ice_area_fraction(G%isd:G%ied,G%jsd:G%jed), source=CS%fallback_ice_area_fraction)
  allocate(CS%day_of_year, source=0.0)
  allocate(CS%mask(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  allocate(CS%bottom_index(G%isd:G%ied,G%jsd:G%jed), source=GV%ke)
  do k=1,GV%ke ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
    CS%mask(i,j,k) = G%mask2dT(i,j)
  enddo ; enddo ; enddo
  call CS%model%set_domain(G%ied-G%isd+1, G%jed-G%jsd+1, GV%ke, seconds_per_time_unit=1.0)
  call CS%model%set_mask(CS%mask, G%mask2dT(G%isd:G%ied,G%jsd:G%jed))
  call CS%model%set_bottom_index(CS%bottom_index)
  do n=1,CS%ntr
    call CS%model%link_interior_state_data(n, CS%tr(:,:,:,n))
  enddo
  ! These links are views, not copies; MOM6 refreshes them before each FABM call.
  call CS%model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, CS%par)
  call CS%model%link_interior_data(fabm_standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
                                    CS%par_attenuation)
  call CS%model%link_interior_data(fabm_standard_variables%temperature, CS%temperature)
  call CS%model%link_interior_data(fabm_standard_variables%practical_salinity, CS%salinity)
  call CS%model%link_interior_data(fabm_standard_variables%density, CS%density)
  call CS%model%link_interior_data(fabm_standard_variables%pressure, CS%pressure)
  call CS%model%link_interior_data(fabm_standard_variables%depth, CS%depth)
  call CS%model%link_interior_data(fabm_standard_variables%cell_thickness, CS%thickness)
  call CS%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux, &
                                     CS%surface_par)
  ! MOM6 geographic coordinates are static for a run. Surface forcing fields
  ! are refreshed immediately before each FABM evaluation below.
  CS%latitude(:,:) = G%geoLatT(G%isd:G%ied,G%jsd:G%jed)
  CS%longitude(:,:) = G%geoLonT(G%isd:G%ied,G%jsd:G%jed)
  call CS%model%link_horizontal_data(fabm_standard_variables%latitude, CS%latitude)
  call CS%model%link_horizontal_data(fabm_standard_variables%longitude, CS%longitude)
  call CS%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_shortwave_flux, &
                                     CS%surface_shortwave)
  call CS%model%link_horizontal_data(fabm_standard_variables%wind_speed, CS%wind_speed)
  call CS%model%link_horizontal_data(fabm_standard_variables%surface_air_pressure, CS%surface_air_pressure)
  call CS%model%link_horizontal_data(fabm_standard_variables%mole_fraction_of_carbon_dioxide_in_air, &
                                     CS%atmospheric_co2)
  call CS%model%link_horizontal_data(fabm_standard_variables%ice_area_fraction, CS%ice_area_fraction)
  call CS%model%link_scalar(fabm_standard_variables%number_of_days_since_start_of_the_year, CS%day_of_year)
  if (CS%model%interior_variable_needs_values(CS%model%get_interior_variable_id(fabm_standard_variables%temperature)) .and. .not.associated(tv%T)) &
    call MOM_error(FATAL, "FABM model requires temperature, but this MOM6 configuration has no temperature tracer.")
  if (CS%model%interior_variable_needs_values(CS%model%get_interior_variable_id(fabm_standard_variables%practical_salinity)) .and. .not.associated(tv%S)) &
    call MOM_error(FATAL, "FABM model requires salinity, but this MOM6 configuration has no salinity tracer.")
  call update_FABM_environment(h, G, GV, US, tv, CS)
  call CS%model%start()
  ! FABM may expose hundreds of diagnostics. Allocate one local 3-D workspace
  ! and reuse it in post_FABM_diagnostics, rather than allocating a complete
  ! 3-D grid for every diagnostic at initialization.
  if (CS%nidiag > 0) allocate(CS%diagnostics(G%isd:G%ied,G%jsd:G%jed,GV%ke), source=0.0)
  if (.not.restart) then
    if (len_trim(CS%initial_conditions_file) > 0) then
      ! Reuse MOM6's tested z-space tracer initializer. This supports source
      ! files on a different horizontal grid as well as files already remapped
      ! to the MOM6 tracer grid. The file must use FABM state names and units.
      do n=1,CS%ntr
        tr_ptr => CS%tr(:,:,:,n)
        call MOM_initialize_tracer_from_Z(h, tr_ptr, G, GV, US, param_file, &
          trim(CS%initial_conditions_file), trim(CS%model%interior_state_variables(n)%name), &
          ongrid=CS%initial_conditions_ongrid)
      enddo
    else
      ! Preserve the YAML-declared uniform initialization used by existing cases.
      do n=1,CS%ntr ; do k=1,GV%ke ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
        CS%tr(i,j,k,n) = CS%model%interior_state_variables(n)%initial_value
      enddo ; enddo ; enddo ; enddo
    endif
  endif
end subroutine initialize_FABM_tracer
!> Refresh calendar forcing that is independent of an individual tracer column.
!! Radiation, wind and atmospheric forcing are refreshed from MOM6 fluxes in
!! FABM_tracer_column_physics immediately before each FABM evaluation.
subroutine FABM_tracer_set_forcing(day_start, G, CS)
  type(time_type), intent(in) :: day_start
  type(ocean_grid_type), intent(in) :: G
  type(FABM_tracer_CS), pointer :: CS
  integer :: year, month, day, hour, minute, second, elapsed_seconds, elapsed_days
  if (.not.associated(CS)) return
  ! FABM defines this scalar as decimal days elapsed since 1 January. Linking a
  ! scalar lets seasonal models use MOM6's actual calendar rather than a
  ! hard-coded annual cycle. This hook is called when MOM6 updates forcing.
  call get_date(day_start, year, month, day, hour, minute, second)
  call get_time(day_start - set_date(year, 1, 1, 0, 0, 0), elapsed_seconds, elapsed_days)
  CS%day_of_year = real(elapsed_days) + real(elapsed_seconds) / 86400.0
end subroutine FABM_tracer_set_forcing
!> Apply one FABM update within MOM6 column physics.
!! The sequence is: refresh environment; construct PAR; obtain FABM sources,
!! boundary fluxes and vertical movement; apply MOM6 vertical diffusion/fluxes;
!! add local biological sources; then transport relative vertical movement.
!! Attached surface/bottom states remain unsupported, but fluxes on interior
!! tracers and interior vertical movement are applied here.
subroutine FABM_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, tv, CS, &
  prediabatic_T, prediabatic_S, evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old, h_new, ea, eb
  type(forcing), intent(in) :: fluxes
  real, intent(in) :: dt
  type(unit_scale_type), intent(in) :: US
  type(thermo_var_ptrs), intent(in) :: tv
  type(FABM_tracer_CS), pointer :: CS
  real, dimension(:,:,:), optional, intent(in) :: prediabatic_T, prediabatic_S
  real, optional, intent(in) :: evap_CFL_limit, minimum_forcing_depth
  real :: depth, par0, shortwave0
  real :: t_start, t_environment, t_prepare, t_sources, t_vertical, t_mom6, t_finish
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work
  integer :: i, j, k, n
  logical :: profile
  if (.not.associated(CS)) return
  ! This opt-in, one-call timing is intended for diagnosis of large FABM
  ! configurations. CPU_TIME is per PE; only the root PE prints to avoid
  ! multiplying the log volume by the MPI task count.
  profile = CS%profile_first_call
  if (profile) call cpu_time(t_start)
  ! Use post-entrainment thickness and, in the modern path, pre-diabatic T/S.
  call update_FABM_environment(h_new, G, GV, US, tv, CS, prediabatic_T, prediabatic_S)
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    ! Prefer MOM6 visible direct+diffuse shortwave. If only total shortwave is
    ! available, 0.42 is MOM6s visible fraction; the fallback is for idealized cases.
    if (associated(fluxes%sw_vis_dir) .and. associated(fluxes%sw_vis_dif)) then
      par0 = US%QRZ_T_to_W_m2 * (fluxes%sw_vis_dir(i,j) + fluxes%sw_vis_dif(i,j))
    elseif (associated(fluxes%sw)) then
      ! MOM6 uses this same visible fraction when only total shortwave is available.
      par0 = 0.42 * US%QRZ_T_to_W_m2 * fluxes%sw(i,j)
    elseif (CS%fallback_surface_par >= 0.0) then
      par0 = US%QRZ_T_to_W_m2 * CS%fallback_surface_par
    else
      call MOM_error(FATAL, 'FABM bridge requires visible or total shortwave, or FABM_SURFACE_PAR.')
    endif
    ! FABM models such as BFM may need both total shortwave and PAR. MOM6
    ! provides total SW directly when available. In idealized cases that use
    ! only FABM_SURFACE_PAR, derive a compatible total-SW surrogate using the
    ! same 0.42 visible fraction as MOM6's total-SW fallback above.
    if (associated(fluxes%sw)) then
      shortwave0 = US%QRZ_T_to_W_m2 * fluxes%sw(i,j)
    else
      shortwave0 = par0 / 0.42
    endif
    CS%surface_par(i,j) = par0 ; depth = 0.0
    CS%surface_shortwave(i,j) = shortwave0
    if (associated(fluxes%u10_sqr)) then
      CS%wind_speed(i,j) = sqrt(max(0.0, fluxes%u10_sqr(i,j) * US%L_to_m**2 * US%s_to_T**2))
    else
      CS%wind_speed(i,j) = CS%fallback_wind_speed
    endif
    ! atm_co2 is allocated by MOM6 only in configurations that request it
    ! (e.g., MARBL). A FABM-specific constant remains available otherwise.
    if (associated(fluxes%atm_co2)) then
      CS%atmospheric_co2(i,j) = fluxes%atm_co2(i,j)
    else
      CS%atmospheric_co2(i,j) = CS%fallback_atmospheric_co2
    endif
    CS%surface_air_pressure(i,j) = CS%fallback_surface_air_pressure
    ! Use instantaneous SIS2 coverage when the coupled forcing supplies it.
    ! Ocean-only cases retain FABM_ICE_AREA_FRACTION (normally zero/open water).
    if (associated(fluxes%ice_fraction)) then
      CS%ice_area_fraction(i,j) = min(1.0, max(0.0, fluxes%ice_fraction(i,j)))
    else
      CS%ice_area_fraction(i,j) = CS%fallback_ice_area_fraction
    endif
    ! Beer--Lambert PAR at layer centres: PAR [W m-2], attenuation [m-1].
    do k=1,GV%ke
      depth = depth + (0.5 * h_new(i,j,k) * GV%H_to_MKS)
      CS%par(i,j,k) = par0 * exp(-CS%light_attenuation * depth)
      depth = depth + (0.5 * h_new(i,j,k) * GV%H_to_MKS)
    enddo
  enddo ; enddo
  if (profile) call cpu_time(t_environment)
  CS%sources(:,:,:,:) = 0.0
  CS%surface_flux(:,:,:) = 0.0
  CS%bottom_flux(:,:,:) = 0.0
  CS%vertical_velocity(:,:,:,:) = 0.0
  ! FABM tendencies returned below have units of state concentration per second.
  call CS%model%prepare_inputs()
  if (profile) call cpu_time(t_prepare)
  do j=G%jsd,G%jed
    call CS%model%get_surface_sources(1, G%ied-G%isd+1, j-G%jsd+1, &
                                       CS%surface_flux(G%isd:G%ied,j,:))
    call CS%model%get_bottom_sources(1, G%ied-G%isd+1, j-G%jsd+1, &
                                      CS%bottom_flux(G%isd:G%ied,j,:), &
                                      CS%bottom_state_sources(G%isd:G%ied,j,:))
  enddo
  do k=1,GV%ke ; do j=G%jsd,G%jed
    call CS%model%get_interior_sources(1, G%ied-G%isd+1, j-G%jsd+1, k, &
                                       CS%sources(G%isd:G%ied,j,k,:))
  enddo ; enddo
  call CS%model%finalize_outputs()
  if (profile) call cpu_time(t_sources)
  do k=1,GV%ke ; do j=G%jsd,G%jed
    call CS%model%get_vertical_movement(1, G%ied-G%isd+1, j-G%jsd+1, k, &
                                         CS%vertical_velocity(G%isd:G%ied,j,k,:))
  enddo ; enddo
  if (profile) call cpu_time(t_vertical)
  do n=1,CS%ntr
    ! Apply FABM boundary fluxes through MOM6's standard vertical tracer operator,
    ! retaining MOM6 evaporation and flux conventions.
    do k=1,GV%ke ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      h_work(i,j,k) = h_old(i,j,k)
    enddo ; enddo ; enddo
    if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) &
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,n), dt, fluxes, h_work, &
                                           evap_CFL_limit, minimum_forcing_depth)
    call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,n), G, GV, &
                         sfc_flux=GV%Rho0*CS%surface_flux(:,:,n), &
                         btm_flux=GV%Rho0*CS%bottom_flux(:,:,n))
  enddo
  if (profile) call cpu_time(t_mom6)
  ! dt [s] times source [concentration s-1] gives a concentration increment.
  do n=1,CS%ntr ; do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (CS%contracted_layer_bug .or. (CS%mask(i,j,k) > 0.5)) &
      CS%tr(i,j,k,n) = CS%tr(i,j,k,n) + (dt * CS%sources(i,j,k,n))
  enddo ; enddo ; enddo ; enddo
  call apply_FABM_vertical_movement(h_new, dt, G, GV, CS)
  call post_FABM_diagnostics(CS)
  if (profile) then
    call cpu_time(t_finish)
    if (mpp_pe() == mpp_root_pe()) then
      write(*,'(a)') 'NOTE: FABM_PROFILE first tracer update (CPU seconds on root PE):'
      write(*,'(a,f12.3)') '  environment/PAR: ', t_environment - t_start
      write(*,'(a,f12.3)') '  prepare_inputs:  ', t_prepare - t_environment
      write(*,'(a,f12.3)') '  sources:         ', t_sources - t_prepare
      write(*,'(a,f12.3)') '  vert_movement:   ', t_vertical - t_sources
      write(*,'(a,f12.3)') '  MOM6 application:', t_mom6 - t_vertical
      write(*,'(a,f12.3)') '  final/post diag: ', t_finish - t_mom6
      write(*,'(a,f12.3)') '  total:           ', t_finish - t_start
    endif
    CS%profile_first_call = .false.
  endif
end subroutine FABM_tracer_column_physics

!> Move FABM interior states relative to the water with conservative first-order upwind transport.
!! FABM velocities are positive upward; MOM6 layer index increases downward.
!! The sign is reversed when forming MOM6 interface fluxes. Surface and seabed
!! interface fluxes are zero, conserving material in the water column but not yet
!! representing export to a sediment/benthic FABM state. Subcycling limits CFL.
subroutine apply_FABM_vertical_movement(h, dt, G, GV, CS)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  real, intent(in) :: dt
  type(FABM_tracer_CS), pointer :: CS
  real :: flux(1:GV%ke+1), tr_new(1:GV%ke), w_interface, dt_sub, max_cfl
  integer :: i, j, k, n, nsub, isub, k_bottom
  character(len=256) :: mesg
  if (.not.associated(CS)) return
  max_cfl = 0.0
  do n=1,CS%ntr ; do k=2,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (CS%contracted_layer_bug .or. &
        ((CS%mask(i,j,k-1) > 0.5) .and. (CS%mask(i,j,k) > 0.5))) then
      w_interface = 0.5 * (CS%vertical_velocity(i,j,k-1,n) + CS%vertical_velocity(i,j,k,n))
      max_cfl = max(max_cfl, abs(w_interface) * dt / &
                    max(h(i,j,k-1)*GV%H_to_MKS, h(i,j,k)*GV%H_to_MKS))
    endif
  enddo ; enddo ; enddo ; enddo
  if (CS%maximum_vertical_substeps > 0) then
    if (max_cfl > 0.5*real(CS%maximum_vertical_substeps)) then
      write(mesg,'(a,es12.4,a,i0,a)') 'vertical-movement CFL is ', max_cfl, &
        ', requiring more than ', CS%maximum_vertical_substeps, ' substeps.'
      call MOM_error(FATAL, 'FABM_tracer.F90, apply_FABM_vertical_movement: '//trim(mesg))
    endif
  endif
  nsub = max(1, ceiling(max_cfl / 0.5))
  dt_sub = dt / real(nsub)
  do isub=1,nsub ; do n=1,CS%ntr ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    flux(:) = 0.0
    tr_new(:) = CS%tr(i,j,:,n)
    k_bottom = GV%ke
    if (.not.CS%contracted_layer_bug) k_bottom = CS%bottom_index(i,j)
    do k=2,k_bottom
      if (CS%contracted_layer_bug .or. &
          ((CS%mask(i,j,k-1) > 0.5) .and. (CS%mask(i,j,k) > 0.5))) then
        w_interface = -0.5 * (CS%vertical_velocity(i,j,k-1,n) + CS%vertical_velocity(i,j,k,n))
        if (w_interface >= 0.0) then
          flux(k) = w_interface * CS%tr(i,j,k-1,n)
        else
          flux(k) = w_interface * CS%tr(i,j,k,n)
        endif
      endif
    enddo
    do k=1,k_bottom
      if (CS%contracted_layer_bug .or. (CS%mask(i,j,k) > 0.5)) &
        tr_new(k) = CS%tr(i,j,k,n) + dt_sub * (flux(k)-flux(k+1)) / &
                    (h(i,j,k)*GV%H_to_MKS)
    enddo
    CS%tr(i,j,:,n) = tr_new(:)
  enddo ; enddo ; enddo ; enddo
end subroutine apply_FABM_vertical_movement
subroutine post_FABM_diagnostics(CS)
  type(FABM_tracer_CS), pointer :: CS
  real, pointer :: fabm_diagnostic(:,:,:) => NULL()
  integer :: n
  if (.not.associated(CS)) return
  do n=1,CS%ntr
    if (CS%id_source(n) > 0) call post_data(CS%id_source(n), CS%sources(:,:,:,n), CS%diag)
    if (CS%id_surface_flux(n) > 0) call post_data(CS%id_surface_flux(n), CS%surface_flux(:,:,n), CS%diag)
    if (CS%id_bottom_flux(n) > 0) call post_data(CS%id_bottom_flux(n), CS%bottom_flux(:,:,n), CS%diag)
    if (CS%id_vertical_velocity(n) > 0) &
      call post_data(CS%id_vertical_velocity(n), CS%vertical_velocity(:,:,:,n), CS%diag)
  enddo
  do n=1,CS%nidiag
    if (CS%id_diagnostic(n) <= 0) cycle
    fabm_diagnostic => CS%model%get_interior_diagnostic_data(n)
    if (.not.associated(fabm_diagnostic)) cycle
    CS%diagnostics(:,:,:) = fabm_diagnostic
    call post_data(CS%id_diagnostic(n), CS%diagnostics(:,:,:), CS%diag)
  enddo
end subroutine post_FABM_diagnostics

!> Refresh the FABM mask and bottom index from the evolving MOM6 layer thickness.
!! MOM6 contracts inactive layers to its minimum thickness. Excluding those
!! layers prevents sinking CFLs from growing spuriously as the layer layout evolves.
subroutine update_FABM_active_layers(h, G, GV, CS)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  type(FABM_tracer_CS), pointer :: CS
  real :: thickness_cutoff
  integer :: i, j, k
  thickness_cutoff = CS%minimum_active_thickness + &
                     max(1.0e-6*CS%minimum_active_thickness, GV%H_subroundoff)
  CS%mask(:,:,:) = 0.0
  CS%bottom_index(:,:) = 1
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    if (G%mask2dT(i,j) > 0.5) then
      do k=1,GV%ke
        if (h(i,j,k) > thickness_cutoff) then
          CS%mask(i,j,k) = 1.0
          CS%bottom_index(i,j) = k
        endif
      enddo
    endif
  enddo ; enddo
end subroutine update_FABM_active_layers

!> Update standard FABM environmental fields from MOM6 state.
!! Thickness/depth are in m, pressure is approximated in dbar, and density is
!! currently the MOM6 reference density rather than in-situ EOS density.
subroutine update_FABM_environment(h, G, GV, US, tv, CS, prediabatic_T, prediabatic_S)
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  type(unit_scale_type), intent(in) :: US
  type(thermo_var_ptrs), intent(in) :: tv
  type(FABM_tracer_CS), pointer :: CS
  real, dimension(:,:,:), optional, intent(in) :: prediabatic_T, prediabatic_S
  real :: z
  integer :: i, j, k
  if (.not.associated(CS)) return
  if (.not.CS%contracted_layer_bug) call update_FABM_active_layers(h, G, GV, CS)
  if (present(prediabatic_T)) then
    CS%temperature(:,:,:) = prediabatic_T(G%isd:G%ied,G%jsd:G%jed,:) * US%C_to_degC
  elseif (associated(tv%T)) then
    CS%temperature(:,:,:) = tv%T(G%isd:G%ied,G%jsd:G%jed,:) * US%C_to_degC
  endif
  if (present(prediabatic_S)) then
    CS%salinity(:,:,:) = prediabatic_S(G%isd:G%ied,G%jsd:G%jed,:) * US%S_to_ppt
  elseif (associated(tv%S)) then
    CS%salinity(:,:,:) = tv%S(G%isd:G%ied,G%jsd:G%jed,:) * US%S_to_ppt
  endif
  ! Adequate for NPZD; density-sensitive FABM models need an in-situ EOS link.
  CS%density(:,:,:) = US%R_to_kg_m3 * GV%Rho0
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
    z = 0.0
    do k=1,GV%ke
      CS%thickness(i,j,k) = h(i,j,k) * GV%H_to_MKS
      z = z + 0.5 * CS%thickness(i,j,k)
      CS%depth(i,j,k) = z
      CS%pressure(i,j,k) = z * CS%density(i,j,k) * GV%g_Earth / 1.0e4
      z = z + 0.5 * CS%thickness(i,j,k)
    enddo
  enddo ; enddo
end subroutine update_FABM_environment
!> Placeholder for FABM contributions to the coupled MOM6 surface state.
!! It is intentionally empty: this bridge has no FABM surface-attached state or
!! direct atmosphere--biogeochemistry exchange implementation yet.
subroutine FABM_tracer_surface_state(sfc_state, h, G, GV, US, CS)
  type(surface), intent(inout) :: sfc_state
  type(ocean_grid_type), intent(in) :: G
  type(verticalGrid_type), intent(in) :: GV
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h
  type(unit_scale_type), intent(in) :: US
  type(FABM_tracer_CS), pointer :: CS
end subroutine FABM_tracer_surface_state
subroutine FABM_tracer_end(CS)
  type(FABM_tracer_CS), pointer, intent(inout) :: CS
  if (.not.associated(CS)) return
  if (associated(CS%model)) then ; call CS%model%finalize() ; deallocate(CS%model) ; endif
  deallocate(CS)
end subroutine FABM_tracer_end
#endif
end module FABM_tracer
