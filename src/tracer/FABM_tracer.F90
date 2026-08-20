! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!> Minimal, column-by-column FABM interior-state bridge for MOM6.
module FABM_tracer
use fabm, only : fabm_create_model, type_fabm_model, fabm_standard_variables
use MOM_diag_mediator, only : diag_ctrl, register_diag_field, post_data
use MOM_error_handler, only : MOM_error, FATAL
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : MOM_restart_CS
use MOM_sponge, only : sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface, thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type
implicit none ; private
#include <MOM_memory.h>
public register_FABM_tracer, initialize_FABM_tracer, FABM_tracer_set_forcing
public FABM_tracer_column_physics, FABM_tracer_surface_state, FABM_tracer_end
type, public :: FABM_tracer_CS ; private
  class(type_fabm_model), pointer :: model => NULL() !< FABM model instance.
  type(diag_ctrl), pointer :: diag => NULL() !< MOM diagnostic control structure.
  real, pointer :: tr(:,:,:,:) => NULL() !< MOM tracer storage [conc].
  real, pointer :: sources(:,:,:,:) => NULL() !< FABM interior sources [conc s-1].
  real, pointer :: diagnostics(:,:,:,:) => NULL() !< FABM interior diagnostics.
  real, pointer :: temperature(:,:,:) => NULL() !< Potential temperature [degC].
  real, pointer :: salinity(:,:,:) => NULL() !< Practical salinity [ppt].
  real, pointer :: density(:,:,:) => NULL() !< Reference density [kg m-3].
  real, pointer :: pressure(:,:,:) => NULL() !< Hydrostatic pressure [dbar].
  real, pointer :: depth(:,:,:) => NULL() !< Depth at layer centers [m].
  real, pointer :: thickness(:,:,:) => NULL() !< Layer thickness [m].
  real, pointer :: par(:,:,:) => NULL() !< Estimated PAR at layer centers [W m-2].
  real, pointer :: surface_par(:,:) => NULL() !< Surface PAR [W m-2].
  real, pointer :: mask(:,:,:) => NULL() !< FABM wet-cell mask [nondim].
  integer, pointer :: bottom_index(:,:) => NULL() !< Bottom layer index [nondim].
  character(len=200) :: config_file = 'fabm.yaml' !< FABM configuration.
  real :: light_attenuation !< PAR attenuation [m-1].
  real :: fallback_surface_par !< Explicit PAR fallback [native heat-flux units].
  integer :: ntr = 0 !< Number of FABM interior states.
  integer :: nidiag = 0 !< Number of FABM interior diagnostic variables.
  integer, allocatable :: id_source(:) !< MOM diagnostic IDs for FABM state sources.
  integer, allocatable :: id_diagnostic(:) !< MOM diagnostic IDs for FABM diagnostics.
end type FABM_tracer_CS
contains
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
  integer :: n
  if (associated(CS)) call MOM_error(FATAL, 'register_FABM_tracer: associated control structure.')
  allocate(CS)
  call log_version(param_file, mdl, version, '')
  call get_param(param_file, mdl, 'FABM_CONFIG_FILE', CS%config_file, &
                 'YAML configuration file for FABM.', default='fabm.yaml')
  call get_param(param_file, mdl, 'FABM_LIGHT_ATTENUATION', CS%light_attenuation, &
                 'Exponential attenuation used to estimate in-water PAR.', units='m-1', &
                 default=0.04, scale=US%m_to_L**(-1))
  call get_param(param_file, mdl, 'FABM_SURFACE_PAR', CS%fallback_surface_par, &
                 'Constant surface PAR used only when MOM6 supplies no shortwave forcing. '//&
                 'A negative value requires shortwave forcing.', units='W m-2', &
                 default=-1.0, scale=US%W_m2_to_QRZ_T)
  CS%model => fabm_create_model(path=trim(CS%config_file))
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
  character(len=96) :: field_name
  if (.not.associated(CS)) return
  allocate(CS%sources(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%ntr), source=0.0)
  CS%diag => diag
  allocate(CS%id_source(CS%ntr), source=0)
  do n=1,CS%ntr
    field_name = "fabm_"//trim(CS%model%interior_state_variables(n)%name)//"_source"
    CS%id_source(n) = register_diag_field("ocean_model", trim(field_name), diag%axesTL, day, &
      "FABM source tendency: "//trim(CS%model%interior_state_variables(n)%long_name), &
      trim(CS%model%interior_state_variables(n)%units)//" s-1")
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
  allocate(CS%surface_par(G%isd:G%ied,G%jsd:G%jed), source=0.0)
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
  call CS%model%link_interior_data(fabm_standard_variables%downwelling_photosynthetic_radiative_flux, CS%par)
  call CS%model%link_interior_data(fabm_standard_variables%temperature, CS%temperature)
  call CS%model%link_interior_data(fabm_standard_variables%practical_salinity, CS%salinity)
  call CS%model%link_interior_data(fabm_standard_variables%density, CS%density)
  call CS%model%link_interior_data(fabm_standard_variables%pressure, CS%pressure)
  call CS%model%link_interior_data(fabm_standard_variables%depth, CS%depth)
  call CS%model%link_interior_data(fabm_standard_variables%cell_thickness, CS%thickness)
  call CS%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux, &
                                     CS%surface_par)
  if (CS%model%interior_variable_needs_values(CS%model%get_interior_variable_id(fabm_standard_variables%temperature)) .and. .not.associated(tv%T)) &
    call MOM_error(FATAL, "FABM model requires temperature, but this MOM6 configuration has no temperature tracer.")
  if (CS%model%interior_variable_needs_values(CS%model%get_interior_variable_id(fabm_standard_variables%practical_salinity)) .and. .not.associated(tv%S)) &
    call MOM_error(FATAL, "FABM model requires salinity, but this MOM6 configuration has no salinity tracer.")
  call update_FABM_environment(h, G, GV, US, tv, CS)
  call CS%model%start()
  if (CS%nidiag > 0) allocate(CS%diagnostics(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%nidiag), source=0.0)
  if (.not.restart) then
    do n=1,CS%ntr ; do k=1,GV%ke ; do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%tr(i,j,k,n) = CS%model%interior_state_variables(n)%initial_value
    enddo ; enddo ; enddo ; enddo
  endif
end subroutine initialize_FABM_tracer
subroutine FABM_tracer_set_forcing(day_start, G, CS)
  type(time_type), intent(in) :: day_start
  type(ocean_grid_type), intent(in) :: G
  type(FABM_tracer_CS), pointer :: CS
end subroutine FABM_tracer_set_forcing
!> Apply FABM interior sources. Boundary fluxes, attached states and sinking are not yet supported.
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
  real :: depth, par0
  integer :: i, j, k, n
  if (.not.associated(CS)) return
  call update_FABM_environment(h_new, G, GV, US, tv, CS, prediabatic_T, prediabatic_S)
  do j=G%jsd,G%jed ; do i=G%isd,G%ied
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
    CS%surface_par(i,j) = par0 ; depth = 0.0
    do k=1,GV%ke
      depth = depth + (0.5 * h_new(i,j,k) * GV%H_to_MKS)
      CS%par(i,j,k) = par0 * exp(-CS%light_attenuation * depth)
      depth = depth + (0.5 * h_new(i,j,k) * GV%H_to_MKS)
    enddo
  enddo ; enddo
  CS%sources(:,:,:,:) = 0.0
  call CS%model%prepare_inputs()
  do k=1,GV%ke ; do j=G%jsd,G%jed
    call CS%model%get_interior_sources(1, G%ied-G%isd+1, j-G%jsd+1, k, &
                                       CS%sources(G%isd:G%ied,j,k,:))
  enddo ; enddo
  do n=1,CS%ntr ; do k=1,GV%ke ; do j=G%jsc,G%jec ; do i=G%isc,G%iec
    CS%tr(i,j,k,n) = CS%tr(i,j,k,n) + (dt * CS%sources(i,j,k,n))
  enddo ; enddo ; enddo ; enddo
  call CS%model%finalize_outputs()
  call post_FABM_diagnostics(CS)
end subroutine FABM_tracer_column_physics
subroutine post_FABM_diagnostics(CS)
  type(FABM_tracer_CS), pointer :: CS
  real, pointer :: fabm_diagnostic(:,:,:) => NULL()
  integer :: n
  if (.not.associated(CS)) return
  do n=1,CS%ntr
    if (CS%id_source(n) > 0) call post_data(CS%id_source(n), CS%sources(:,:,:,n), CS%diag)
  enddo
  do n=1,CS%nidiag
    if (CS%id_diagnostic(n) <= 0) cycle
    fabm_diagnostic => CS%model%get_interior_diagnostic_data(n)
    if (.not.associated(fabm_diagnostic)) cycle
    CS%diagnostics(:,:,:,n) = fabm_diagnostic
    call post_data(CS%id_diagnostic(n), CS%diagnostics(:,:,:,n), CS%diag)
  enddo
end subroutine post_FABM_diagnostics
!> Update standard FABM environmental fields from MOM6 state.
!! Density is the MOM6 reference density; in-situ density will be added with the EOS stage.
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
end module FABM_tracer
