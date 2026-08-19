! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0
!> Minimal, column-by-column FABM interior-state bridge for MOM6.
module FABM_tracer
use fabm, only : fabm_create_model, type_fabm_model, fabm_standard_variables
use MOM_diag_mediator, only : diag_ctrl
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
  real, pointer :: tr(:,:,:,:) => NULL() !< MOM tracer storage [conc].
  real, pointer :: sources(:,:,:,:) => NULL() !< FABM interior sources [conc s-1].
  real, pointer :: par(:,:,:) => NULL() !< Estimated PAR at layer centers [W m-2].
  real, pointer :: surface_par(:,:) => NULL() !< Surface PAR [W m-2].
  real, pointer :: mask(:,:,:) => NULL() !< FABM wet-cell mask [nondim].
  integer, pointer :: bottom_index(:,:) => NULL() !< Bottom layer index [nondim].
  character(len=200) :: config_file = 'fabm.yaml' !< FABM configuration.
  real :: light_attenuation !< PAR attenuation [m-1].
  real :: fallback_surface_par !< Explicit PAR fallback [native heat-flux units].
  integer :: ntr = 0 !< Number of FABM interior states.
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
  if (.not.associated(CS)) return
  allocate(CS%sources(G%isd:G%ied,G%jsd:G%jed,GV%ke,CS%ntr), source=0.0)
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
  call CS%model%link_horizontal_data(fabm_standard_variables%surface_downwelling_photosynthetic_radiative_flux, &
                                     CS%surface_par)
  call CS%model%start()
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
end subroutine FABM_tracer_column_physics
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
