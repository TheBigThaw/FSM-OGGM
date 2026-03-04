!----------------------------------------------------------------------!
! Factorial Snow Model adapted for OGGM mass balance                   !
!                                                                      !
! Richard Essery                                                       !
! School of GeoSciences                                                !
! University of Edinburgh                                              !
!----------------------------------------------------------------------!
program FSM

implicit none

! Model layers
integer, parameter :: &
  Nsmx = 3,          &! Maximum number of snow layers
  Nice = 10           ! Number of ice layers
real :: &
  Dmin(Nsmx),        &! Minimum snow layer thicknesses (m)
  Dice(Nice)          ! Ice layer thicknesses (m)
data Dmin / 0.1, 0.2, 0.4 /
data Dice / 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0 /

! Meteorological variables at reference elevation and in bands
integer :: &
  year,              &! Year
  month,             &! Month of year
  day,               &! Day of month
  hour                ! Hour of day
real :: &
  LW,LWz,            &! Incoming longwave radiation (W/m^2)
  Ps,Psz,            &! Surface pressure (Pa)
  Qa,Qaz,            &! Specific humidity (kg/kg)
  Rf,Rfz,            &! Rainfall rate (kg/m^2/s)
  Sf,Sfz,            &! Snowfall rate (kg/m^2/s)
  SW,SWz,            &! Shortwave radiation (W/m^2)
  Ta,Taz,            &! Air temperature (K)
  Ua,Uaz              ! Wind speed (m/s)
real, allocatable :: &  
  dz(:)               ! Elevation bands relative to reference height (m)

! Model state variables  
integer, allocatable :: &
  Nsnw(:)             ! Number of snow layers
real, allocatable :: &
  albs(:),           &! Snow albedo
  Tsrf(:),           &! Snow/ground surface temperature (K)
  Dsnw(:,:),         &! Snow layer thicknesses (m)
  Sice(:,:),         &! Ice content of snow layers (kg/m^2)
  Sliq(:,:),         &! Liquid content of snow layers (kg/m^2)
  Tice(:,:),         &! Ice layer temperatures (K)
  Tsnw(:,:)           ! Snow layer temperatures (K)
logical :: dump       ! True if dump file exists

! Diagnostics
real, allocatable :: &
  Mice(:),           &! Ice mass change (kg/m^2)
  RofI(:),           &! Ice Runoff (kg/m^2)
  RofS(:),           &! Snow Runoff (kg/m^2)
  snd(:),            &! Snow depth (m)
  SWE(:),            &! Snow water equivalent (kg/m^2) 
  SWE0(:),           &! Initial snow water equivalent (kg/m^2) 
  Mice_tot(:)         ! Cumulative ice mass change (kg/m^2)

! Counters
integer :: &
  k,                 &! Elevation band counter
  Nbnd                ! Number of elevation bands
  
call SET_PARAMETERS

! Elevation bands
open(9,file='FSM_bands')
read(9,*) Nbnd
allocate(dz(Nbnd))
read(9,*) dz
close(9)

! Allocate and initialize state variables
allocate(albs(Nbnd))
allocate(Nsnw(Nbnd))
allocate(Tsrf(Nbnd))
allocate(Dsnw(Nsmx,Nbnd))
allocate(Sice(Nsmx,Nbnd))
allocate(Sliq(Nsmx,Nbnd))
allocate(Tice(Nice,Nbnd))
allocate(Tsnw(Nsmx,Nbnd))
albs(:) = 0.8
Nsnw(:) = 0
Tsrf(:) = 273.15
Dsnw(:,:) = 0
Sice(:,:) = 0
Sliq(:,:) = 0
Tice(:,:) = 273.15
Tsnw(:,:) = 273.15

! Initialize state variables from FSM_dump if it exists
inquire(file='FSM_dump',exist=dump)
if (dump) then
  open(9,file='FSM_dump') 
  read(9,*) albs
  read(9,*) Dsnw
  read(9,*) Nsnw
  read(9,*) Sice
  read(9,*) Sliq
  read(9,*) Tice
  read(9,*) Tsnw
  read(9,*) Tsrf
end if
close(9)

! Allocate and initialize diagnostics
allocate(Mice(Nbnd))
allocate(Mice_tot(Nbnd))
allocate(RofI(Nbnd))
allocate(RofS(Nbnd))
allocate(snd(Nbnd))
allocate(SWE(Nbnd))
allocate(SWE0(Nbnd))
Mice_tot(:) = 0
SWE0(:) = 0

! Run the model with meteorlogical data from FSM_met
open(9,file='FSM_met')
open(10,file='FSM_out')
do
  read(9,*,end=1) year,month,day,hour,SW,LW,Rf,Sf,Ta,Qa,Ua,Ps
  do k = 1, Nbnd
    call DOWNSCALE(LW,Ps,Qa,Rf,Sf,SW,Ta,Ua,                            &
                   dz(k),LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz)
    call FSM_TIMESTEP(Nice,Nsmx,                                       &
                      Dice,Dmin,LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz,       &
                      albs(k),Dsnw(:,k),Nsnw(k),Sice(:,k),Sliq(:,k),   &
                      Tice(:,k),Tsnw(:,k),Tsrf(k),                     &
                      Mice(k),RofI(k),RofS(k),snd(k),SWE(k))
  end do
  Mice_tot = Mice_tot + Mice
  if ((month==12) .and. (day==31) .and. (hour==23)) then
    write(10,100) year,SWE(:) - SWE0(:) - Mice_tot(:)
    Mice_tot(:) = 0
    do k = 1, Nbnd
      SWE0(k) = sum(Sice(:,k)) + sum(Sliq(:,k))
    end do
  endif 
end do
1 continue
close(9)

close(10)
100 format(i4,*(f12.2))

! Write state variables at end of run to FSM_dump
open(9,file='FSM_dump')
write(9,*) albs
write(9,*) Dsnw
write(9,*) Nsnw
write(9,*) Sice
write(9,*) Sliq
write(9,*) Tice
write(9,*) Tsnw
write(9,*) Tsrf
close(9)

end program FSM

!-----------------------------------------------------------------------
! FSM parameters -- now declared before fsmpy
!-----------------------------------------------------------------------
module PARAMETERS

! Snow parameters
real :: &
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  bstb,              &! Stability slope parameter
  bthr,              &! Snow thermal conductivity exponent
  hfsn,              &! Snow cover fraction depth scale (m)
  rhof,              &! Fresh snow density (kg/m^3)
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  tcld,              &! Cold snow albedo decay timescale (h)
  tmlt,              &! Melting snow albedo decay timescale (h)
  trho,              &! Snow compaction time scale (h)
  Wirr,              &! Irreducible liquid water content of snow
  z0sn                ! Snow surface roughness length (m)

! Ice parameters
real :: &
  aice,              &! Ice albedo
  z0ic,              & ! Ice surface roughness length (m)
  rho_ice              ! ice density
  
! Metorology downscaling parameters
real :: &
  elapse,            &! Vapour pressure lapse rate (1/m)
  Plapse,            &! Precipitation adjustment factor (1/m)
  Tlapse,            &! Temperature laspe rate (K/m)
  Pf                  ! Precipitation multiplier

integer:: &
  sigmoidDscale       ! sigmoid fn for solid fraction (1 or 0)


  
end module PARAMETERS


!-----------------------------------------------------------------------
! Landing routine for calling FSM from Python
!-----------------------------------------------------------------------

! DNG: some arguments added: Nroff (number of runoff records to return),
!      as well as RoffGl and RoffSn output arrays. Declarations below
!      modified similarly

subroutine FSMpy(Nbnd,Nice,Nsmx,Ntim,Nseg,Nroff,                       &
                 Dice,Dmin,dz,months,LW,Ps,Qa,Rf,Sf,SW,Ta,Ua,          &
                 areas, heights,topo,                                  &
                 albs,Dsnw,Nsnw,Sice,Sliq,Tice,Tsnw,Tsrf,massb,        &
                 RoffGl,RoffSn)

use PARAMETERS, only: &
  rho_ice              ! ice density (kg/m3)

implicit none
integer, intent(in) :: Nbnd,Nice,Nsmx,Ntim,Nseg,Nroff
                                                     ! Nroff: number of runoff records
                                                     ! should be ntim / (# of time steps in roff freq)
real, dimension(Nice), intent(in) :: Dice
real, dimension(Nsmx), intent(in) :: Dmin
real, dimension(Nbnd), intent(in) :: dz
integer, dimension(Ntim), intent(in) :: months
real, dimension(Ntim), intent(in) :: LW,Ps,Qa,Rf,Sf,SW,Ta,Ua
real, dimension(Nseg), intent(in) :: areas, heights  ! represent ice-covered area (m^2) of segments
real, dimension(Nseg), intent(in) :: topo            ! and surface height and bed elev (m) of each segment
                                                     ! (even non-ice covered)
real, dimension(Nbnd), intent(inout) :: albs,Tsrf
integer, dimension(Nbnd), intent(inout) :: Nsnw
real, dimension(Nsmx,Nbnd), intent(inout) :: Dsnw,Sice,Sliq,Tsnw
real, dimension(Nice,Nbnd), intent(inout) :: Tice
real, dimension(12,Nbnd), intent(out) :: massb
real, dimension(Nroff), intent(out) :: RoffGl
real, dimension(Nroff), intent(out) :: RoffSn
integer :: k,n,n_roff
real, dimension(Nbnd) :: Mice,RofI,RofS,snd,SWE,SWE0,Hice
real :: LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz
integer :: recalc_swe

call SET_PARAMETERS

! DNG wasnt sure if this needed
massb(:,:) = 0

! DNG zeroing new arrays
RoffGl(:) = 0
RofI(:) = 0
RofS(:) = 0
n_roff = Ntim / Nroff
do k = 1, Nbnd
  ! DNG we calculate and keep track of ice thickness
  !  solely for the purpose of truncating runoff. We 
  !  do NOT truncate SMB -- this will be important in
  !  melting ice that flows from the terminus to prevent
  !  spurious advance
  Hice(k) = MAX(heights(k)-topo(k),0.0)
end do
do n = 1, Ntim
  recalc_swe = 0
  if (n.eq.1) recalc_swe = 1
  if (n.gt.1) then
          if (months(n).ne.months(n-1)) recalc_swe = 1
  endif
  do k = 1, Nbnd
    if (recalc_swe == 1) then
      SWE0(k) = sum(Sice(:,k)) + sum(Sliq(:,k))
    endif
    call DOWNSCALE(LW(n),Ps(n),Qa(n),Rf(n),Sf(n),SW(n),Ta(n),Ua(n),    &
                   dz(k),LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz) 
    call FSM_TIMESTEP(Nice,Nsmx,                                       &
                      Dice,Dmin,LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz,       &
                      albs(k),Dsnw(:,k),Nsnw(k),Sice(:,k),Sliq(:,k),   &
                      Tice(:,k),Tsnw(:,k),Tsrf(k),                     &
                      Mice(k),RofI(k),RofS(k), snd(k),SWE(k))
! DNG my understanding is that the Runoff found is the amount per time step
!     so here accumulate in time and by band
    RoffSn(1+(n-1)/n_roff) = RoffSn(1+(n-1)/n_roff) + RofS(k) * areas(k)
    ! DNG stop accumulating runoff where thickness is zero
    if (Hice(k) .gt. 0) then
      RoffGl(1+(n-1)/n_roff) = RoffGl(1+(n-1)/n_roff) + RofI(k) * areas(k)
    endif
  end do
  massb(months(n),:) = massb(months(n),:) - Mice
  if (n.eq.Ntim) then
      massb(months(n),:) = massb(months(n),:) + SWE - SWE0
  else if (months(n+1).ne.months(n)) then
      massb(months(n),:) = massb(months(n),:) + SWE - SWE0
  endif
  do k = 1, Nbnd
    Hice(k) = MAX(Hice(k) - Mice(k)/rho_ice,0.0)
  end do
end do

! -- below is to test that arrays are being passed and is commented
! -- floats are only correct to about 1 in 10^-8
!print *, Nseg
!print *, areas(1)
!print *, heights(1)
end subroutine FSMpy

!-----------------------------------------------------------------------
! Physical constants
!-----------------------------------------------------------------------
module CONSTANTS
real, parameter :: &
  cp = 1005,         &! Specific heat capacity of dry air (J/K/kg)
  eps = 0.622,       &! Ratio of molecular weights of water and dry air
  e0 = 611.213,      &! Saturation vapour pressure at Tm (Pa)
  g = 9.81,          &! Acceleration due to gravity (m/s^2)
  hcap_ice = 2100.,  &! Specific heat capacity of ice (J/K/kg)
  hcap_wat = 4180.,  &! Specific heat capacity of water (J/K/kg)
  hcon_ice = 2.24,   &! Thermal conductivity of ice (W/m/K)
  Lc = 2.501e6,      &! Latent heat of condensation (J/kg)
  Lf = 0.334e6,      &! Latent heat of fusion (J/kg)
  Ls = Lc + Lf,      &! Latent heat of sublimation (J/kg)
  Rgas = 287,        &! Gas constant for dry air (J/K/kg)
  Rwat = 462,        &! Gas constant for water vapour (J/K/kg)
  rho_wat = 1000.,   &! Density of water (kg/m^3)
  sb = 5.67e-8,      &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm = 273.15,       &! Melting point (K)
  vkman = 0.4         ! Von Karman constant
end module CONSTANTS

!-----------------------------------------------------------------------
! Driving data characteristics
!-----------------------------------------------------------------------
module DRIVING
real, parameter :: &
  dt = 3600,         &! Timestep (s)
  zT = 2,            &! Temperature and humidity measurement height (m)
  zU = 10             ! Wind speed measurement height (m)
end module DRIVING

!-----------------------------------------------------------------------
! set FSM parameters
!-----------------------------------------------------------------------
subroutine SET_PARAMETERS

use PARAMETERS

implicit none

namelist /params/ asmx,asmn,bstb,bthr,hfsn,rhof,rcld,rmlt,Salb,tcld,   &
                  tmlt,trho,Wirr,z0sn,                                 &
                  aice,z0ic, rho_ice,                                  &
                  elapse,Plapse,Tlapse, Pf, sigmoidDscale

! Snow parameters
  asmx = 0.85         ! Maximum albedo for fresh snow
  asmn = 0.6          ! Minimum albedo for melting snow
  bstb = 5            ! Stability slope parameter
  bthr = 2            ! Snow thermal conductivity exponent
  hfsn = 0.1          ! Snow cover fraction depth scale (m)
  rhof = 100          ! Fresh snow density (kg/m^3)
  rcld = 300          ! Maximum density for cold snow (kg/m^3)
  rmlt = 500          ! Maximum density for melting snow (kg/m^3)
  Salb = 10           ! Snowfall to refresh albedo (kg/m^2)
  tcld = 1000         ! Cold snow albedo decay timescale (h)
  tmlt = 100          ! Melting snow albedo decay timescale (h)
  trho = 200          ! Snow compaction time scale (h)
  Wirr = 0.03         ! Irreducible liquid water content of snow
  z0sn = 0.001        ! Snow surface roughness length (m)

! Ice parameters
  aice = 0.5          ! Ice albedo
  z0ic = 0.01         ! Ice surface roughness length (m)
  rho_ice = 917.0     ! Ice density (kg/m^3)
  
! Metorology downscaling parameters
  elapse = 0.         ! Vapour pressure lapse rate (1/m)
  Plapse = 0.0        ! Precipitation adjustment factor (1/m)
  Tlapse = 6.5e-3     ! Temperature laspe rate (K/m)
  Pf = 1.2            ! Precipitation multiplier
  sigmoidDscale = 0   ! sigmoid fn for solid fraction (1 or 0)
  
open(8,file='nlst') 
read(8,params)
close(8)

end subroutine SET_PARAMETERS

!-----------------------------------------------------------------------
! Downscale reference meteorological data to elevation band
!-----------------------------------------------------------------------
subroutine DOWNSCALE(LW,Ps,Qa,Rf,Sf,SW,Ta,Ua,                          &
                     dz,LWz,Psz,Qaz,Rfz,Sfz,SWz,Taz,Uaz)

use CONSTANTS, only: &
  Tm                  ! Melting point (K)
  
use PARAMETERS, only: &
  elapse,            &! Vapour pressure laps rate (1/m)
  Plapse,            &! Precipitation adjustment factor (1/m)
  Tlapse,            &! Temperature lapse rate (K/m)
  Pf,                &! precip fac
  sigmoidDscale       ! sigmoid fn for solid fraction (1 or 0)

implicit none

! Meteorological variables at reference elevation
real, intent(in) :: &
  dz,                &! Elevation relative to reference height (m)
  LW,                &! Incoming longwave radiation (W/m^2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  SW,                &! Shortwave radiation (W/m^2)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)

! Meteorological variables in elevation band
real, intent(out) :: &
  LWz,               &! Incoming longwave radiation (W/m^2)
  Psz,               &! Surface pressure (Pa)
  Qaz,               &! Specific humidity (kg/kg)
  Rfz,               &! Rainfall rate (kg/m^2/s)
  Sfz,               &! Snowfall rate (kg/m^2/s)
  SWz,               &! Shortwave radiation (W/m^2)
  Taz,               &! Air temperature (K)
  Uaz                 ! Wind speed (m/s)
  
real :: &
  fs,fs0,fs1,        &! Snow fraction
  Pr,                &! Precipitation rate (kg/m2/s)
  Qs                  ! Saturation specific humidity

! Radiation fluxes and wind speed copied without downscaling       
LWz = LW
SWz = SW
Uaz = Ua

! Downscaling to elevation band  
Psz = Ps*exp(-dz/8000)
Taz = Ta - Tlapse*dz
call QSAT(Psz,Taz,Qs)
Qaz = min(Qa,Qs)
Pr = Rf + Sf
Pr = Pr*Pf*(1 + Plapse*dz)/(1 - Plapse*dz)
fs1 = 1 / (1 + exp((Ta - Tm - 1.6)/1))
fs0 = 1 - 0.5*(Taz - 273.15)
fs0 = min(1.,max(fs0,0.))
fs = sigmoidDscale * fs1 + (1-sigmoidDscale) * fs0
Rfz = (1 - fs)*Pr
Sfz = fs*Pr

end subroutine DOWNSCALE

!-----------------------------------------------------------------------
! Call physics subroutines
!-----------------------------------------------------------------------
subroutine FSM_TIMESTEP(Nice,Nsmx,Dice,Dmin,LW,Ps,Qa,Rf,Sf,SW,Ta,Ua,   &
                        albs,Dsnw,Nsnw,Sice,Sliq,Tice,Tsnw,Tsrf,       &
                        Mice,RofI,RofS,snd,SWE                         )

implicit none

! Layers
integer, intent(in) :: &
  Nice,             &! Number of soil layers
  Nsmx               ! Maximum number of snow layers
real :: &
  Dice(Nice),       &! Ice layer thicknesses (m)
  Dmin(Nsmx)         ! Minimum snow layer thicknesses (m)

! Meteoroeogical variables
real, intent(in) :: &
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Rf,                &! Rainfall rate (kg/m2/s)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW,                &! Shortwave radiation (W/m^2)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)
  
! State variables
integer, intent(inout) :: &
  Nsnw                ! Number of snow layers
real, intent(inout) :: &
  albs,              &! Snow albedo
  Tsrf,              &! Surface temperature (K)
  Dsnw(Nsmx),        &! Snow layer thicknesses (m)
  Sice(Nsmx),        &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmx),        &! Liquid content of snow layers (kg/m^2)
  Tice(Nice),        &! Ice layer temperatures (K)
  Tsnw(Nsmx)          ! Snow layer temperatures (K)
  
! Diagnostics  
real, intent(out) :: &
  Mice,              &! Ice mass change (kg/m^2)
  RofI,              &! Cumulative runoff from Ice (kg/m^2)
  RofS,              &! Cumulative runoff from Snow (kg/m^2)
  snd,               &! Snow depth (m)
  SWE                 ! Snow water equivalent (kg/m^2) 

! Fluxes
real :: &
  Esrf,              &! Sublimation rate (kg/m^2/s)
  Gice,              &! Heat flux into ice (W/m^2)
  Gsrf,              &! Heat flux into surface (W/m^2)
  Melt                ! Melt rate (kg/m^2/s)

! DNG in SNOW() and ICE() definitions we separate the snow runoff RofS 
! and ice runoff RofI by call.

call SURFACE(Nice,Nsmx,Dice,Dsnw,LW,Ps,Qa,Sf,Sice,Sliq,SW,Ta,Tice,     &
             Tsnw,Ua,albs,Tsrf,Esrf,Gsrf,Melt)

call SNOW(Nice,Nsmx,Dice,Dmin,Gsrf,Rf,Sf,Ta,Tice,Dsnw,Esrf,Melt,Nsnw,  &
          Sice,Sliq,Tsnw,Gice,RofS,snd,SWE)

call ICE(Nice,Dice,Esrf,Gice,Melt,RofI,Tice,Mice)

end subroutine FSM_TIMESTEP

!-----------------------------------------------------------------------
! Update ice temperatures
!-----------------------------------------------------------------------
subroutine ICE(Nice,Dice,Esrf,Gice,Melt,RofI,Tice,Mice)

! DNG everywhere in ICE() replace Roff by RofI

use CONSTANTS, only : &
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcon_ice,          &! Thermal conductivity of ice (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  Tm                  ! Melting point (K)

use PARAMETERS, only : & 
  rho_ice            ! Density of ice (kg/m^3)

use DRIVING, only : &
  dt                  ! Timestep (s)

implicit none

integer, intent(in) :: &
  Nice                ! Number of ice layers

real, intent(in) :: &
  Dice(Nice),        &! Ice layer thicknesses (m)
  Esrf,              &! Sublimation rate (kg/m^2/s)
  Gice                ! Heat flux into ice (W/m^2)
  
real, intent(inout) :: &
  Melt,              &! Surface melt rate (kg/m^2/s)
  RofI,              &! Runoff (kg/m^2)
  Tice(Nice)          ! Ice layer temperatures (K)
  
real, intent(out) :: &
  Mice                ! Ice mass change (kg/m^2)

integer :: &
  k                   ! Ice layer counter

real :: &
  a(Nice),           &! Below-diagonal matrix elements
  b(Nice),           &! Diagonal matrix elements
  c(Nice),           &! Above-diagonal matrix elements
  dTice(Nice),       &! Temperature increments (K)
  Gs(Nice),          &! Thermal conductivity between layers (W/m^2/k)
  rhs(Nice)           ! Matrix equation rhs

do k = 1, Nice - 1
  Gs(k) = 2*hcon_ice / (Dice(k) + Dice(k+1))
end do
a(1) = 0
b(1) = rho_ice*hcap_ice*Dice(1) + Gs(1)*dt
c(1) = - Gs(1)*dt
rhs(1) = (Gice - Gs(1)*(Tice(1) - Tice(2)))*dt
do k = 2, Nice - 1
  a(k) = c(k-1)
  b(k) = rho_ice*hcap_ice*Dice(k) + (Gs(k-1) + Gs(k))*dt
  c(k) = - Gs(k)*dt
  rhs(k) = Gs(k-1)*(Tice(k-1) - Tice(k))*dt  &
           + Gs(k)*(Tice(k+1) - Tice(k))*dt 
end do
k = Nice
Gs(k) = hcon_ice/Dice(k)
a(k) = c(k-1)
b(k) = rho_ice*hcap_ice*Dice(k) + (Gs(k-1) + Gs(k))*dt
c(k) = 0
rhs(k) = Gs(k-1)*(Tice(k-1) - Tice(k))*dt
call TRIDIAG(Nice,Nice,a,b,c,rhs,dTice)
Mice = (Melt + Esrf)*dt
! DNG in Main branch, Roff is not set to Melt as it is here
RofI = Melt*dt

do k = 1, Nice
 Tice(k) = Tice(k) + dTice(k)
! DNG from conv with Richard, the lines below is not needed
! if (Tice(k) > Tm) then
!    Melt = rho_ice*hcap_ice*Dice(k)*(Tice(k) - Tm)/Lf
!    Mice = Mice + Melt*dt
!    RofI = RofI + Melt*dt
!    Tice(k) = Tm
! end if
end do

end subroutine ICE

!-----------------------------------------------------------------------
! Saturation specific humidity
!-----------------------------------------------------------------------
subroutine QSAT(P,T,Qs)

use Constants, only : &
  eps,               &! Ratio of molecular weights of water and dry air
  e0,                &! Saturation vapour pressure at Tm (Pa)
  Tm                  ! Melting point (K)

implicit none

real, intent(in) :: &
  P,                 &! Air pressure (Pa)
  T                   ! Temperature (K)

real, intent(out) :: &
  Qs                  ! Saturation specific humidity

real :: &
  Tc,                &! Temperature (C)
  es                  ! Saturation vapour pressure (Pa)

Tc = T - Tm
if (Tc > 0) then
  es = e0*exp(17.5043*Tc / (241.3 + Tc))
else
  es = e0*exp(22.4422*Tc / (272.186 + Tc))
end if
Qs = eps*es / P

end subroutine QSAT

!-----------------------------------------------------------------------
! Snow thermodynamics and hydrology
!-----------------------------------------------------------------------
subroutine SNOW(Nice,Nsmx,Dice,Dmin,Gsrf,Rf,Sf,Ta,Tice,                &
                Dsnw,Esrf,Melt,Nsnw,Sice,Sliq,Tsnw,                    &
                Gice,RofS,snd,SWE)
 
! DNG everywhere in ICE() replace Roff by RofS

use CONSTANTS, only : &
  hcap_ice,          &! Specific heat capacity of ice (J/K/kg)
  hcap_wat,          &! Specific heat capacity of water (J/K/kg)
  hcon_ice,          &! Thermal conductivity of ice (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  rho_wat,           &! Density of water (kg/m^3)
  Tm                  ! Melting point (K)

use DRIVING, only : &
  dt                  ! Timestep (s)

use PARAMETERS, only : &
  bthr,              &! Snow thermal conductivity exponent
  rcld,              &! Maximum density for cold snow (kg/m^3)
  rhof,              &! Fresh snow density (kg/m^3)
  rmlt,              &! Maximum density for melting snow (kg/m^3)
  trho,              &! Snonw compaction time scale (h)
  rho_ice,           &! Density of ice (kg/m^3)
  Wirr                ! Irreducible liquid water content of snow

implicit none

integer, intent(in) :: &
  Nice,              &! Number of ice layers
  Nsmx                ! Maximum number of snow layers

real, intent(in) :: &
  Dice(Nice),        &! Ice layer thicknesses (m)
  Dmin(Nsmx),        &! Minimum snow layer thicknesses (m)
  Tice(Nice),        &! Ice layer temperatures (K)
  Gsrf,              &! Heat flux into surface (W/m^2)
  Rf,                &! Rainfall rate (kg/m^2/s)
  Sf,                &! Snowfall rate (kg/m^2/s)
  Ta                  ! Air temperature (K)
  
integer, intent(inout) :: &
  Nsnw                ! Number of snow layers
  
real, intent(inout) :: &
  Dsnw(Nsmx),        &! Snow layer thicknesses (m)
  Sice(Nsmx),        &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmx),        &! Liquid content of snow layers (kg/m^2)
  Tsnw(Nsmx),        &! Snow layer temperatures (K)
  Esrf,              &! Sublimation rate (kg/m^2/s)
  Melt                ! Surface melt rate (kg/m^2/s)

real, intent(out) :: &
  Gice,              &! Heat flux into ice (W/m^2)
  RofS,              &! Runoff (kg/m^2)
  snd,               &! Snow depth (m)
  SWE                 ! Snow water equivalent (kg/m^2) 

real :: &
  a(Nsmx),           &! Below-diagonal matrix elements
  b(Nsmx),           &! Diagonal matrix elements
  c(Nsmx),           &! Above-diagonal matrix elements
  csnw(Nsmx),        &! Areal heat capacity of snow (J/K/m^2)
  dTs(Nsmx),         &! Temperature increments (k)
  D(Nsmx),           &! Layer thickness before adjustment (m)
  E(Nsmx),           &! Energy contents before adjustment (J/m^2)
  Gs(Nsmx),          &! Thermal conductivity between layers (W/m^2/k)
  ksnw(Nsmx),        &! Thermal conductivity of snow (W/m/K)
  rhs(Nsmx),         &! Matrix equation rhs
  S(Nsmx),           &! Ice contents before adjustment (kg/m^2)
  U(Nsmx),           &! Layer internal energy contents (J/m^2)
  W(Nsmx)             ! Liquid contents before adjustment (kg/m^2)

real :: &
  coldcont,          &! Layer cold content (J/m^2)
  dnew,              &! New snow layer thickness (m)
  dSice,             &! Change in layer ice content (kg/m^2)
  phi,               &! Porosity
  rhos,              &! Density of snow layer (kg/m^3)
  SliqMax,           &! Maximum liquid content for layer (kg/m^2)
  tau,               &! Snow compaction timescale (s)
  wt                  ! Layer weighting

integer :: & 
  k,                 &! Snow layer pointer
  knew,              &! New snow layer pointer
  kold,              &! Old snow layer pointer
  Nold                ! Previous number of snow layers

Gice = Gsrf
RofS = Rf*dt

if (Nsnw > 0) then   ! Existing snwpack

! Heat capacity
  do k = 1, Nsnw
    csnw(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
    rhos = rhof
    if (Dsnw(k) > epsilon(Dsnw)) rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
    ksnw(k) = hcon_ice*(rhos/rho_ice)**bthr
  end do

! Heat conduction
  if (Nsnw == 1) then
    Gs(1) = 2 / (Dsnw(1)/ksnw(1) + Dice(1)/hcon_ice)
    dTs(1) = (Gsrf + Gs(1)*(Tice(1) - Tsnw(1)))*dt /  &
             (csnw(1) + Gs(1)*dt)
  else
    do k = 1, Nsnw - 1
      Gs(k) = 2 / (Dsnw(k)/ksnw(k) + Dsnw(k+1)/ksnw(k+1))
    end do
    a(1) = 0
    b(1) = csnw(1) + Gs(1)*dt
    c(1) = - Gs(1)*dt
    rhs(1) = (Gsrf - Gs(1)*(Tsnw(1) - Tsnw(2)))*dt
    do k = 2, Nsnw - 1
      a(k) = c(k-1)
      b(k) = csnw(k) + (Gs(k-1) + Gs(k))*dt
      c(k) = - Gs(k)*dt
      rhs(k) = Gs(k-1)*(Tsnw(k-1) - Tsnw(k))*dt  &
               + Gs(k)*(Tsnw(k+1) - Tsnw(k))*dt 
    end do
    k = Nsnw
    Gs(k) = 2 / (Dsnw(k)/ksnw(k) + Dice(1)/hcon_ice)
    a(k) = c(k-1)
    b(k) = csnw(k) + (Gs(k-1) + Gs(k))*dt
    c(k) = 0
    rhs(k) = Gs(k-1)*(Tsnw(k-1) - Tsnw(k))*dt  &
             + Gs(k)*(Tice(1) - Tsnw(k))*dt
    call TRIDIAG(Nsnw,Nsmx,a,b,c,rhs,dTs)
  end if 
  do k = 1, Nsnw
    Tsnw(k) = Tsnw(k) + dTs(k)
  end do
  Gice = Gs(Nsnw)*(Tsnw(Nsnw) - Tice(1))

! Convert melting ice to liquid water
  dSice = Melt*dt
  Melt = Melt - sum(Sice(1:Nsnw))/dt
  Melt = max(Melt, 0.)
  do k = 1, Nsnw
    coldcont = csnw(k)*(Tm - Tsnw(k))
    if (coldcont < 0) then
      dSice = dSice - coldcont/Lf
      Tsnw(k) = Tm
    end if
    if (dSice > 0) then
      if (dSice > Sice(k)) then  ! Layer melts completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sliq(k) = Sliq(k) + Sice(k)
        Sice(k) = 0
      else                       ! Layer melts partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        Sliq(k) = Sliq(k) + dSice
        dSice = 0                ! Melt exhausted
      end if
    end if
  end do

! Remove snow by sublimation 
  dSice = max(Esrf, 0.)*dt
  if (dSice > 0) then
    Esrf = Esrf - sum(Sice(1:Nsnw))/dt
    Esrf = max(Esrf, 0.)
    do k = 1, Nsnw
      if (dSice > Sice(k)) then  ! Layer sublimates completely
        dSice = dSice - Sice(k)
        Dsnw(k) = 0
        Sice(k) = 0
      else                       ! Layer sublimates partially
        Dsnw(k) = (1 - dSice/Sice(k))*Dsnw(k)
        Sice(k) = Sice(k) - dSice
        dSice = 0                ! Sublimation exhausted
      end if
    end do
  end if

! Snow hydraulics
  do k = 1, Nsnw
    phi = 0
    if (Dsnw(k) > epsilon(Dsnw)) phi = 1 - Sice(k)/(rho_ice*Dsnw(k))
    SliqMax = rho_wat*Dsnw(k)*phi*Wirr
    Sliq(k) = Sliq(k) + RofS
    RofS = 0
    if (Sliq(k) > SliqMax) then  ! Liquid capacity exceeded
      RofS = Sliq(k) - SliqMax   ! so drainage to next layer
      Sliq(k) = SliqMax
    end if
    coldcont = csnw(k)*(Tm - Tsnw(k))
    if (coldcont > 0) then       ! Liquid can freeze
      dSice = min(Sliq(k), coldcont/Lf)
      Sliq(k) = Sliq(k) - dSice
      Sice(k) = Sice(k) + dSice
      Tsnw(k) = Tsnw(k) + Lf*dSice/csnw(k)
    end if
  end do

! Snow compaction
  tau = 3600*trho
  do k = 1, Nsnw
    if (Dsnw(k) > epsilon(Dsnw)) then
      rhos = (Sice(k) + Sliq(k)) / Dsnw(k)
      if (Tsnw(k) >= Tm) then
          if (rhos < rmlt) rhos = rmlt + (rhos - rmlt)*exp(-dt/tau)
      else
          if (rhos < rcld) rhos = rcld + (rhos - rcld)*exp(-dt/tau)
      end if
      Dsnw(k) = (Sice(k) + Sliq(k)) / rhos
    end if
  end do

end if  ! Existing snowpack

! Add snowfall and frost to layer 1
dSice = Sf*dt - min(Esrf, 0.)*dt
Dsnw(1) = Dsnw(1) + dSice / rhof
Sice(1) = Sice(1) + dSice
Esrf = max(Esrf, 0.)

! New snowpack
if (Nsnw == 0 .and. Sice(1) > 0) then
  Nsnw = 1
  Tsnw(1) = min(Ta, Tm)
end if

! Calculate snow depth and SWE
snd = 0
SWE = 0
do k = 1, Nsnw
  snd = snd + Dsnw(k)
  SWE = SWE + Sice(k) + Sliq(k)
end do

! Store state of old layers
D(:) = Dsnw(:)
S(:) = Sice(:)
W(:) = Sliq(:)
do k = 1, Nsnw
  csnw(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
  E(k) = csnw(k)*(Tsnw(k) - Tm)
end do
Nold = Nsnw

! Initialise new layers
Dsnw(:) = 0
Sice(:) = 0
Sliq(:) = 0
Tsnw(:) = Tm
U(:) = 0
Nsnw = 0

if (SWE > 0) then  ! Existing or new snowpack

! Re-assign and count snow layers
  dnew = snd
  Dsnw(1) = dnew
  k = 1
  if (Dsnw(1) > Dmin(1)) then 
    do k = 1, Nsmx
      Dsnw(k) = Dmin(k)
      dnew = dnew - Dmin(k)
      if (dnew <= Dmin(k) .or. k == Nsmx) then
        Dsnw(k) = Dsnw(k) + dnew
        exit
      end if
    end do
  end if
  Nsnw = k

! Fill new layers from the top downwards
  knew = 1
  dnew = Dsnw(1)
  do kold = 1, Nold
    do
      if (D(kold) < dnew) then
! Transfer all snow from old layer and move to next old layer
        Sice(knew) = Sice(knew) + S(kold)
        Sliq(knew) = Sliq(knew) + W(kold)
        U(knew) = U(knew) + E(kold)
        dnew = dnew - D(kold)
        exit
      else
! Transfer some snow from old layer and move to next new layer
        wt = dnew / D(kold)
        Sice(knew) = Sice(knew) + wt*S(kold) 
        Sliq(knew) = Sliq(knew) + wt*W(kold)
        U(knew) = U(knew) + wt*E(kold)
        D(kold) = (1 - wt)*D(kold)
        E(kold) = (1 - wt)*E(kold)
        S(kold) = (1 - wt)*S(kold)
        W(kold) = (1 - wt)*W(kold)
        knew = knew + 1
        if (knew > Nsnw) exit
        dnew = Dsnw(knew)
      end if
    end do
  end do

! Diagnose snow layer temperatures
  do k = 1, Nsnw
    csnw(k) = Sice(k)*hcap_ice + Sliq(k)*hcap_wat
    if (csnw(k) > epsilon(csnw)) Tsnw(k) = Tm + U(k) / csnw(k)
  end do

end if  ! Existing or new snowpack

end subroutine SNOW

!-----------------------------------------------------------------------
! Solve surface energy balance
!-----------------------------------------------------------------------
subroutine SURFACE(Nice,Nsmx,Dice,Dsnw,LW,Ps,Qa,Sf,Sice,Sliq,SW,Ta,    &
                   Tice,Tsnw,Ua,albs,Tsrf,Esrf,Gsrf,Melt)

use CONSTANTS, only: &
  cp,                &! Specific heat capacity of dry air (J/K/kg)
  g,                 &! Acceleration due to gravity (m/s^2)
  hcon_ice,          &! Thermal conductivity of ice (W/m/K)
  Lf,                &! Latent heat of fusion (J/kg)
  Ls,                &! Latent heat of sublimation (J/kg)
  Rgas,              &! Gas constant for dry air (J/K/kg)
  Rwat,              &! Gas constant for water vapour (J/K/kg)
  sb,                &! Stefan-Boltzmann constant (W/m^2/K^4)
  Tm,                &! Melting point (K)
  vkman               ! Von Karman constant

use DRIVING, only: &
  dt,                &! Timestep (s)
  zT,                &! Temperature and humidity measurement height (m)
  zU                  ! Wind speed measurement height (m)

use PARAMETERS, only: &
  aice,              &! Ice albedo
  asmx,              &! Maximum albedo for fresh snow
  asmn,              &! Minimum albedo for melting snow
  bstb,              &! Stability slope parameter
  bthr,              &! Snow thermal conductivity exponent
  hfsn,              &! Snow cover fraction depth scale (m)
  rhof,              &! Fresh snow density (kg/m^3)
  Salb,              &! Snowfall to refresh albedo (kg/m^2)
  tcld,              &! Cold snow albedo decay timescale (h)
  tmlt,              &! Melting snow albedo decay timescale (h)
  z0ic,              &! Ice surface roughness length (m)  
  rho_ice,           &! Density of ice (kg/m^3)
  z0sn                ! Snow surface roughness length (m)

implicit none

integer, intent(in) :: &
  Nice,              &! Number of ice layers
  Nsmx                ! Maximum number of snow layers  

real, intent(in) :: &
  LW,                &! Incoming longwave radiation (W/m2)
  Ps,                &! Surface pressure (Pa)
  Qa,                &! Specific humidity (kg/kg)
  Sf,                &! Snowfall rate (kg/m2/s)
  SW,                &! Incoming shortwave radiation (W/m2)
  Ta,                &! Air temperature (K)
  Ua                  ! Wind speed (m/s)
  
real, intent(in) :: &
  Dice(Nice),        &! Ice layer thicknesses (m)
  Dsnw(Nsmx),        &! Snow layer thicknesses (m)
  Sice(Nsmx),        &! Ice content of snow layers (kg/m^2)
  Sliq(Nsmx),        &! Liquid content of snow layers (kg/m^2)  
  Tice(Nice),        &! Ice layer temperatures (K)
  Tsnw(Nsmx)          ! Snow layer temperatures (K)

real, intent(inout) :: &
  albs,              &! Snow albedo
  Tsrf                ! Surface skin temperature (K)

real, intent(out) :: &
  Esrf,              &! Sublimation rate (kg/m^2/s)
  Gsrf,              &! Heat flux into surface (W/m^2)
  Melt                ! Melt rate (kg/m^2/s)

real :: &
  alb,               &! Albedo
  alim,              &! Limiting albedo
  CD,                &! Drag coefficient
  CH,                &! Transfer coefficient for heat and moisture
  D,                 &! dQsat/dT (1/K)
  dE,                &! Change in surface moisture flux (kg/m^2/s)
  dG,                &! Change in surface heat flux (W/m^2)
  dTs,               &! Change in surface skin temperatures (K)
  Dz1,               &! Surface layer thickness (m)
  fh,                &! Stability correction
  fsnw,              &! Snow cover fraction
  Hsrf,              &! Surface sensible heat flux (W/m^2)
  ksnw,              &! Thermal conductivity of snow (W/m/K)
  ksrf,              &! Surface layer thermal conductivity (W/m/K)
  Qs,                &! Saturation humidity at surface layer temperature
  rho,               &! Air density (kg/m^3)
  rhos,              &! Snow density (kg/m^3)
  RiB,               &! Bulk Richardson number
  rKH,               &! rho*CH*Ua (kg/m^2/s)
  Rnet,              &! Surface net radiation (W/m^2)
  rt,                &! Reciprocal timescale for albedo adjustment (1/s)
  snd,               &! Snow depth (m)
  tau,               &! Snow albedo decay timescale (s) 
  Ts1,               &! Surface layer temperature (K)
  z0,                &! Roughness length for momentum (m)
  z0h                 ! Roughness length for heat and moisture (m)
  
! Partial snow cover
snd = sum(Dsnw)
fsnw = tanh(snd/hfsn)
  
! Albedo
tau = 3600*tcld
if (Tsrf >= Tm) tau = 3600*tmlt
rt = 1/tau + Sf/Salb
alim = (asmn/tau + Sf*asmx/Salb)/rt
albs = alim + (albs - alim)*exp(-rt*dt)
alb = fsnw*albs + (1 - fsnw)*aice

! Surface layer temperature and thermal conductivity
rhos = rhof
if (Dsnw(1) > epsilon(Dsnw)) rhos = (Sice(1) + Sliq(1)) / Dsnw(1)
ksnw = hcon_ice*(rhos/rho_ice)**bthr
Dz1 = max(Dice(1), Dsnw(1))
Ts1 = Tice(1) + (Tsnw(1) - Tice(1))*Dsnw(1)/Dice(1)
ksrf = Dice(1) / (2*Dsnw(1)/ksnw + (Dice(1) - 2*Dsnw(1))/hcon_ice)
if (Dsnw(1) > 0.5*Dice(1)) ksrf = ksnw
if (Dsnw(1) > Dice(1)) Ts1 = Tsnw(1)
 
! Surface exchange coefficient
z0 = (z0sn**fsnw) * (z0ic**(1 - fsnw))
z0h = 0.1*z0
CD = (vkman / log(zU/z0))**2
CH = vkman**2 / (log(zU/z0)*log(zT/z0h))
RiB = g*(Ta - Tsrf)*zU**2 / (zT*Ta*Ua**2)
if (RiB > 0) then 
  fh = 1/(1 + 3*bstb*RiB*sqrt(1 + bstb*RiB))
else
  fh = 1 - 3*bstb*RiB / (1 + 3*bstb**2*CD*sqrt(-RiB*zU/z0))
end if
CH = fh*CH  
rho = Ps / (Rgas*Ta)
rKH = rho*CH*Ua

! Surface energy balance without melt
call QSAT(Ps,Tsrf,Qs)
D = Ls*Qs/(Rwat*Tsrf**2)
Esrf = rKH*(Qs - Qa)
Gsrf = 2*ksrf*(Tsrf - Ts1)/Dz1
Hsrf = cp*rKH*(Tsrf - Ta)
Melt = 0
Rnet = (1 - alb)*SW + LW - sb*Tsrf**4
dTs = (Rnet - Hsrf - Ls*Esrf - Gsrf) / &
      ((cp + Ls*D)*rKH + 2*ksrf/Dz1 + 4*sb*Tsrf**3)
dE = rKH*D*dTs
dG = 2*ksrf*dTs/Dz1

! Surface melting
if (Tsrf + dTs > Tm) then
  call QSAT(Ps,Tm,Qs)
  Esrf = rKH*(Qs - Qa)  
  Gsrf = 2*ksrf*(Tm - Ts1)/Dz1
  Hsrf = cp*rKH*(Tm - Ta)
  Rnet = (1 - alb)*SW + LW - sb*Tm**4
  Melt = (Rnet - Hsrf - Ls*Esrf - Gsrf) / Lf
  dE = 0
  dG = 0
  dTs = Tm - Tsrf
end if

! Update surface temperature and fluxes
Tsrf = Tsrf + dTs
Esrf = Esrf + dE
Gsrf = Gsrf + dG

end subroutine SURFACE

!-----------------------------------------------------------------------
! Solve tridiagonal matrix equation
!-----------------------------------------------------------------------
subroutine TRIDIAG(Nvec,Nmax,a,b,c,r,x)

implicit none

integer, intent(in) :: &
  Nvec,              &! Vector length
  Nmax                ! Maximum vector length

real, intent(in) :: &
  a(Nmax),           &! Below-diagonal matrix elements
  b(Nmax),           &! Diagonal matrix elements
  c(Nmax),           &! Above-diagonal matrix elements
  r(Nmax)             ! Matrix equation rhs

real, intent(out) :: &
  x(Nmax)             ! Solution vector

integer :: n          ! Loop counter 

! Work space   
real :: beta, gamma(Nvec) 

beta = b(1)
x(1) = r(1) / beta

do n = 2, Nvec
  gamma(n) = c(n-1) / beta
  beta = b(n) - a(n)*gamma(n)
  x(n) = (r(n) - a(n)*x(n-1)) / beta
end do

do n = Nvec - 1, 1, -1
  x(n) = x(n) - gamma(n+1)*x(n+1)
end do
  
end subroutine TRIDIAG

