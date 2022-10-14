module hrweno_kinds
!! Aux module to define real precision.
   use iso_fortran_env, only: real32, real64, real128
   implicit none
   private

   public :: rk

!#ifdef REAL32
   integer, parameter :: rk = real32
! #elif REAL64
!    integer, parameter :: rk = real64
! #elif REAL128
!    integer, parameter :: rk = real128
! #else
!    integer, parameter :: rk = real64
! #endif

end module hrweno_kinds
