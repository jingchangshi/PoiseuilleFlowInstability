program main
use mydata
! use output
use ode
! use root_find
use intersection
implicit none

real(kind = digit) :: Alpha, Re
integer :: iter_C_R, iter_C_I, NMAX_C_R, NMAX_C_I
real(kind = digit) :: C_R0, C_I0, C_R_Range, C_I_Range
real(kind = digit), dimension(:), pointer :: C_R_Array, C_I_Array
real(kind = digit), dimension(:), pointer :: Params_In
complex(kind = digit) :: Des
complex(kind = digit), dimension(:,:), pointer :: Des_Mat
real(kind = digit) :: GrowMax_C_R, GrowMax_C_I

integer :: istat
real(kind = digit) :: Time_Start, Time_Stop

call cpu_time(Time_Start)

allocate(Params_In(1:4), stat=istat)

Alpha = 1.0_digit
Re = 1E4_digit
Params_In(1) = Alpha
Params_In(2) = Re
C_R0 = 0.1_digit
C_I0 = -0.1_digit
C_R_Range = 0.2_digit
C_I_Range = 0.2_digit
NMAX_C_R = 100 + 1
NMAX_C_I = 100 + 1
allocate(C_R_Array(1:NMAX_C_R), &
     C_I_Array(1:NMAX_C_I), &
     Des_Mat(1:NMAX_C_R, 1:NMAX_C_I), stat=istat)

do iter_C_R = 1, NMAX_C_R
   write(*, *) iter_C_R
   C_R_Array(iter_C_R) = C_R0 + (iter_C_R-1) * C_R_Range / (NMAX_C_R-1)
   Params_In(3) = C_R_Array(iter_C_R)
   do iter_C_I = 1, NMAX_C_I
      ! write(*, *) iter_C_I
      C_I_Array(iter_C_I) = C_I0 + (iter_C_I-1) * C_I_Range / (NMAX_C_I-1)
      Params_In(4) = C_I_Array(iter_C_I)
      call CompMat_ODE_Eval(Params_In, Des)
!       write(*, 1000) Des
! 1000 format(ES22.15, ' ', ES22.15)
      Des_Mat(iter_C_R, iter_C_I) = Des
   end do
end do

call Contour_Cross(C_R_Array, C_I_Array, Des_Mat, 0.0_digit, &
     GrowMax_C_R, GrowMax_C_I)
write(*, *) GrowMax_C_R, GrowMax_C_I

deallocate(Params_In, stat=istat)
deallocate(C_R_Array, C_I_Array, Des_Mat, stat=istat)

call cpu_time(Time_Stop)
write(*, *) 'Elapsed time is: ', Time_Stop - Time_Start, ' seconds.'

end program main
