module ode
use mydata
implicit none
contains
subroutine CompMat_ODE_Eval(Params_In, Des)
  real(kind = digit), dimension(:), pointer, intent(in) :: Params_In
  complex(kind = digit), intent(out) :: Des

  real(kind = digit) :: t
  complex(kind = digit), dimension(:), pointer :: VectorY, VectorY_Dif
  integer :: Neq
  real(kind = digit) :: t0, t_end
  complex(kind = digit), dimension(:), pointer :: VectorY0
  real(kind = digit) :: EPSABS

  complex(kind = digit), dimension(4) :: VectorPsi1, VectorPsi2
  integer :: istat

  Neq = 6
  allocate(VectorY(1:Neq), VectorY0(1:Neq), VectorY_Dif(1:Neq), stat=istat)

  VectorPsi1(1) = (1.0_digit, 0.0_digit)
  VectorPsi1(2) = (0.0_digit, 0.0_digit)
  VectorPsi1(3) = (0.0_digit, 0.0_digit)
  VectorPsi1(4) = (0.0_digit, 0.0_digit)
  VectorPsi2(1) = (0.0_digit, 0.0_digit)
  VectorPsi2(2) = (0.0_digit, 0.0_digit)
  VectorPsi2(3) = (1.0_digit, 0.0_digit)
  VectorPsi2(4) = (0.0_digit, 0.0_digit)
  call Psi2Y(VectorPsi1, VectorPsi2, VectorY0)
  t0 = 1.0_digit
  t_end = 0.0_digit
  EPSABS = 1E-6_digit

  call C_RK45(t, VectorY, VectorY_Dif, Neq, t0, t_end, VectorY0, EPSABS, Params_In)
  Des = VectorY(1) / maxval(abs(VectorY)) ! normalize

  deallocate(VectorY, VectorY0, VectorY_Dif, stat=istat)

contains
  subroutine Psi2Y(VectorPsi1, VectorPsi2, VectorY)
    complex(kind = digit), dimension(4) :: VectorPsi1, VectorPsi2
    complex(kind = digit), dimension(6) :: VectorY

    VectorY(1) = VectorPsi1(1) * VectorPsi2(2) &
         - VectorPsi1(2) * VectorPsi2(1)
    VectorY(2) = VectorPsi1(1) * VectorPsi2(3) &
         - VectorPsi1(3) * VectorPsi2(1)
    VectorY(3) = VectorPsi1(1) * VectorPsi2(4) &
         - VectorPsi1(4) * VectorPsi2(1)
    VectorY(4) = VectorPsi1(2) * VectorPsi2(3) &
         - VectorPsi1(3) * VectorPsi2(2)
    VectorY(5) = VectorPsi1(2) * VectorPsi2(4) &
         - VectorPsi1(4) * VectorPsi2(2)
    VectorY(6) = VectorPsi1(3) * VectorPsi2(4) &
         - VectorPsi1(4) * VectorPsi2(3)
  end subroutine Psi2Y

end subroutine CompMat_ODE_Eval


subroutine CompMat_Func(t, Y, Y_Dif, Params_In)
  real(kind = digit) :: t
  complex(kind = digit), dimension(:), pointer :: Y
  complex(kind = digit), dimension(:), pointer :: Y_Dif
  real(kind = digit), dimension(:), pointer :: Params_In

  real(kind = digit) :: Alpha, Re
  complex(kind = digit) :: C
  complex(kind = digit), dimension(6, 6) :: MatrixB

  Alpha = Params_In(1)
  Re = Params_In(2)
  C = cmplx(Params_In(3), Params_In(4), digit)

  call SetMatrixB(Alpha, Re, C, t, MatrixB)
  Y_Dif = matmul(MatrixB, Y)

end subroutine CompMat_Func

subroutine SetMatrixB(Alpha, Re, C, y, MatrixB)
  real(kind = digit) :: Alpha, Re
  complex(kind = digit) :: C
  complex(kind = digit), dimension(6, 6) :: MatrixB
  complex(kind = digit) :: A1, A2, A3, A4

  real(kind = digit) :: y
  real(kind = digit) :: Vel, Vel_Dif2
!!!!!!!!!!!!!!!!!!!!
!! Poiseuille Flow Velocity Profile
!!!!!!!!!!!!!!!!!!!!
  Vel = 2.0_digit * y - y**2
  Vel_Dif2 = -2.0_digit
!!!!!!!!!!!!!!!!!!!!
  A1 = (0.0_digit, 0.0_digit)
  A2 = 2.0_digit * Alpha**2 &
       + (0.0_digit, 1.0_digit) * Alpha * Re * (Vel - C)
  A3 = (0.0_digit, 0.0_digit)
  A4 = (-1.0_digit) &
       * ( Alpha**4 + (0.0_digit, 1.0_digit) * Alpha*Re &
       * ( Alpha**2*(Vel-C) + Vel_Dif2 ) )

  MatrixB = (0.0_digit, 0.0_digit)

  MatrixB(1, 2) = (1.0_digit, 0.0_digit)

  MatrixB(2, 3) = (1.0_digit, 0.0_digit)
  MatrixB(2, 4) = (1.0_digit, 0.0_digit)

  MatrixB(3, 5) = (1.0_digit, 0.0_digit)
  MatrixB(3, 1) = A3
  MatrixB(3, 2) = A2
  MatrixB(3, 3) = A1

  MatrixB(4, 5) = (1.0_digit, 0.0_digit)

  MatrixB(5, 6) = (1.0_digit, 0.0_digit)
  MatrixB(5, 1) = -A4
  MatrixB(5, 4) = A2
  MatrixB(5, 5) = A1

  MatrixB(6, 2) = -A4
  MatrixB(6, 4) = -A3
  MatrixB(6, 6) = A1
end subroutine SetMatrixB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Complex Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine C_RK45(t, y, f, Neq, t0, t_end, y0, TOL, Params_In)

  real(kind = digit) :: t
  complex(kind = digit), dimension(:), pointer :: y
  complex(kind = digit), dimension(:), pointer :: f
  integer :: Neq
  real(kind = digit) :: t0, t_end
  complex(kind = digit), dimension(:), pointer :: y0
  real(kind = digit) :: TOL
  real(kind = digit), dimension(:), pointer :: Params_In

  real(kind = digit), dimension(:), pointer :: y_scal
  real(kind = digit) :: t_step_try, t_step, t_step_min, t_step_next
  real(kind = digit) :: Bound_RK

  integer :: istat, iter, NMAX_iter
  
  allocate(y_scal(1:Neq), stat=istat)

  t_step_try = 0.1_digit
  t_step_min = 0.0_digit
  t = t0
  t_step = sign(t_step_try, t_end - t0)
  y = y0
  Bound_RK = 1E50_digit

  NMAX_iter = 10000
  do iter = 1, NMAX_iter
     call CompMat_Func(t, y, f, Params_In)
     y_scal = abs(y) + abs(t_step * f) + TINY
     if ( (t + t_step - t_end)*(t + t_step - t0) > 0 ) then ! if stepsize can overshoot, decrease.
        t_step = t_end - t
     end if
     call C_RKQS(t, y, y_scal, f, Neq, t_step, t_step_next, TOL, Params_In)
!!!!!!!!!!!!!!!!!!!!
!! Output intermediate results
!!!!!!!!!!!!!!!!!!!!
!      write(*, 1001) t, y
! 1001 format("t is: ", ES14.7, "; y is: ", 6ES14.7)
!!!!!!!!!!!!!!!!!!!!
     if ( minval(abs(y)) > Bound_RK ) then ! sometimes, the ODE integration blows up due to the improper parameters.
        write(*, *) "The ODE integration is finished."
        write(*, *) "However, the result blows up!"
        return
     end if
     if ( (t - t_end)*(t_end - t0) >= 0 ) then ! reach the end and so the integration is done
        ! write(*, *) "The ODE integration is finished successfully!"
        deallocate(y_scal, stat=istat)
        return ! normal exit
     end if
     if ( abs(t_step_next) < t_step_min ) then
        write(*, *) 'stepsize is smaller than minimum'
        read(*, *)
     end if
     t_step = t_step_next
  end do
  write(*, *) 'too many steps!'
  read(*, *)
  return
   
end subroutine C_RK45

subroutine C_RKQS(t, y, y_scal, f, Neq, t_step_try, t_step_next, TOL, Params_In)

  real(kind = digit) :: t
  complex(kind = digit), dimension(:), pointer :: y, f
  real(kind = digit), dimension(:), pointer :: y_scal
  integer :: Neq
  real(kind = digit) :: t_step_try, t_step_next, TOL
  real(kind = digit), dimension(:), pointer :: Params_In

  real(kind = digit) :: t_step, t_step_temp
  real(kind = digit) :: ErrMax
  complex(kind = digit), dimension(:), pointer :: y_err, y_temp
  
  real(kind = digit), parameter :: SAFETY=0.9_digit, PGROW=-0.2_digit, PSHRNK=-0.25_digit, ERRCON=1.89e-4_digit
  integer :: i, istat
  
  allocate(y_err(1:Neq), y_temp(1:Neq), stat=istat)
  
  t_step = t_step_try
01 call C_RKCK(t, y, y_temp, y_err, f, Neq, t_step, Params_In)
  ErrMax = 0.0_digit
  do i = 1, Neq
     ErrMax = max(ErrMax, abs(y_err(i)/y_scal(i)))
  end do
  ErrMax = ErrMax / TOL ! scale relative to required tolerance
  if ( ErrMax > 1.0_digit ) then ! truncation error too large, reduce stepsize
     t_step_temp = SAFETY * t_step * (ErrMax**PSHRNK)
     t_step = sign( max( abs(t_step_temp), 0.1_digit*abs(t_step) ), t_step ) ! no more than a factor of 10
     if (t + t_step == t) then
        write(*, *) t, t_step
        write(*, *) y
        write(*, *) 'stepsize underflow in C_RKQS'
        read(*, *)
     end if
     goto 01
  else
     if ( ErrMax > ERRCON ) then
        t_step_next = SAFETY * t_step * (ErrMax**PGROW)
     else
        t_step_next = 5.0_digit * t_step
     end if
     t = t + t_step
     y = y_temp
     deallocate(y_err, y_temp, stat=istat)
     return
  end if
  
end subroutine C_RKQS

subroutine C_RKCK(t, y, y_out, y_err, f, Neq, t_step, Params_In)
  real(kind = digit) :: t
  complex(kind = digit), dimension(:), pointer :: y, y_err, y_out, f
  integer :: Neq
  real(kind = digit) :: t_step
  real(kind = digit), dimension(:), pointer :: Params_In

  complex(kind = digit), dimension(:), pointer :: K1, K2, K3, K4, K5, K6
  real(kind = digit), parameter :: &
       A2 = 0.2_digit, A3 = 0.3_digit, A4 = 0.6_digit, &
       A5= 1.0_digit, A6 = 0.875_digit
  real(kind = digit), parameter :: &
       B21 = 0.2_digit
  real(kind = digit), parameter :: &
       B31 = 3.0_digit / 40.0_digit, B32 = 9.0_digit / 40.0_digit
  real(kind = digit), parameter :: &
       B41 = 0.3_digit, B42 = -0.9_digit, B43 = 1.2_digit
  real(kind = digit), parameter :: &
       B51 = -11.0_digit / 54.0_digit, B52 = 2.5_digit, &
       B53 = -70.0_digit / 27.0_digit, B54 = 35.0_digit / 27.0_digit
  real(kind = digit), parameter :: &
       B61 = 1631.0_digit / 55296.0_digit, B62 = 175.0_digit / 512.0_digit, &
       B63 = 575.0_digit / 13824.0_digit, B64 = 44275.0_digit / 110592.0_digit, &
       B65 = 253.0_digit / 4096.0_digit
  real(kind = digit), parameter :: &
       C1 = 37.0_digit / 378.0_digit, C3 = 250.0_digit / 621.0_digit, &
       C4 = 125.0_digit / 594.0_digit, C6 = 512.0_digit / 1771.0_digit
  real(kind = digit), parameter :: &
       DC1 = C1 - 2825.0_digit / 27648.0_digit, DC3 = C3 - 18575.0_digit / 48384.0_digit, &
       DC4 = C4 - 13525.0_digit / 55296.0_digit, DC5 = -277.0_digit / 14336.0_digit, &
       DC6 = C6 - 0.25_digit

  real(kind = digit) :: t_temp
  complex(kind = digit), dimension(:), pointer :: y_temp
  integer :: istat

  allocate(K1(1:Neq), K2(1:Neq), K3(1:Neq), K4(1:Neq), K5(1:Neq), K6(1:Neq), stat=istat)
  allocate(y_temp(1:Neq), stat=istat)

  ! first step
  K1 = f
  ! second step
  t_temp = t + A2 * t_step
  y_temp = y + B21 * t_step * K1 !f is K1
  call CompMat_Func(t_temp, y_temp, K2, Params_In)
  ! third step
  t_temp = t + A3 * t_step
  y_temp = y + t_step * (B31 * K1 + B32 * K2)
  call CompMat_Func(t_temp, y_temp, K3, Params_In)
  ! fourth step
  t_temp = t + A4 * t_step
  y_temp = y + t_step * (B41 * K1 + B42 * K2 + B43 * K3)
  call CompMat_Func(t_temp, y_temp, K4, Params_In)
  ! fifth step
  t_temp = t + A5 * t_step
  y_temp = y + t_step * (B51 * K1 + B52 * K2 + B53 * K3 + B54 * K4)
  call CompMat_Func(t_temp, y_temp, K5, Params_In)
  ! sixth step
  t_temp = t + A6 * t_step
  y_temp = y + t_step * (B61 * K1 + B62 * K2 + B63 * K3 + B64 * K4 + B65 * K5)
  call CompMat_Func(t_temp, y_temp, K6, Params_In)
  ! accumulate increments with proper weights
  y_out = y + t_step * (C1 * K1 + C3 * K3 + C4 * K4 + C6 * K6)
  ! extimate error as different between fourth and fifth order methods
  y_err = t_step * (DC1*K1 + DC3*K3 + DC4*K4 + DC5*K5 + DC6*K6)

  deallocate(K1, K2, K3, K4, K5, K6, stat=istat)
  deallocate(y_temp, stat=istat)

  return
end subroutine C_RKCK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Real Version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine R_RK45(t, y, f, Neq, t0, t_end, y0, TOL)
!   real(kind = digit) :: t
!   real(kind = digit), dimension(:), pointer :: y
!   real(kind = digit), dimension(:), pointer :: f
!   integer :: Neq
!   real(kind = digit) :: t0, t_end
!   real(kind = digit), dimension(:), pointer :: y0
!   real(kind = digit) :: TOL
!   ! real(kind = digit), dimension(:), pointer :: Params_In

!   real(kind = digit), dimension(:), pointer :: y_scal
!   real(kind = digit) :: t_step_try, t_step, t_step_min, t_step_next
!   real(kind = digit) :: Bound_RK

!   integer :: istat, iter, NMAX_iter
  
!   allocate(y_scal(1:Neq), stat=istat)

!   t_step_try = 0.1_digit
!   t_step_min = 0.0_digit
!   t = t0
!   t_step = sign(t_step_try, t_end - t0)
!   y = y0
!   Bound_RK = 1E4_digit

!   NMAX_iter = 10000
!   do iter = 1, NMAX_iter
! !!!!!!!!!!!!!!!!!!!!
! !! ODE function
! !!!!!!!!!!!!!!!!!!!!
!      call Mean_Profile_Func(t, y, f)
! !!!!!!!!!!!!!!!!!!!!
!      y_scal = abs(y) + abs(t_step * f) + TINY
!      if ( (t + t_step - t_end)*(t + t_step - t0) > 0 ) then ! if stepsize can overshoot, decrease.
!         t_step = t_end - t
!      end if
!      call R_RKQS(t, y, y_scal, f, Neq, t_step, t_step_next, TOL)
! !!!!!!!!!!!!!!!!!!!!
! !! Output intermediate results
! !!!!!!!!!!!!!!!!!!!!
! !      write(*, 1001) t, y
! ! 1001 format("t is: ", F10.7, "; y is: ", ES14.7)
! !!!!!!!!!!!!!!!!!!!!
!      if ( minval(abs(y)) > Bound_RK ) then ! sometimes, the ODE integration blows up due to the improper parameters.
!         write(*, *) "The ODE integration is finished."
!         write(*, *) "However, the result blows up!"
!         return
!      end if
!      if ( (t - t_end)*(t_end - t0) >= 0 ) then ! reach the end and so the integration is done
!         ! write(*, *) "The ODE integration is finished successfully!"
!         deallocate(y_scal, stat=istat)
!         return ! normal exit
!      end if
!      if ( abs(t_step_next) < t_step_min ) then
!         write(*, *) 'stepsize is smaller than minimum'
!         read(*, *)
!      end if
!      t_step = t_step_next
!   end do
!   write(*, *) 'too many steps!'
!   read(*, *)

!   return  
! end subroutine R_RK45

! subroutine R_RKQS(t, y, y_scal, f, Neq, t_step_try, t_step_next, TOL)

!   real(kind = digit) :: t
!   real(kind = digit), dimension(:), pointer :: y, f
!   real(kind = digit), dimension(:), pointer :: y_scal
!   integer :: Neq
!   real(kind = digit) :: t_step_try, t_step_next, TOL
!   ! real(kind = digit), dimension(:), pointer :: Params_In

!   real(kind = digit) :: t_step, t_step_temp
!   real(kind = digit) :: ErrMax
!   real(kind = digit), dimension(:), pointer :: y_err, y_temp
  
!   real(kind = digit), parameter :: SAFETY=0.9_digit, PGROW=-0.2_digit, PSHRNK=-0.25_digit, ERRCON=1.89e-4_digit
!   integer :: i, istat
  
!   allocate(y_err(1:Neq), y_temp(1:Neq), stat=istat)
  
!   t_step = t_step_try
! 01 call R_RKCK(t, y, y_temp, y_err, f, Neq, t_step)
!   ErrMax = 0.0_digit
!   do i = 1, Neq
!      ErrMax = max(ErrMax, abs(y_err(i)/y_scal(i)))
!   end do
!   ErrMax = ErrMax / TOL ! scale relative to required tolerance
!   if ( ErrMax > 1.0_digit ) then ! truncation error too large, reduce stepsize
!      t_step_temp = SAFETY * t_step * (ErrMax**PSHRNK)
!      t_step = sign( max( abs(t_step_temp), 0.1_digit*abs(t_step) ), t_step ) ! no more than a factor of 10
!      if (t + t_step == t) then
!         write(*, *) t, t_step
!         write(*, *) y
!         write(*, *) 'stepsize underflow in R_RKQS'
!         read(*, *)
!      end if
!      goto 01
!   else
!      if ( ErrMax > ERRCON ) then
!         t_step_next = SAFETY * t_step * (ErrMax**PGROW)
!      else
!         t_step_next = 5.0_digit * t_step
!      end if
!      t = t + t_step
!      y = y_temp
!      deallocate(y_err, y_temp, stat=istat)
!      return
!   end if

! end subroutine R_RKQS

! subroutine R_RKCK(t, y, y_out, y_err, f, Neq, t_step)
!   real(kind = digit) :: t
!   real(kind = digit), dimension(:), pointer :: y, y_err, y_out, f
!   integer :: Neq
!   real(kind = digit) :: t_step
!   ! real(kind = digit), dimension(:), pointer :: Params_In

!   real(kind = digit), dimension(:), pointer :: K1, K2, K3, K4, K5, K6
!   real(kind = digit), parameter :: &
!        A2 = 0.2_digit, A3 = 0.3_digit, A4 = 0.6_digit, &
!        A5= 1.0_digit, A6 = 0.875_digit
!   real(kind = digit), parameter :: &
!        B21 = 0.2_digit
!   real(kind = digit), parameter :: &
!        B31 = 3.0_digit / 40.0_digit, B32 = 9.0_digit / 40.0_digit
!   real(kind = digit), parameter :: &
!        B41 = 0.3_digit, B42 = -0.9_digit, B43 = 1.2_digit
!   real(kind = digit), parameter :: &
!        B51 = -11.0_digit / 54.0_digit, B52 = 2.5_digit, &
!        B53 = -70.0_digit / 27.0_digit, B54 = 35.0_digit / 27.0_digit
!   real(kind = digit), parameter :: &
!        B61 = 1631.0_digit / 55296.0_digit, B62 = 175.0_digit / 512.0_digit, &
!        B63 = 575.0_digit / 13824.0_digit, B64 = 44275.0_digit / 110592.0_digit, &
!        B65 = 253.0_digit / 4096.0_digit
!   real(kind = digit), parameter :: &
!        C1 = 37.0_digit / 378.0_digit, C3 = 250.0_digit / 621.0_digit, &
!        C4 = 125.0_digit / 594.0_digit, C6 = 512.0_digit / 1771.0_digit
!   real(kind = digit), parameter :: &
!        DC1 = C1 - 2825.0_digit / 27648.0_digit, DC3 = C3 - 18575.0_digit / 48384.0_digit, &
!        DC4 = C4 - 13525.0_digit / 55296.0_digit, DC5 = -277.0_digit / 14336.0_digit, &
!        DC6 = C6 - 0.25_digit

!   real(kind = digit) :: t_temp
!   real(kind = digit), dimension(:), pointer :: y_temp
!   integer :: istat

!   allocate(K1(1:Neq), K2(1:Neq), K3(1:Neq), K4(1:Neq), K5(1:Neq), K6(1:Neq), stat=istat)
!   allocate(y_temp(1:Neq), stat=istat)

! !!!!!!!!!!!!!!!!!!!!
! !! ODE function calls
! !!!!!!!!!!!!!!!!!!!!
!   ! first step
!   K1 = f
!   ! second step
!   t_temp = t + A2 * t_step
!   y_temp = y + B21 * t_step * K1 !f is K1
!   call Mean_Profile_Func(t_temp, y_temp, K2)
!   ! third step
!   t_temp = t + A3 * t_step
!   y_temp = y + t_step * (B31 * K1 + B32 * K2)
!   call Mean_Profile_Func(t_temp, y_temp, K3)
!   ! fourth step
!   t_temp = t + A4 * t_step
!   y_temp = y + t_step * (B41 * K1 + B42 * K2 + B43 * K3)
!   call Mean_Profile_Func(t_temp, y_temp, K4)
!   ! fifth step
!   t_temp = t + A5 * t_step
!   y_temp = y + t_step * (B51 * K1 + B52 * K2 + B53 * K3 + B54 * K4)
!   call Mean_Profile_Func(t_temp, y_temp, K5)
!   ! sixth step
!   t_temp = t + A6 * t_step
!   y_temp = y + t_step * (B61 * K1 + B62 * K2 + B63 * K3 + B64 * K4 + B65 * K5)
!   call Mean_Profile_Func(t_temp, y_temp, K6)
!   ! accumulate increments with proper weights
!   y_out = y + t_step * (C1 * K1 + C3 * K3 + C4 * K4 + C6 * K6)
!   ! extimate error as different between fourth and fifth order methods
!   y_err = t_step * (DC1*K1 + DC3*K3 + DC4*K4 + DC5*K5 + DC6*K6)
! !!!!!!!!!!!!!!!!!!!!

!   deallocate(K1, K2, K3, K4, K5, K6, stat=istat)
!   deallocate(y_temp, stat=istat)

!   return
! end subroutine R_RKCK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! subroutine Profile_ODE_Eval(t_target, y_target, f_target)
!   real(kind = digit) :: t_target
!   real(kind = digit), dimension(:), pointer :: y_target
!   real(kind = digit), dimension(:), pointer :: f_target

!   real(kind = digit) :: t
!   real(kind = digit), dimension(:), pointer :: y
!   real(kind = digit), dimension(:), pointer :: f
!   integer :: Neq
!   real(kind = digit) :: t0, t_end
!   real(kind = digit), dimension(:), pointer :: y0
!   real(kind = digit) :: EPSABS

!   integer :: istat

!   Neq = 3
!   allocate(y(1:Neq), f(1:Neq), y0(1:Neq), stat=istat)

!   t0 = 0.0_digit
!   y0(1) = 0.0_digit
!   y0(2) = 0.58787095546722412_digit
!   y0(3) = 0.19993467877941015_digit
!   EPSABS = 1E-6_digit
!   t_end = t_target

!   call R_RK45(t, y, f, Neq, t0, t_end, y0, EPSABS)

!   y_target = y
!   f_target = f

!   deallocate(y, f, y0, stat=istat)

! end subroutine Profile_ODE_Eval

! subroutine Mean_Profile_Func(t, y, f)
!   real(kind = digit) :: t
!   real(kind = digit), dimension(:), pointer :: y
!   real(kind = digit), dimension(:), pointer :: f
  
!   f(1) = y(2)
!   f(2) = y(3)
!   f(3) = -0.5_digit * y(1) * y(3)
  
! end subroutine Mean_Profile_Func



end module ode


