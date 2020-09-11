module root_find
use mydata
use ode
implicit none
contains

subroutine C_RootMuller(Root_Found)
  complex(kind = digit) :: Root_Found

  complex(kind = digit), dimension(3) :: root
  real(kind = digit), dimension(3) :: func_value
  complex(kind = digit) :: q, A, B, C, D, D1, D2
  complex(kind = digit) :: root_init, root_init_perturb
  real(kind = digit) :: TOL_rootfind
  integer :: iter, init_iter, N_Muller_Iter

  N_Muller_Iter = 200
  TOL_rootfind = 1.0E-6_digit
! NOTE: the imaginary part should also be changed!
  root_init = Root_Found
  root_init_perturb = cmplx(1.0E-3_digit, 1.0E-2_digit, digit)
  root(1) = root_init - root_init_perturb
  root(2) = root_init
  root(3) = root_init + root_init_perturb

  ! func_value(1) = cmplx(1.0_digit, 0.0_digit, digit)
  ! func_value(2) = cmplx(1.0_digit, 0.0_digit, digit)
  ! func_value(3) = cmplx(1.0_digit, 0.0_digit, digit)
  
  iter = 1
  init_iter = 1
  do
     ! Muller method succeeds only if the initial guess is relatively good.
     ! So, if the iteration takes too many steps, it would have probably failed.
     ! Then just reset the initial guess and restart the iteration.
!      if ( iter > N_Muller_Iter ) then
!         init_iter = init_iter + 1
!         root_init = root_init + 1.0E-3_digit
!         root(1) = root_init - 1.0E-4_digit
!         root(2) = root_init
!         root(3) = root_init + 1.0E-4_digit
!         iter = 1
!         write(*, *) "Change initial guess of root. Restart the iteration."
!         write(*, 1002) root_init
! 1002    format("Current initial guess of root is: ", 2ES14.6)
!         read(*, *)
!      end if

     call C_In_R_Out(real(root(1)), imag(root(1)), func_value(1))
     call C_In_R_Out(real(root(2)), imag(root(2)), func_value(2))
     call C_In_R_Out(real(root(3)), imag(root(3)), func_value(3))
!      write(*, 1000) iter, root(3), abs(func_value(3))
! 1000 format(I3, 'Step: Temp root is: ', 2ES14.6, "; result value is: ", ES14.6)
     
     if ( abs(func_value(1)) < TOL_rootfind ) then
        ! write(*, *) "Converged. Information is as follows."
        ! write(*, 1001) root(1), abs(func_value(1))
        Root_Found = root(1)
        exit
        else if ( abs(func_value(2)) < TOL_rootfind ) then
           ! write(*, *) "Converged. Information is as follows."
           ! write(*, 1001) root(2), abs(func_value(2))
           Root_Found = root(2)
           exit
        else if ( abs(func_value(3)) < TOL_rootfind ) then
           ! write(*, *) "Converged. Information is as follows."
           ! write(*, 1001) root(3), abs(func_value(3))
           Root_Found = root(3)
           exit
        end If
! 1001    format("C is: ", 2ES20.12, "; Error range of result is: ", ES14.6)

     q = (root(3) - root(2)) / (root(2) - root(1))
     
     A = q * func_value(3) &
          - q * (1.0_digit + q) * func_value(2) &
          + q * q * func_value(1)
     B = (2.0_digit * q + 1.0_digit) * func_value(3) &
          - (1.0_digit + q) * (1.0_digit + q) * func_value(2) &
          + q * q * func_value(1)
     C = (1.0_digit + q) * func_value(3)
     ! D1 = B + sqrt(B*B-4*A*C)
     ! D2 = B - sqrt(B*B-4*A*C)
     D1 = B + sqrt(abs(B*B-4.0_digit*A*C))
     D2 = B - sqrt(abs(B*B-4.0_digit*A*C))
     !~ Select the one leading next value to be more close to the real one.
     if (abs(D1) > abs(D2)) then
        D = D1
     else
        D = D2
     end if
     Root_Found = root(3) - (root(3) - root(2)) * (2.0_digit * C / D)
     !~ If it blow up or 2 of all 3 roots are close enough to each other, we
     !~ need to know it.
     if(isnan(real(Root_Found))) Then 
        write(*, *) 'Root_Found is NaN!'
        read(*, *)
        exit
     end if
     !~ Update the guess of roots.
     root(1) = root(2)
     root(2) = root(3)
     root(3) = Root_Found
     
     iter = iter + 1
  end do

end subroutine C_RootMuller


end module root_find
