module mydata
implicit none

integer, parameter :: digit = selected_real_kind(15, 100)

real(kind = digit), parameter :: PI = 4.0_digit * atan(1.0_digit)
real(kind = digit), parameter :: TINY = 1E-15_digit

real(kind = digit) :: Alpha_Global

end module mydata
