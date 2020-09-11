module interpolation
use mydata
implicit none
contains

subroutine R_Tab_Ref(X, Y, Ref_Flag)
! Make sure that the minimum step size is larger than the minimum step size in the input data file.
! If step size is integer times of step size in source data file, then the ref module can be further simplified.
  real(kind = digit), intent(in) :: X
  real(kind = digit), intent(out) :: Y
  integer, intent(in) :: Ref_Flag
  
  real(kind = digit), dimension(NMAX_Mean_Profile_ROW) :: X_Tab, Y_Tab
  integer :: Loc_last, Loc_next, Loc

  X_Tab = Source_X_Tab
  if ( Ref_Flag == 0 ) then
     Y_Tab = Source_Vel_Tab
  else if ( Ref_Flag == 1 ) then
     Y_Tab = Source_Vel_Dif_Tab
  else
     write(*, *) 'Ref_Flag is wrong!'
     return
  end if

  Loc = minloc(abs(X_Tab - X), 1) ! find the most approximate value to the index

  if ( X > X_Tab(Loc) ) then
     Loc_last = Loc
     Loc_next = Loc_last + 1
  else if ( X < X_Tab(Loc) ) then
     Loc_next = Loc
     Loc_last = Loc_next - 1
  else
     Loc_next = Loc
     Loc_last = Loc
  end if
  ! Use midpoint interpolation
  Y = 0.5_digit * ( Y_Tab(Loc_last) + Y_Tab(Loc_next) )

end subroutine R_Tab_Ref


end module interpolation
