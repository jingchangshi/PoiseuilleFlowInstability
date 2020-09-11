module Intersection
use mydata
use dislin
implicit none
contains

subroutine Contour(X_Array, Y_Array, Z_Mat, Z_Level)
  real(kind = digit), dimension(:), pointer, intent(in) :: X_Array, Y_Array
  real(kind = digit), dimension(:, :), pointer, intent(in) :: Z_Mat
  real(kind = digit), intent(in) :: Z_Level
  
  real(kind = digit) :: XA, XE, XOR, XSTP, YA, YE, YOR, YSTP
  integer, parameter :: NMAX_PTS = 1000000, NMAX_CURVES = 10000
  real(kind = digit), dimension(NMAX_PTS) :: X_PT_Array, Y_PT_Array
  integer, dimension(NMAX_CURVES) :: IRAY
  integer :: N_CURVE, Index_Curve, Start_Curve

  XA = minval(X_Array)
  XE = maxval(X_Array)
  YA = minval(Y_Array)
  YE = maxval(Y_Array)
  XOR = XA
  YOR = YA
  XSTP = (XE - XA) * 1E-1_digit ! improve
  YSTP = (YE - YA) * 1E-1_digit ! improve
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call PAGE(3000, 3000) ! the size of the page, dislin level 0
  call METAFL('EPS') ! metafile format: EPS; deslin level 0
  call DISINI() ! dislin level 0, initialization
  call GRAF(XA, XE, XOR, XSTP, YA, YE, YOR, YSTP) ! dislin level 1, set 2D axis
  
  call SETCLR(150) ! the foreground colour, dislin level 1,2,3
  call MYLINE ( (/10, 10/), 2 ) ! global line style: dash
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 1. X_PT_Array, Y_PT_Array: are returned arrays containing the calculated contour.
!! The arrays can contain several curves.
!! 2. NMAX_PTS: is the maximal number of points that can be passed to 
!! X_PT_Array, Y_PT_Array.
!! 3. IRAY: is a returned integer array that contains the number of points 
!! for each generated contour curve.
!! 4. NMAX_CURVES: is the maximal number of entries that can be passed to IRAY.
!! 5. N_CURVE: is the returned number of generated curves.
  call CONPTS (X_Array, size(X_Array), Y_Array, size(Y_Array), Z_Mat, Z_Level, &
       X_PT_Array, Y_PT_Array, NMAX_PTS, IRAY, NMAX_CURVES, N_CURVE)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Start_Curve = 1
  do Index_Curve = 1, N_CURVE
     call CURVE(X_PT_Array(Start_Curve),Y_PT_Array(Start_Curve),IRAY(Index_Curve))
     Start_Curve = Start_Curve + IRAY(Index_Curve)
  end do
  call COLOR('FORE')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call ENDGRF ! dislin level 1
  ! call RESET('ALL')
  call DISFIN() ! dislin level 0, end
end subroutine Contour
!!==============================================================================
subroutine Contour_Cross(X_Array, Y_Array, Z_Mat, Z_Level, GrowMax_C_R, GrowMax_C_I)
  real(kind = digit), dimension(:), pointer, intent(in) :: X_Array, Y_Array
  complex(kind = digit), dimension(:, :), pointer, intent(in) :: Z_Mat
  real(kind = digit), intent(in) :: Z_Level
  real(kind = digit), intent(out) :: GrowMax_C_R, GrowMax_C_I

  real(kind = digit) :: XA, XE, XOR, XSTP, YA, YE, YOR, YSTP
  integer, parameter :: NMAX_PTS = 10000000, NMAX_CURVES = 10000
  real(kind = digit), dimension(NMAX_PTS) :: &
       X_PT_Array_R, Y_PT_Array_R, &
       X_PT_Array_I, Y_PT_Array_I
  integer, dimension(NMAX_CURVES) :: IRAY_R, IRAY_I
  integer :: N_CURVE_R, N_CURVE_I, Index_R, Index_I, Index_Cross
  integer :: Index_Curve, Start_Curve
  integer, parameter :: NMAX_CrossArray = 1000
  real(kind = digit), dimension(NMAX_CrossArray) :: CrossArray_X, CrossArray_Y
  real(kind = digit) :: MAX_CrossArray_X, MAX_CrossArray_Y
  real(kind = digit) :: TOL_CROSS

  TOL_CROSS = 1E-3_digit
  X_PT_Array_R = -1E3_digit
  X_PT_Array_I = -1E3_digit
  Y_PT_Array_R = -1E3_digit
  Y_PT_Array_I = -1E3_digit

  XA = minval(X_Array)
  XE = maxval(X_Array)
  YA = minval(Y_Array)
  YE = maxval(Y_Array)
  XOR = XA
  YOR = YA
  ! XSTP = (XE - XA) * 5E-1_digit ! improve
  ! YSTP = (YE - YA) * 5E-1_digit ! improve
  XSTP = 1E-1_digit
  YSTP = 1E-1_digit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------
  call PAGE(3000, 3000) ! the size of the page, dislin level 0
  call METAFL('EPS') ! metafile format: EPS; deslin level 0
!!------------------
  call DISINI() ! dislin level 0, initialization
  call LABDIG (2, 'XYZ') ! the number of digits is automatically calculated
  call GRAF(XA, XE, XOR, XSTP, YA, YE, YOR, YSTP) ! dislin level 1, set 2D axis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! CONPTS: level 0, 1, 2, 3
!! 1. X_PT_Array, Y_PT_Array: are returned arrays containing the calculated contour.
!! The arrays can contain several curves.
!! 2. NMAX_PTS: is the maximal number of points that can be passed to 
!! X_PT_Array, Y_PT_Array.
!! 3. IRAY: is a returned integer array that contains the number of points 
!! for each generated contour curve.
!! 4. NMAX_CURVES: is the maximal number of entries that can be passed to IRAY.
!! 5. N_CURVE: is the returned number of generated curves.
  call CONPTS (X_Array, size(X_Array), Y_Array, size(Y_Array), real(Z_Mat), Z_Level, &
       X_PT_Array_R, Y_PT_Array_R, NMAX_PTS, IRAY_R, NMAX_CURVES, N_CURVE_R)
  write(*, 1001) count((X_PT_Array_R >= XA))
1001 format('Total number of points of the green line on which real(Des_Mat) is 0: ', I10)
!!------------------
  call SETRGB(0.0_digit, 1.0_digit, 0.0_digit) ! the foreground colour, dislin level 1,2,3
  call MYLINE ( (/10, 10/), 2 ) ! global line style: dash
  Start_Curve = 1
  do Index_Curve = 1, N_CURVE_R
     call CURVE(X_PT_Array_R(Start_Curve),Y_PT_Array_R(Start_Curve),IRAY_R(Index_Curve))
     Start_Curve = Start_Curve + IRAY_R(Index_Curve)
  end do
  call COLOR('FORE')
!!------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call CONPTS (X_Array, size(X_Array), Y_Array, size(Y_Array), imag(Z_Mat), Z_Level, &
       X_PT_Array_I, Y_PT_Array_I, NMAX_PTS, IRAY_I, NMAX_CURVES, N_CURVE_I)
  write(*, 1002) count((X_PT_Array_I >= XA))
1002 format('Total number of points of the red line on which imag(Des_Mat) is 0: ', I10)
!!------------------
  call SETRGB(1.0_digit, 0.0_digit, 0.0_digit)
  call MYLINE ((/1/), 1)
  Start_Curve = 1
  do Index_Curve = 1, N_CURVE_I
     call CURVE(X_PT_Array_I(Start_Curve),Y_PT_Array_I(Start_Curve),IRAY_I(Index_Curve))
     Start_Curve = Start_Curve + IRAY_I(Index_Curve)
  end do
  call COLOR('FORE')
!!------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get the cross points, stored in CrossArray
!! Actually, cross points are original points.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CrossArray_X = XA - 1E3_digit
  CrossArray_Y = XA - 1E3_digit
  Index_Cross = 1
  do Index_R = 1, count((X_PT_Array_R >= XA))
     do Index_I = 1, count((X_PT_Array_I >= XA))
        if ( (X_PT_Array_R(Index_R) - X_PT_Array_I(Index_I))**2 &
             + (Y_PT_Array_R(Index_R) - Y_PT_Array_I(Index_I))**2 &
             < TOL_CROSS**2 ) then
           CrossArray_X(Index_Cross) = X_PT_Array_I(Index_I)
           CrossArray_Y(Index_Cross) = Y_PT_Array_I(Index_I)
           Index_Cross = Index_Cross + 1
           exit
        end if
     end do
  end do
  write(*, 1003) count((CrossArray_X >= XA))
1003 format('Total number of points of intersection: ', I7)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------------------
  call MARKER (15) ! the symbols used to plot points, 15 is circle
!! INCMRK(NMRK) selects line or symbol mode for CURVE.
!! NMRK = - n means that CURVE plots only symbols. Every n-th point will be marked by a symbol.
  call INCMRK (-1)
  call CURVE(CrossArray_X(1:count((CrossArray_X >= XA))), &
       CrossArray_Y(1:count((CrossArray_X >= XA))), &
       count((CrossArray_X >= XA)))
!!------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  MAX_CrossArray_Y = maxval(CrossArray_Y(1:count((CrossArray_X >= XA))))
  MAX_CrossArray_X = CrossArray_X(maxloc(CrossArray_Y(1:count((CrossArray_X >= XA))), 1))
  GrowMax_C_R = MAX_CrossArray_X
  GrowMax_C_I = MAX_CrossArray_Y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call ENDGRF ! dislin level 1
  call DISFIN() ! dislin level 0, end

end subroutine Contour_Cross


end module
