module output
use mydata
implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintReal2D(Real2D, RowMax, ColMax)

integer :: RowMax, ColMax
real(kind = digit), dimension(RowMax, ColMax) :: Real2D
integer :: i, j

outer: do i = 1, RowMax
  inner: do j = 1, ColMax
    write(*, 1000, advance = 'no') Real2D(i, j)
    1000 format(ES9.2)
    if (j .eq. ColMax) then
      write(*, *)
    end if
  end do inner
end do outer

end subroutine PrintReal2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine PrintComplex2D(Complex2D, RowMax, ColMax)

integer :: RowMax, ColMax
complex(kind = digit), dimension(RowMax, ColMax) :: Complex2D
integer :: i, j

outer: do i = 1, RowMax
  inner: do j = 1, ColMax
    write(*, 1000, advance = 'no') Complex2D(i, j)
    1000 format(('(', ES9.2, ',', ES9.2, ')'))
    if (j .eq. ColMax) then
      write(*, *)
    end if
  end do inner
end do outer

end subroutine PrintComplex2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteReal2D(Real2D, RowMax, ColMax, filename)

integer :: RowMax, ColMax
real(kind = digit), dimension(RowMax, ColMax) :: Real2D
character(len = *) :: filename
integer :: i, j, ierror

open(unit = 01, file = filename, status = 'replace', action = 'write', iostat = ierror)

outer: do i = 1, RowMax
  inner: do j = 1, ColMax
    write(01, 1000, advance = 'no') Real2D(i, j)
    1000 format(ES23.16, ' ')
    if (j .eq. ColMax) then
      write(01, *) !~begin new line
    end if
  end do inner
end do outer

close(unit = 01)

end subroutine WriteReal2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteComplex2D(Complex2D, RowMax, ColMax, filename)

integer :: RowMax, ColMax
complex(kind = digit), dimension(RowMax, ColMax) :: Complex2D
character(len = *) :: filename
integer :: i, j, ierror

open(unit = 01, file = filename, status = 'replace', action = 'write', iostat = ierror)

outer: do i = 1, RowMax
  inner: do j = 1, ColMax
    write(01, 1000, advance = 'no') Complex2D(i, j)
    1000 format('(', ES17.10, ',', ES17.10, ')', ';')
    if (j .eq. ColMax) then
      write(01, *) !~begin new line
    end if
  end do inner
end do outer

close(unit = 01)

end subroutine WriteComplex2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteInteger2D(Integer2D, RowMax, ColMax, filename)

integer :: RowMax, ColMax
integer, dimension(RowMax, ColMax) :: Integer2D
character(len = *) :: filename
integer :: i, j, ierror

open(unit = 01, file = filename, status = 'replace', action = 'write', iostat = ierror)

outer: do i = 1, RowMax
  inner: do j = 1, ColMax
    write(01, 1000, advance = 'no') Integer2D(i, j)
    1000 format(I3, ',')
    if (j .eq. ColMax) then
      write(01, *) !~begin new line
    end if
  end do inner
end do outer

close(unit = 01)

end subroutine WriteInteger2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ReadReal2D(Real2D, RowMax, ColMax, filename)

  integer :: RowMax, ColMax
  real(kind = digit), dimension(RowMax, ColMax) :: Real2D
  real(kind = digit), dimension(ColMax) :: temp
  character(len = *) :: filename
  integer :: ierror, Row
  Row = 1
  open(unit = 01, file = filename, status = 'old', action = 'read', iostat = ierror)

  openif: if(ierror == 0) then
     readloop: do
        read(01, *, iostat = ierror) temp
        if(ierror /= 0) then
           exit
        end if
        Real2D(Row, :) = temp
        Row = Row + 1
     end do readloop
     readif: if(ierror > 0) then
        write(*, *) 'Read error occured!'
        stop
     else readif
        write(*, *) 'End of the line reached!'
     end if readif
  else openif
     write(*, *) 'Error opening file!'
     stop
  end if openif

  close(unit = 01)

end subroutine ReadReal2D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gnuplot_box(data_file_name, output_format)

character(len = *) :: data_file_name
character(len = 32) :: sys_command
character(len = *) :: output_format
integer :: ierror

open(unit = 01, file = 'gnuplot_command_file.txt', status = 'replace', action = 'write', iostat = ierror)
write(01, *) 'set term ', output_format
write(01, *) 'set boxwidth ', '0.75', ' relative'
write(01, *) 'set style fill solid 1.0'
write(01, *) 'set output ', "'", data_file_name, '.', output_format, "'"
write(01, *) 'plot ', "'", data_file_name, '.dat', "'", 'with boxes'
write(01, *) 'pause -1'

close(unit = 01)
sys_command = 'gnuplot gnuplot_command_file.txt'
call system(sys_command)

end subroutine gnuplot_box
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gnuplot_line(data_file_name, output_format)
character(len = *) :: data_file_name
character(len = 32) :: sys_command
character(len = *) :: output_format
integer :: ierror

open(unit = 01, file = 'gnuplot_command_file.txt', status = 'replace', action = 'write', iostat = ierror)
write(01, *) 'set term ', output_format
write(01, *) 'set style fill solid 1.0'
write(01, *) 'set output ', "'", data_file_name, '.', output_format, "'"
write(01, *) 'plot ', "'", data_file_name, '.dat', "'"
write(01, *) 'pause -1'

close(unit = 01)
sys_command = 'gnuplot gnuplot_command_file.txt'
call system(sys_command)

end subroutine gnuplot_line
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gnuplot_contour(data_file_name, output_format)
character(len = *) :: data_file_name
character(len = 32) :: sys_command
character(len = *) :: output_format
integer :: ierror

open(unit = 01, file = 'gnuplot_command_file.txt', status = 'replace', action = 'write', iostat = ierror)
write(01, *) 'set term ', output_format
write(01, *) 'set output ', "'", data_file_name, '.', output_format, "'"
! write(01, *) 'set auto'
! write(01, *) 'set zrange [-0.3:0.3]'
! write(01, *) 'set parametric'
write(01, *) 'set contour base'
write(01, *) 'set cntrparam levels 0.05'
write(01, *) 'set style data lines'
write(01, *) 'splot ', "'", data_file_name, '.dat', "'"
write(01, *) 'pause -1'

close(unit = 01)
sys_command = 'gnuplot gnuplot_command_file.txt'
call system(sys_command)

end subroutine gnuplot_contour
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module output
