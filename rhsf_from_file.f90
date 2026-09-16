! rhsf_from_file.f90
! Recomputes the state u_k of Appendix C (Fig. 3) from a published
! smoothing-parameter file, without any random-number generator.
!
! Input file (real_high_s.txt or real_low_s.txt): 841 lines,
!   column 1 = normalised threshold tau_j, column 2 = s_j,
!   rows ordered j = 841, 840, ..., 1  (tau_841 = 1 closes the nest).
!
! For each present layer k = 840, ..., 1 the present is substituted
! first, t = tau_k, and the nest k, ..., 841 is evaluated:
!   G_841 = -t + tau_841
!   G_j   = -t + tau_j * H_{s_{j+1}}(G_{j+1}),   j = 840, ..., k
!   u_k   = 1 - H_{s_k}(G_k),   H_s(x) = 1/(1+exp(-2x/s)).
! Arithmetic is quadruple precision (kind=16).
!
! Usage:  ./rhsf_from_file  real_high_s.txt  real_high_s_u.asc
program rhsf_from_file
  implicit none
  integer, parameter :: qp = 16
  integer :: n, i, j, k, ios
  real(kind=qp), allocatable :: tau(:), s(:), u(:)
  real(kind=qp) :: t, h, a, b
  character(len=256) :: fin, fout

  call get_command_argument(1, fin)
  call get_command_argument(2, fout)
  if (len_trim(fin) == 0 .or. len_trim(fout) == 0) then
     print *, 'usage: rhsf_from_file <s-file> <output-file>'
     stop 1
  end if

  ! count lines
  open(11, file=trim(fin), status='old', action='read')
  n = 0
  do
     read(11, *, iostat=ios) a, b
     if (ios /= 0) exit
     n = n + 1
  end do
  rewind(11)

  allocate(tau(n), s(n), u(n))
  do i = n, 1, -1                 ! file rows run j = n, n-1, ..., 1
     read(11, *) tau(i), s(i)
  end do
  close(11)
  print *, 'layers read =', n, '  tau(n) =', real(tau(n)), '  s(n) =', real(s(n))

  do k = n - 1, 1, -1
     t = tau(k)                   ! the present is substituted first
     h = hs(-t + tau(n), s(n))    ! closing layer
     do j = n - 1, k, -1
        h = hs(-t + tau(j)*h, s(j))
     end do
     u(k) = 1.0_qp - h
  end do

  open(12, file=trim(fout), status='replace', action='write')
  do k = n - 1, 1, -1
     write(12, '(2(1x,e14.7))') real(tau(k)), real(u(k))
  end do
  close(12)

contains

  pure function hs(x, sv) result(y)
    real(kind=qp), intent(in) :: x, sv
    real(kind=qp) :: y
    y = 1.0_qp / (1.0_qp + exp(-2.0_qp*x/sv))
  end function hs

end program rhsf_from_file
