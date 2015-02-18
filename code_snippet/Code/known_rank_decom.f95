!Find the value of D required for matrix B to be formed
!Reversing the kronecker product
do i = 1, o
  do j = 1, n
    if(C(j) /= 0.0_dp) then
      D(i) = (1.0_dp/(C(j))) * B(j, i)
      exit
    else
      continue
    end if
  end do
end do
    