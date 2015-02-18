!Perform an inverse of the vec() operator by reshaping
!vector A into a matrix of dimensions n x o 
B = TRANSPOSE(RESHAPE(A, (/ o, n /)))
    
    
do i = 1, o
  do j = 1, n
    if(B(j, i) /= 0.0_dp) then
      C(j) = B(j, i)
      exit
    else
      continue
    end if
  end do
end do
    
!Normalise 
    
!Calculate length
do i = 1, n
 magnitude = magnitude + C(i)*C(i)
end do    
    
magnitude = sqrt(magnitude)

    
!Divide components by length
C = C / magnitude
    
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