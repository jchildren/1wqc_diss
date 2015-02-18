!Loops over values 1 to 4 with J
!Reflective of the unitary operator which has 4 terms
do j = 1, 4

  !Allocates the secondary work matrix into the initial size
  !for Kronecker product multiplication
  Allocate(out_matrix(2, 2))
  
  !Performs a check to see if the first matrix in the order of
  !multiplication corresponds to a control or target 
  !if they do, assigns the appropriate matrix dependant
  !upon the value of j.
  !Otherwise assigns the identity matrix
  if((trgt.eq.1).and.((j.eq.2).or.(j.eq.4))) then
    out_matrix = z_matrix
  elseif((ctrl.eq.1).and.((j.eq.3).or.(j.eq.4))) then
    out_matrix = z_matrix
  else
    out_matrix = identity
  end if

  !Loops over each individual qubit
  !This will produce a matrix of appropriate size for each
  !term of the unitary operator
  do i = 2, n
      
    !Assigns a size value to the work matrix dependent
    !on the step in the loop, thus allowing it
    !to contain the appropriate size of matrix at this step
    Allocate(work_matrix(2**(i-1), 2**(i-1)))
      
    !Assigns the work matrix the value of the output matrix
    !Takes the value from the output of last step
    !Allowing output matrix to be deallocated
    work_matrix = out_matrix
      
    !Deallocate output matrix
    Deallocate(out_matrix)
       
    !Allocates new size to the output matrix
    !New size is appropriate for storage of
    !Kronecker product between work matrix and a 2x2 matrix
    Allocate(out_matrix(2**i, 2**i))
      
    !Determines whether the next operator of the multiplication
    !Will be a control or target qubit, depending on the term of U
    !Assigns the appropriate value if so for kronecker products
    !Otherwise uses the identity matrix for the kronecker product
    if((ctrl.eq.i).and.((j.eq.3).or.(j.eq.4))) then
	  call kronecker_product_complex(work_matrix, z_matrix, out_matrix)
    elseif((trgt.eq.i).and.((j.eq.2).or.(j.eq.4))) then
	  call kronecker_product_complex(work_matrix, z_matrix, out_matrix)
    else
	  call kronecker_product_complex(work_matrix, identity, out_matrix)
    end if
      
    !Deallocates the work matrix ready
    !For allocation in next loop
    Deallocate(work_matrix)
      
  !Ends do loop for the jth term of the Operator
  end do
    
  !Determines which term of operator the loop is on
  !If j is four, removes output from cz_matrix
  !otherwise adds output to cz_matrix
  !Reflective of signs of unitary operator
  if(j.eq.4) then
    cz_matrix = cz_matrix - out_matrix
  else
    cz_matrix = cz_matrix + out_matrix
  end if
    
  !Deallocates out matrix for next loop of j
  Deallocate(out_matrix)
  
end do
