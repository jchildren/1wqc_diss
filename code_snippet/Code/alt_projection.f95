    do i = 1, 2**n
      operator_sum(i, i) = 1.0_dp
    end do
    
    do j = 1, total_m
      !READ(200, *) m_number, m_type, m_phase
      
      m_number = j
      m_type = 'X'
      m_phase = 0.0_dp
    
      do i = 1, n
    
	m_outer_product(:,:) = 0.0_dp
          
	if(m_number == i) then
	  call measurement_type(m_vector, m_type, m_phase, m_result)
	  call measurement_to_file(m_number, m_type, i, m_result, m_phase)
	  call outer_product_complex_vector(m_vector, m_vector, m_outer_product)
	else
	  m_outer_product(1, 1) = 1.0_dp
	  m_outer_product(2, 2) = 1.0_dp
	end if

      
	Allocate(operator_out(2**i, 2**i))
      
	if(i == 1) then
	  operator_out = m_outer_product
	else
	  call kronecker_product_complex(operator_in, m_outer_product, operator_out)
	  Deallocate(operator_in)
	end if
      
	if(i < n) then
	  Allocate(operator_in(2**i, 2**i))
	  operator_in = operator_out
	  Deallocate(operator_out)
	end if
      
      
      end do
      
      operator_sum = MATMUL(operator_sum, SQRT(2.0_dp) * operator_out)
      Deallocate(operator_out)

    end do

    !call for double precision complex matrix vector multiply from BLAS
    !Applies cz_matrix operator to state vector and outputs it
    call ZGEMV( 'N', 2**n, 2**n, 1.0_dp, operator_sum, 2**n, &
    state_vector, 1, 0.0_dp, state_vector, 1)

    Deallocate(operator_sum)
    
    return