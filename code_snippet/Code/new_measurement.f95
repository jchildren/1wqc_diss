!Gets the state vector of the target qubit and outputs it to work vector
call get_vector_complex(state_in, work_vector, state_out, n, trgt, U, S, VT)
    
state_out(:) = 0.0_dp
    
do i = 1, 2**n

  state_out = state_out + (SQRT(2.0_dp) * ZDOTC(2, m_vector, 1, U(:, i), 1) * S(i) * VT(i, :))

end do
