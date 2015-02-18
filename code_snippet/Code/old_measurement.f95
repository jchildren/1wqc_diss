!Gets the state vector of the target qubit and outputs it to work vector
call get_vector_complex(state_in, work_vector, state_out, n, trgt)
    
!Calculates inner product between measurement result and target qubit
!Uses either BLAS ZDOTC or intrinsic function depending on preprocessor
#ifdef lblas

 m_value = ZDOTC(2**n, m_vector, 1, work_vector, 1)

#else

  m_value = DOT_PRODUCT(m_vector, work_vector)

#endif

!Multiplies result of dot product with rest of state vector
state_out = m_value * state_out
