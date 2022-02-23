macro(add_mpi_test tname np)
add_test(NAME ${tname}_${np}
  COMMAND ${CMAKE_COMMAND} -E env
  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
  ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${tname}> ${MPIEXEC_POSTFLAGS})
endmacro()
