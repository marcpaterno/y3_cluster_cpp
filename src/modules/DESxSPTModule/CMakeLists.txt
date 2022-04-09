include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

add_library(
  AbundanceIntegralMLamtrue_DESxSPTModule MODULE
  AbundanceIntegralMLamtrue_DESxSPTModule.cc
  )
target_link_libraries(
  AbundanceIntegralMLamtrue_DESxSPTModule
  PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  AbundanceIntegralMLamtrue_DESxSPTModule
  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS}
  )

add_library(
  SelectionIntegralDESxSPTModule MODULE SelectionIntegralDESxSPTModule.cc
  )
target_link_libraries(
  SelectionIntegralDESxSPTModule PRIVATE cosmosis utils models ${GSL_LIBRARIES}
                                         ${CUBA_LIBRARIES}
  )
target_include_directories(
  SelectionIntegralDESxSPTModule PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR}
                                         ${GSL_INCLUDE_DIRS}
  )

add_library(
  AbundanceIntegralZtrue_DESxSPTModule MODULE
  AbundanceIntegralZtrue_DESxSPTModule.cc
  )
target_link_libraries(
  AbundanceIntegralZtrue_DESxSPTModule
  PRIVATE cosmosis utils models ${GSL_LIBRARIES} ${CUBA_LIBRARIES}
  )
target_include_directories(
  AbundanceIntegralZtrue_DESxSPTModule
  PRIVATE ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS}
  )
