include_directories(${COSMOSIS_BASE_DIR} ${CUBA_INCLUDE_DIR})
link_directories(${COSMOSIS_BASE_DIR}/cosmosis/datablock)

# Richness-selection N-operator modules.  Each .cc is a thin driver
# (10 lines) that typedefs NOperatorSelScalar<W> or NOperatorSelRadial<W>
# from src/models/n_operator_sel_t.hh and invokes the module macro.
# Lensing weights live in lensing_weights.hh (header-only).

set(_NC_SEL_SCALAR_MODULES NumCounts MassWeighted BiasWeighted)
set(_NC_SEL_RADIAL_MODULES Sigma1h SigmaTot DSigma1h DSigmaTot
                            Shear1h ShearTot)

foreach(_mod IN LISTS _NC_SEL_SCALAR_MODULES _NC_SEL_RADIAL_MODULES)
  add_library(${_mod}Sel MODULE ${_mod}.cc)
  target_link_libraries(${_mod}Sel PRIVATE cosmosis utils models
                        ${GSL_LIBRARIES} ${CUBA_LIBRARIES})
  target_include_directories(${_mod}Sel PRIVATE
    ${CMAKE_SOURCE_DIR}/src ${CUBACPP_DIR} ${GSL_INCLUDE_DIRS})
endforeach()
