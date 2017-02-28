add_executable(NL_elasticity NL_elasticity.cpp)
target_link_libraries(NL_elasticity ${LINK_LIBS})
install(TARGETS NL_elasticity RUNTIME DESTINATION bin)

