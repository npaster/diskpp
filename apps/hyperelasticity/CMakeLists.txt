add_executable(hyperelasticity hyperelasticity.cpp)
target_link_libraries(hyperelasticity ${LINK_LIBS})
install(TARGETS hyperelasticity RUNTIME DESTINATION bin)
