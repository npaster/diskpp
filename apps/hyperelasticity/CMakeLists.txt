 add_executable(hyperelasticity hyperelasticity.cpp)
 target_link_libraries(hyperelasticity ${LINK_LIBS})
 install(TARGETS hyperelasticity RUNTIME DESTINATION bin)

add_executable(hyperelasticity_test_aurrichio hyperelasticity_test_aurrichio.cpp)
target_link_libraries(hyperelasticity_test_aurrichio ${LINK_LIBS})
 install(TARGETS hyperelasticity_test_aurrichio RUNTIME DESTINATION bin)

 add_executable(hyperelasticity_test_convergence_rates hyperelasticity_test_convergence_rates.cpp)
 target_link_libraries(hyperelasticity_test_convergence_rates ${LINK_LIBS})
 install(TARGETS hyperelasticity_test_convergence_rates RUNTIME DESTINATION bin)


add_executable(hyperelasticity_test_locking_2d hyperelasticity_test_locking_2d.cpp)
target_link_libraries(hyperelasticity_test_locking_2d ${LINK_LIBS})
install(TARGETS hyperelasticity_test_locking_2d RUNTIME DESTINATION bin)


add_executable(hyperelasticity_test_cook hyperelasticity_test_cook.cpp)
target_link_libraries(hyperelasticity_test_cook ${LINK_LIBS})
install(TARGETS hyperelasticity_test_cook RUNTIME DESTINATION bin)
