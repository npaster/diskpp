 add_executable(hyperelasticity2 hyperelasticity2.cpp)
 target_link_libraries(hyperelasticity2 ${LINK_LIBS})
 install(TARGETS hyperelasticity2 RUNTIME DESTINATION bin)

 add_executable(hyperelasticity2_test_aurrichio hyperelasticity2_test_aurrichio.cpp)
 target_link_libraries(hyperelasticity2_test_aurrichio ${LINK_LIBS})
 install(TARGETS hyperelasticity2_test_aurrichio RUNTIME DESTINATION bin)

add_executable(hyperelasticity2_test_convergence_rates hyperelasticity2_test_convergence_rates.cpp)
target_link_libraries(hyperelasticity2_test_convergence_rates ${LINK_LIBS})
install(TARGETS hyperelasticity2_test_convergence_rates RUNTIME DESTINATION bin)


 add_executable(hyperelasticity2_test_locking_2d hyperelasticity2_test_locking_2d.cpp)
 target_link_libraries(hyperelasticity2_test_locking_2d ${LINK_LIBS})
 install(TARGETS hyperelasticity2_test_locking_2d RUNTIME DESTINATION bin)



add_executable(hyperelasticity2_test_cook hyperelasticity2_test_cook.cpp)
target_link_libraries(hyperelasticity2_test_cook ${LINK_LIBS})
install(TARGETS hyperelasticity2_test_cook RUNTIME DESTINATION bin)
