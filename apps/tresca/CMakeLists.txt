add_executable(tresca tresca.cpp)
target_link_libraries(tresca ${LINK_LIBS})
install(TARGETS tresca RUNTIME DESTINATION bin)

add_executable(tresca_test tresca_test.cpp)
target_link_libraries(tresca_test ${LINK_LIBS})
install(TARGETS tresca_test RUNTIME DESTINATION bin)
