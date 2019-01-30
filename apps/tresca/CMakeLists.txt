add_executable(tresca tresca.cpp)
target_link_libraries(tresca ${LINK_LIBS})
install(TARGETS tresca RUNTIME DESTINATION bin)
