add_executable(CH_mixte CH_mixte.cpp)
target_link_libraries(CH_mixte ${LINK_LIBS})
install(TARGETS CH_mixte RUNTIME DESTINATION bin)