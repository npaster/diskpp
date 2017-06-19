add_executable(pLaplace pLaplace.cpp)
target_link_libraries(pLaplace ${LINK_LIBS})
install(TARGETS pLaplace RUNTIME DESTINATION bin)
