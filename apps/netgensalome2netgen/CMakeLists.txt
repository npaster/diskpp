add_executable(netgensalome2netgen netgensalome2netgen.cpp)
target_link_libraries(netgensalome2netgen ${LINK_LIBS})
install(TARGETS netgensalome2netgen RUNTIME DESTINATION bin)
