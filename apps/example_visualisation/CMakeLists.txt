add_executable(example_visualisation  example_visualisation.cpp)
target_link_libraries(example_visualisation ${LINK_LIBS})
install(TARGETS example_visualisation RUNTIME DESTINATION bin)
