add_executable(nonlinear_diffusion nonlinear_diffusion.cpp)
target_link_libraries(nonlinear_diffusion ${LINK_LIBS})
install(TARGETS nonlinear_diffusion RUNTIME DESTINATION bin)

