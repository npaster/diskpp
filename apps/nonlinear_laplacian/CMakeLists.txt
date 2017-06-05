add_executable(NL_laplacian NL_laplacian.cpp)
target_link_libraries(NL_laplacian ${LINK_LIBS})
install(TARGETS NL_laplacian RUNTIME DESTINATION bin)

