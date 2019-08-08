add_executable(cubemesh cubemesh.cpp)
target_link_libraries(cubemesh ${LINK_LIBS})
install(TARGETS cubemesh RUNTIME DESTINATION bin)

add_executable(createnonconforming createnonconforming.cpp)
target_link_libraries(createnonconforming ${LINK_LIBS})
install(TARGETS createnonconforming RUNTIME DESTINATION bin)

add_executable(salomeConverter salomeConverter)
target_link_libraries(salomeConverter ${LINK_LIBS})
install(TARGETS salomeConverter RUNTIME DESTINATION bin)