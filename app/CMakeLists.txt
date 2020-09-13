# Build all of the applications here.

# Add a program to test writing the objets.
add_executable(cubeWriteTest cubeWriteTest.cxx)
target_link_libraries(cubeWriteTest LINK_PUBLIC cuberecon)
install(TARGETS cubeWriteTest RUNTIME DESTINATION bin)

# Add a program to translate an ERepSim event into simple hits.
add_executable(cubeERepTranslate cubeERepTranslate.cxx)
target_link_libraries(cubeERepTranslate LINK_PUBLIC cuberecon)
install(TARGETS cubeERepTranslate RUNTIME DESTINATION bin)

# Add a program to generate voxel hits.
add_executable(cubeMakeHits3D cubeMakeHits3D.cxx)
target_link_libraries(cubeMakeHits3D LINK_PUBLIC cuberecon)
install(TARGETS cubeMakeHits3D RUNTIME DESTINATION bin)
