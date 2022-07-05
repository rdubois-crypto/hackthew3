add_definitions(-DST31)

link_libraries(emu -lcmocka-static)

add_executable(test_pbc test_pbc.c)


foreach(target test_swu_hash)
    add_test(NAME ${target} COMMAND qemu-arm-static ${target} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
