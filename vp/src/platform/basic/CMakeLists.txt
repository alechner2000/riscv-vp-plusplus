file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.h)

add_library(platform-basic
ethernet.cpp
display.cpp
${HEADERS})

target_include_directories(platform-basic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(platform-basic systemc)

add_executable(riscv-vp
        main.cpp)

target_link_libraries(riscv-vp rv32 platform-basic platform-common gdb-mc ${Boost_LIBRARIES} systemc pthread)

INSTALL(TARGETS riscv-vp RUNTIME DESTINATION bin)


add_executable(riscv-vp4Canny
        main4Canny.cpp)

target_link_libraries(riscv-vp4Canny rv32 platform-basic platform-common gdb-mc ${Boost_LIBRARIES} systemc pthread)

INSTALL(TARGETS riscv-vp4Canny RUNTIME DESTINATION bin)

