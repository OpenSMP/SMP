cmake_minimum_required(VERSION 3.2)

include_directories(${GTEST_INCLUDE_DIRS})

set(RUNTIME_LIBS   
    helib
    cryptgmm
    ntl
    pthread
    boost_system
    gmp)

set(CPP_ITESTS
    test_net_io
    fast_decryption
	test_mult_in_slots
    )

#The integration tests must be single source code, and are compiled as a standalone application
message(STATUS ${GTEST_BOTH_LIBRARIES})
foreach (CPP_ITEST ${CPP_ITESTS})
    add_executable(${CPP_ITEST} ${CPP_ITEST}.cpp)
    target_link_libraries(${CPP_ITEST} ${RUNTIME_LIBS} ${GTEST_BOTH_LIBRARIES})
endforeach(CPP_ITEST)
