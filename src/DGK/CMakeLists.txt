SET(DGK_SRC
    dgk.cpp
    powmod.cpp)

set(RUNTIME_LIBS   
	helib
	cryptgmm
    ntl
    pthread
    boost_system
    gmp)

add_executable(bench_dgk_mt ben_dgk_mt.cpp ${DGK_SRC})
target_link_libraries(bench_dgk_mt ${RUNTIME_LIBS})

add_executable(bench_dgk_concurrent ben_dgk_concurrent.cpp ${DGK_SRC})
target_link_libraries(bench_dgk_concurrent ${RUNTIME_LIBS})


