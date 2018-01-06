SET(SRC
    Matrix.cpp
    literal.cpp
    HElib.cpp
	net_io.cpp
    )

add_library(cryptgmm ${SRC})

set(RUNTIME_LIBS   
    helib
	cryptgmm
    ntl
    pthread
    boost_system
    gmp)

add_executable(GMM GMM.cpp)
target_link_libraries(GMM ${RUNTIME_LIBS})

add_executable(GMM2 GMM-2.cpp)
target_link_libraries(GMM2 ${RUNTIME_LIBS})


add_executable(MiniONN MiniONN.cpp)
target_link_libraries(MiniONN ${RUNTIME_LIBS})

add_executable(FBPacking FBPacking.cpp)
target_link_libraries(FBPacking ${RUNTIME_LIBS})

add_executable(FindParams FindParams.cpp)
target_link_libraries(FindParams ${RUNTIME_LIBS})

add_executable(GetCtxtSize GetCtxtSize.cpp)
target_link_libraries(GetCtxtSize ${RUNTIME_LIBS})

