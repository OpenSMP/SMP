SET(SRC
    Matrix.cpp
    literal.cpp
    HElib.cpp
    )

add_library(cryptgmm ${SRC})

set(RUNTIME_LIBS   
    cryptgmm 
    helib
    ntl
    pthread
    boost_system
    gmp)

add_executable(GMM GMM.cpp)
target_link_libraries(GMM ${RUNTIME_LIBS})

add_executable(MiniONN MiniONN.cpp)
target_link_libraries(MiniONN ${RUNTIME_LIBS})

add_executable(FindParams FindParams.cpp)
target_link_libraries(FindParams ${RUNTIME_LIBS})
