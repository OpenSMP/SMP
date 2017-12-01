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
