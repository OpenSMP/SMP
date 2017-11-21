SET(SRC
    DoublePacking.cpp
    Matrix.cpp
    literal.cpp)

add_library(cryptgmm ${SRC})

add_executable(main GMM.cpp)
target_link_libraries(main cryptgmm ntl helib-b pthread gmp boost_system)
