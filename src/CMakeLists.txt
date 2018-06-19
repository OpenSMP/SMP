SET(SRC
    Matrix.cpp
    literal.cpp
    HElib.cpp
    net_io.cpp
	SMPServer.cpp
	DoublePacking.cpp
    )

add_library(cryptgmm ${SRC})

set(RUNTIME_LIBS   
    fhe
	cryptgmm
    ntl
    pthread
    boost_system
    gmp)

macro (add_exec _name)
	add_executable(${_name} "${_name}.cpp")
    target_link_libraries(${_name} ${RUNTIME_LIBS})
endmacro()

add_exec(GetCtxtSize)
add_exec(FindParams)
add_exec(MiniONN)
add_exec(SMP)
add_exec(SMP_local)
add_exec(faster_unpacking)
add_exec(fast_decryption)
add_subdirectory(DGK)
