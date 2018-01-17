## Secure Matrix Multiplication Protocol
A two-party secure computation protocol for computing matrix multiplication.

* Use HElib as the basic FHE implementation.
    * Install NTL, GMP
	* plz use the fork version in the submodule.
* Use Boost network library

* Build
	* pull the submodule, `git submodule init && git submodule update`
	* create the build dir, `mkdir build && cd build`
	* cmake and make `cmake -DCMAKE_BUILD_TYPE=Release && make`
	* check out the executable in `build/bin`

* Main Executable
	* `bin/SMP` with parameters `R, a, N, M` and `D`
	   ```
	   Usage: bin/SMP [ name=value ]...
        N       n1  [ default=8 ]
        M       n2  [ default=8 ]
        D       n3  [ default=8 ]
        R       role. 0 for server and 1 for client  [ default=-1 ]
        a       server address  [ default=127.0.0.1 ]
        p       port  [ default=12345 ]

	   ```
	* The client holds a `N*M` matrix, and the server holds a `M*D` matrix. 
	* The computed `N*D` matrix will be sent back to the client.
