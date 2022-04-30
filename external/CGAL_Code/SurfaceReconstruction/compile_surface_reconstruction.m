function compile_surface_reconstruction
    
    mex -v surface_reconstruction.cpp ...
        CXXOPTIMFLAGS="-O3" ...
        CXXFLAGS="$CXXFLAGS -march=native -fopenmp -fPIC" ...
        LDFLAGS="$LDFLAGS -fopenmp" ...
        -I/usr/include ...
        -I/usr/local/include ...
        -I/usr/local/include/libigl/include ...
        -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lmpfr -lboost_thread -lboost_system
    
end