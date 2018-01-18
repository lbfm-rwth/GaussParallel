Contents of repo:

    read.g: Loads all necessary functions into standard GAP, assuming access to packages "IO" and "GAUSS"
    read_hpc.g : Corresponding version for HPCGAP. (My current versions of HPCGAP/ the GAUSS pkg are not compatible, hence the repo at the moment
                 uses a rather obscure looking work-around..)
    utils.g: Collection of small basic functions used in subfunctions of the algorithm
    subfunctions.g: Collection of larger subfunctions used in the Gaussian elimination alg.
    main_seq_trafo.g: Contains a version of the elimination alg. for standard GAP computing RREF and a transformation
    main_semi_par_trafo.g: Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, where the second step of the
                            algorithm runs in parallel ( using HPCGAP's task arch. )
    main_par_trafo.g: Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, where the first step of the
                            algorithm runs in parallel ( using HPCGAP's task arch. )
    main_full_par_trafo.g: Contains a version of the elimination alg. for HPCGAP computing RREF and a transformation, running completely in parallel.
            
