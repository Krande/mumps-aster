!
!  This file is part of MUMPS 5.6.2, released
!  on Wed Oct 11 09:36:25 UTC 2023
!
!
!  Copyright 1991-2023 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
!  Mumps Technologies, University of Bordeaux.
!
!  This version of MUMPS is provided to you free of charge. It is
!  released under the CeCILL-C license 
!  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
!  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
!
!     This file includes various internal datastructures
!     passed through the main MUMPS structure between successive
!     phases of the solver. The main one is root information for
!     the multifrontal tree.
      TYPE DMUMPS_ROOT_STRUC
        SEQUENCE
        INTEGER(4) :: MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER(4) :: MYROW, MYCOL
        INTEGER(4) :: SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER(4) :: RHS_NLOC
        INTEGER(4) :: ROOT_SIZE, TOT_ROOT_SIZE
!       descriptor for scalapack
        INTEGER(4), DIMENSION( 9 ) :: DESCRIPTOR
        INTEGER(4) :: CNTXT_BLACS, LPIV, rootpad0
        INTEGER(4), DIMENSION(:), POINTER :: RG2L
        INTEGER(4), DIMENSION(:), POINTER :: IPIV
!       Centralized master of root
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_CNTR_MASTER_ROOT
!       Used to access Schur easily from root structure
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_POINTER
!       for try_null_space preprocessing constant only:
        DOUBLE PRECISION, DIMENSION(:), POINTER :: QR_TAU, rootpad1
!       Fwd in facto: 
!           case of scalapack root: to store RHS in 2D block cyclic
!           format compatible with root distribution
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: RHS_ROOT, rootpad2
!       for try_nullspace preprocessing constant only:
        DOUBLE PRECISION :: QR_RCOND, rootpad3
        LOGICAL(4) :: yes, gridinit_done
!       for SVD on root (#define try_null_space)
        DOUBLE PRECISION, DIMENSION(:,:), POINTER :: SVD_U, SVD_VT
!       for RR on root (#define try_null_space)
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SINGULAR_VALUES
        INTEGER(4) :: NB_SINGULAR_VALUES,rootpad4
!
      END TYPE DMUMPS_ROOT_STRUC
!     multicore
      TYPE DMUMPS_L0OMPFAC_T
         SEQUENCE
         DOUBLE PRECISION, POINTER, DIMENSION(:) :: A
         INTEGER(8) :: LA
      END TYPE DMUMPS_L0OMPFAC_T
