program run_wigner
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, bg
    USE constants_epw, ONLY : eps6
    USE wigner,        ONLY : wigner_seitz_wrap
    implicit none

    integer :: nk1, nk2, nk3, nq1, nq2, nq3
    integer :: nbndsub, nat
    integer :: i, j, k

  INTEGER :: dims
  !! Dims is either nbndsub if use_ws or 1 if not
  INTEGER :: dims2
  !! Dims is either nat if use_ws or 1 if not
  INTEGER, ALLOCATABLE :: irvec_k(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point in the basis
  !! of the lattice vectors for electrons
  INTEGER, ALLOCATABLE :: irvec_q(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point for phonons
  INTEGER, ALLOCATABLE :: irvec_g(:,:)
  !! integer components of the ir-th Wigner-Seitz grid point for electron-phonon
  INTEGER, ALLOCATABLE :: ndegen_k (:,:,:)
  !! Wigner-Seitz number of degenerescence (weights) for the electrons grid
  INTEGER, ALLOCATABLE :: ndegen_q (:,:,:)
  !! Wigner-Seitz weights for the phonon grid that depend on 
  !! atomic positions $R + \tau(nb) - \tau(na)$
  INTEGER, ALLOCATABLE :: ndegen_g (:,:,:,:)
  !! Wigner-Seitz weights for the electron-phonon grid that depend on 
  !! atomic positions $R - \tau(na)$
  INTEGER :: nrr_k 
  !! Number of WS points for electrons
  INTEGER :: nrr_q
  !! Number of WS points for phonons
  INTEGER :: nrr_g
  !! Number of WS points for electron-phonons
  REAL(kind=DP), ALLOCATABLE :: wslen_k(:)
  !! real-space length for electrons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: wslen_q(:)
  !! real-space length for phonons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: wslen_g(:)
  !! real-space length for electron-phonons, in units of alat
  REAL(kind=DP), ALLOCATABLE :: w_centers(:,:)
  !! Wannier centers  
  REAL(DP), ALLOCATABLE :: tau(:,:)     !  initial positions read from stdin (in bohr)
  REAL(DP) :: dummy(3)
  LOGICAL :: use_ws

   

  use_ws = .false.
  nk1 = 3
  nk2 = 3
  nk3 = 1
  nq1 = 3
  nq2 = 3
  nq3 = 1
  nbndsub = 1
  nat = 1
  ALLOCATE(w_centers(3,nbndsub))
  ALLOCATE(tau(3,nat))
  w_centers = 0.0_dp
  tau = 0.0_dp
  
! simple cubic
  at(:,1) = (/ 1.0, 0.0, 0.0/)
  at(:,2) = (/ 0.0, 1.0, 0.0/)
  at(:,3) = (/ 0.0, 0.0, 1.0/)
  bg(:,1) = (/ 1.0, 0.0, 0.0/)
  bg(:,2) = (/ 0.0, 1.0, 0.0/)
  bg(:,3) = (/ 0.0, 0.0, 1.0/)

! fcc
!  at(:,1) = (/-1.0, 0.0, 1.0/) / 2.0
!  at(:,2) = (/ 0.0, 1.0, 1.0/) / 2.0
!  at(:,3) = (/-1.0, 1.0, 0.0/) / 2.0
!  bg(:,1) = (/-1.0,-1.0, 1.0/)
!  bg(:,2) = (/ 1.0, 1.0, 1.0/)
!  bg(:,3) = (/-1.0, 1.0,-1.0/)


! input: nk1, nk2, nk3, nq1, nq2, nq3, dims, dims2, w_centers, tau
  IF (use_ws) THEN
    ! Use Wannier-centers to contstruct the WS for electonic part and el-ph part
    ! Use atomic position to contstruct the WS for the phonon part
    dims  = nbndsub
    dims2 = nat
    CALL wigner_seitz_wrap ( nk1, nk2, nk3, nq1, nq2, nq3, irvec_k, irvec_q, irvec_g, &
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                             w_centers, dims, tau, dims2 )
  ELSE
    ! Center the WS at Gamma for electonic part, the phonon part and el-ph part
    dims  = 1
    dims2 = 1
    dummy(:) = (/0.0,0.0,0.0/)
    CALL wigner_seitz_wrap ( nk1, nk2, nk3, nq1, nq2, nq3, irvec_k, irvec_q, irvec_g, &
                             ndegen_k, ndegen_q, ndegen_g, wslen_k, wslen_q, wslen_g, &
                             dummy, dims, dummy, dims2 )
  ENDIF
  ! 
  ! Determine the size of the respective WS sets based on the length of the matrices
  nrr_k = SIZE(irvec_k(1,:))
  nrr_q = SIZE(irvec_q(1,:))
  nrr_g = SIZE(irvec_g(1,:))


  write(*,*) '========================================================'
  write(*,'(x,a,I7)') 'Wigner-Seitz cell for ELECTRON: nrr_k = ', nrr_k
  write(*,'(x,a)') '(i_R-1)       irvec    ndegen       wslen'
  do i = 1, nrr_k
    write(*,'(x,I7, 3I5, I7, 1F12.7)') i-1, irvec_k(:,i), ndegen_k(i,1,1), wslen_k(i)
  end do
  write(*,*) '========================================================'
  write(*,'(x,a,I7)') 'Wigner-Seitz cell for PHONON: nrr_q = ', nrr_q
  write(*,'(x,a)') '(i_R-1)       irvec    ndegen       wslen'
  do i = 1, nrr_k
    write(*,'(x,I7, 3I5, I7, 1F12.7)') i-1, irvec_q(:,i), ndegen_q(i,1,1), wslen_q(i)
  end do
  write(*,*) '========================================================'
  write(*,'(x,a,I7)') 'Wigner-Seitz cell for ELECTRON-PHONON: nrr_g = ', nrr_g
  write(*,'(x,a)') '(i_R-1)       irvec    ndegen       wslen'
  do i = 1, nrr_k
    write(*,'(x,I7, 3I5, I7, 1F12.7)') i-1, irvec_g(:,i), ndegen_g(i,1,1,1), wslen_g(i)
  end do
end program
