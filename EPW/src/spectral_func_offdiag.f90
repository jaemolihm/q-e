  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE spectral_func_q_offdiag (iqq, iq, totq, cufkk_all)
  !-----------------------------------------------------------------------
  !!  jmlim: include the off-diagonal part of self-energy to calculate the
  !!         spectral function
  !!
  !!  Compute the electron spectral function including the  electron-
  !!  phonon interaction in the Migdal approximation. 
  !!  
  !!  We take the trace of the spectral function to simulate the photoemission
  !!  intensity. I do not consider the c-axis average for the time being.
  !!  The main approximation is constant dipole matrix element and diagonal
  !!  selfenergy. The diagonality can be checked numerically. 
  !!
  !!  Use matrix elements, electronic eigenvalues and phonon frequencies
  !!  from ep-wannier interpolation
  !!
  !-----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iospectral_sup ,iospectral
  USE phcom,         ONLY : nmodes
  USE epwcom,        ONLY : nbndsub, eps_acustic, &
                            fsthick, eptemp, ngaussw, degaussw, wmin_specfun,&
                            wmax_specfun, nw_specfun, shortrange, &
                            efermi_read, fermi_energy
  USE pwcom,         ONLY : nelec, ef
  USE klist_epw,     ONLY : isk_dummy
  USE elph2,         ONLY : etf, ibndmin, ibndmax, nkqf, xqf, &
                            epf17, wkf, nkf, wf, wqf, xkf, nkqtotf,&
                            esigmar_all_offd, esigmai_all_offd, a_all_offd, &
                            esigmar_all, esigmai_all
  USE constants_epw, ONLY : ryd2mev, one, ryd2ev, two, zero, pi, ci, eps8, czero, cone
  USE mp,            ONLY : mp_barrier, mp_sum
  USE mp_global,     ONLY : me_pool, inter_pool_comm
  USE division,      ONLY : fkbounds
  !
  implicit none
  !
  INTEGER, INTENT (in) :: iqq
  !! Current q-point index in selecq
  INTEGER, INTENT (in) :: iq
  !! Current q-point index
  INTEGER, INTENT (in) :: totq
  !! Total number of q-point in window
  complex(DP) :: cufkk_all(nbndsub, nbndsub, nkf)
  !! jmlim: electron eigenvector. To rotate self-energy to Wannier basis.
  !
  ! Local variables
  !
  INTEGER :: iw
  !! Counter on the frequency
  INTEGER :: ik
  !! Counter on the k-point index
  INTEGER :: ikk
  !! k-point index
  INTEGER :: ikq
  !! q-point index
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: imode
  !! Counter on mode
  INTEGER :: fermicount
  !! Number of states on the Fermi surface
  INTEGER :: nksqtotf
  !! Total number of k+q points
  INTEGER :: lower_bnd
  !! Lower bounds index after k or q paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k or q paral
  !
  REAL(kind=DP) :: g2
  !! Electron-phonon matrix elements squared in Ry^2
  REAL(kind=DP) :: ekk
  !! Eigen energy on the fine grid relative to the Fermi level
  REAL(kind=DP) :: ekq
  !! Eigen energy of k+q on the fine grid relative to the Fermi level
  REAL(kind=DP) :: wq
  !! Phonon frequency on the fine grid
  REAL(kind=DP) :: ef0
  !! Fermi energy level
  REAL(kind=DP) :: wgq
  !! Bose occupation factor $n_{q\nu}(T)$
  REAL(kind=DP) :: wgkq
  !! Fermi-Dirac occupation factor $f_{nk+q}(T)$
  REAL(kind=DP) :: weight
  !! Self-energy factor
  !!$$ N_q \Re( \frac{f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} + \omega_{q\nu} - i\delta }) $$
  !!$$ + N_q \Re( \frac{1- f_{mk+q}(T) + n_{q\nu}(T)}{ \varepsilon_{nk} - \varepsilon_{mk+q} - \omega_{q\nu} - i\delta }) $$
  REAL(kind=DP) :: inv_wq
  !! $frac{1}{2\omega_{q\nu}}$ defined for efficiency reasons
  REAL(kind=DP) :: inv_eptemp0
  !! Inverse of temperature define for efficiency reasons
  REAL(kind=DP) :: g2_tmp
  !! If the phonon frequency is too small discart g
  REAL(kind=DP) :: inv_degaussw
  !! Inverse of the smearing for efficiency reasons
  REAL(kind=DP) :: ww
  !! Current frequency
  REAL(kind=DP) :: dw
  !! Frequency intervals
  real(kind=DP) :: specfun_sum, esigmar0
  real(kind=DP) :: fermi(nw_specfun)
  real(kind=DP), external :: efermig, dos_ef, wgauss
  !
  ! variables for collecting data from all pools in parallel case
  !
  real(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:)
  ! begin jmlim off-diagonal
  INTEGER :: ibnd2
  complex(DP), allocatable :: mat_temp(:,:) ! for matrix inversion
  complex(DP), allocatable :: green_diag(:,:) ! for matrix inversion
  complex(DP), allocatable :: mat_temp2(:,:) ! for unitary rotation
  complex(DP), allocatable :: mat_temp3(:,:) ! for unitary rotation
  REAL(DP), allocatable :: evals_temp(:)
  complex(DP), allocatable :: work(:)
  integer :: lwork
  integer, allocatable :: ipiv(:)
  integer :: info, reclen
  REAL(DP), allocatable :: a_all_offd_resolv(:,:,:)
  REAL(DP), allocatable :: a_all_d_resolv(:,:,:)
  COMPLEX(DP), allocatable :: esigma_part(:,:,:,:)
  COMPLEX(DP), allocatable :: cufkk_all_all(:,:,:)
  COMPLEX(DP), allocatable :: ham_wannier(:,:)
  COMPLEX(kind=DP), ALLOCATABLE :: weight_arr(:,:)
  ! end jmlim off-diagonal
  !
  ! SP: Define the inverse so that we can efficiently multiply instead of
  ! dividing
  !
  inv_eptemp0 = 1.0/eptemp
  inv_degaussw = 1.0/degaussw
  !
  ! energy range and spacing for spectral function
  !
  dw = ( wmax_specfun - wmin_specfun ) / dble (nw_specfun-1)
  !
  IF (iqq == 1) THEN
    !
    WRITE(stdout,'(/5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Electron Spectral Function in the Migdal Approximation")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !
    IF ( fsthick < 1.d3 ) &
       WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick * ryd2ev, ' eV'
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp * ryd2ev, ' eV'
    !
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  IF (efermi_read) THEN
    !
    ef0 = fermi_energy
    !
  ELSE
    !
    ef0 = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk_dummy)
    ! if some bands are skipped (nbndskip /= 0), nelec has already been recalculated
    ! in ephwann_shuffle
    !
  ENDIF
  !
  IF (iq == 1) THEN
    WRITE (stdout, 100) degaussw * ryd2ev, ngaussw
    WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! The total number of k points
  !
  nksqtotf = nkqtotf/2 ! odd-even for k,k+q
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds( nksqtotf, lower_bnd, upper_bnd )
  !
  ! SP: Sum rule added to conserve the number of electron.
  IF (iq == 1) THEN
    WRITE (stdout,'(5x,a)') 'The sum rule to conserve the number of electron is enforced.'
    WRITE (stdout,'(5x,a)') 'The self energy is rescaled so that its real part is zero at the Fermi level.'
    WRITE (stdout,'(5x,a)') 'The sum rule replace the explicit calculation of the Debye-Waller term.'
    WRITE (stdout,'(a)') ' '
  ENDIF
  !
  ! loop over all k points of the fine mesh
  !
  ! jmlim optimization
  ALLOCATE(weight_arr(nw_specfun, ibndmax-ibndmin+1))
  !
  allocate(esigma_part(ibndmax-ibndmin+1, ibndmax-ibndmin+1, nkf, nw_specfun))
  esigma_part = czero
  !
  fermicount = 0
  DO ik=1, nkf
!    write(stdout, *) 'ik, iq', ik, iq
    !
    ikk = 2 * ik - 1
    ikq = ikk + 1
    !
    ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
    ! (but in this case they are the same)
    !
    IF ((MINVAL(ABS(etf (:, ikk) - ef) ) < fsthick) .AND. &
        (MINVAL(ABS(etf (:, ikq) - ef) ) < fsthick)) THEN
      !
      fermicount = fermicount + 1
      DO imode=1, nmodes
        !
        ! the phonon frequency and Bose occupation
        wq = wf (imode, iq)
        ! SP: Define the inverse for efficiency
        inv_wq = 1.0/( two * wq )
        wgq = wgauss( -wq*inv_eptemp0, -99)
        wgq = wgq / ( one - two * wgq )
        !
        ! SP: Avoid if statement in inner loops
        IF (wq > eps_acustic) THEN
          g2_tmp = 1.0
        ELSE
          g2_tmp = 0.0
        ENDIF
        !
        ! pre-compute weight, which are independent of ibnd
        DO jbnd = 1, ibndmax-ibndmin+1
          !
          !  the fermi occupation for k+q
          ekq = etf (ibndmin-1+jbnd, ikq) - ef0
          wgkq = wgauss( -ekq/eptemp, -99)
          !
          DO iw = 1, nw_specfun
            !
            ww = wmin_specfun + dble (iw-1) * dw
            !
            ! jml: real and imaginary done at once
            weight_arr(iw, jbnd) = wqf(iq) * (                                   &
              ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) - ci * degaussw ) ) )
            !
          ENDDO !iw
          !
        ENDDO !jbnd
        !
        ! end pre-computing weight
        !
        DO ibnd = 1, ibndmax-ibndmin+1
          DO ibnd2 = 1, ibndmax-ibndmin+1 ! off-diagonal
            !
            !  the energy of the electron at k (relative to Ef)
            ekk = etf (ibndmin-1+ibnd, ikk) - ef0
            !
            DO jbnd = 1, ibndmax-ibndmin+1
              !
              !  the fermi occupation for k+q
              ekq = etf (ibndmin-1+jbnd, ikq) - ef0
              wgkq = wgauss( -ekq/eptemp, -99)
              !
              ! here we take into account the zero-point sqrt(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              IF ( shortrange .AND. ( abs(xqf (1, iq))> eps8 .OR. abs(xqf (2, iq))> eps8 &
                 .OR. abs(xqf (3, iq))> eps8 )) THEN
                ! SP: The abs has to be removed. Indeed the epf17 can be a pure imaginary
                !     number, in which case its square will be a negative number.
                ! g2 = REAL( (epf17 (jbnd, ibnd, imode, ik)**two)*inv_wq*g2_tmp )
                g2 = REAL(epf17 (jbnd, ibnd, imode, ik) &
                        * epf17 (jbnd, ibnd2, imode, ik), DP) * inv_wq * g2_tmp
                ! jmlim: is this correct?
              ELSE
                ! g2 = (abs(epf17 (jbnd, ibnd, imode, ik))**two)*inv_wq*g2_tmp
                g2 = CONJG(epf17 (jbnd, ibnd, imode, ik)) &
                   * epf17 (jbnd, ibnd2, imode, ik) * inv_wq * g2_tmp
              ENDIF
              !
              DO iw = 1, nw_specfun
                !
                esigma_part(ibnd,ibnd2,ik,iw) = esigma_part(ibnd,ibnd2,ik,iw) &
                                              + g2 * weight_arr(iw, jbnd)
                !
  ! jml: ignore this correction (no Debye-Waller term in my model, anyway.)
  !              ! SP : Application of the sum rule
  !              esigmar0 =  g2 *  wqf(iq) * real (                                   &
  !                ( (       wgkq + wgq ) / ( -( ekq - wq ) - ci * degaussw )  +  &
  !                  ( one - wgkq + wgq ) / ( -( ekq + wq ) - ci * degaussw ) ) )
  !              esigmar_all_offd(ibnd,ik+lower_bnd-1,iw)=esigmar_all_offd(ibnd,ik+lower_bnd-1,iw)-esigmar0
                !
              ENDDO
              !
            ENDDO !jbnd
            !
          ENDDO ! ibnd2
        ENDDO !ibnd
        !
      ENDDO !imode
      !
    ENDIF ! endif  fsthick
    !
  ENDDO ! end loop on k
  ! jmlim: rotate self-energy from energy basis to Wannier basis
  ! This is needed to add self-energy at different iq values
  ! esigmar_all_offd(ibnd,ibnd2,ik+lower_bnd-1,iw) = esigmar_all_offd(ibnd,ibnd2,ik+lower_bnd-1,iw) + g2 * weight
  if (ibndmax-ibndmin+1 /= nbndsub) &
    CALL errore('It is assumed that ibndmin-ibndmin+1 == nbndsub ', ibndmax - ibndmin+ 1)
  ALLOCATE(mat_temp(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
  ALLOCATE(mat_temp2(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
  DO ik = 1, nkf
    DO iw = 1, nw_specfun
      CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
          cone, cufkk_all(1,1,ik), nbndsub, esigma_part(1,1,ik,iw), nbndsub, &
          czero, mat_temp, nbndsub)
      CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
          cone, mat_temp, nbndsub, cufkk_all(1,1,ik), nbndsub, &
          czero, mat_temp2, nbndsub)
      esigmar_all_offd(:,:,ik+lower_bnd-1,iw) = esigmar_all_offd(:,:,ik+lower_bnd-1,iw) +  REAL( mat_temp2 )
      esigmai_all_offd(:,:,ik+lower_bnd-1,iw) = esigmai_all_offd(:,:,ik+lower_bnd-1,iw) + AIMAG( mat_temp2 )
    END DO ! iw
  END DO ! ik
  DEALLOCATE(mat_temp)
  DEALLOCATE(mat_temp2)
  deallocate(esigma_part)
  !
  ! The k points are distributed among pools: here we collect them
  !
  IF (iqq == totq) THEN
    !
    ! jmlim: gather electron wavefunction
    allocate(cufkk_all_all(nbndsub, nbndsub, nksqtotf))
    cufkk_all_all = czero
    do ik = 1, nkf
      cufkk_all_all(:,:,ik+lower_bnd-1) = cufkk_all(:,:,ik)
    end do
    !
    !
    ALLOCATE (xkf_all(3,       nkqtotf))
    ALLOCATE (etf_all(nbndsub, nkqtotf))
    xkf_all(:,:) = zero
    etf_all(:,:) = zero
    !
#if defined(__MPI)
    !
    ! note that poolgather2 works with the doubled grid (k and k+q)
    !
    CALL poolgather2 ( 3,       nkqtotf, nkqf, xkf,    xkf_all  )
    CALL poolgather2 ( nbndsub, nkqtotf, nkqf, etf,    etf_all  )
    CALL mp_sum( esigmar_all_offd, inter_pool_comm )
    CALL mp_sum( esigmai_all_offd, inter_pool_comm )
    CALL mp_sum( fermicount, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    call mp_sum( cufkk_all_all, inter_pool_comm ) ! jmlim
    !
#else
    !
    xkf_all = xkf
    etf_all = etf
    !
#endif
    !
    ! Output electron spectral function here after looping over all q-points 
    ! (with their contributions summed in a etc.)
    !
    WRITE(stdout,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")')
    !
    ! construct the trace of the spectral function (assume diagonal selfenergy
    ! and constant matrix elements for dipole transitions)
    !
!     IF (me_pool == 0) then
!       OPEN(UNIT=iospectral,FILE='specfun_offdiag.elself')
! !      OPEN(UNIT=iospectral_sup,FILE='specfun_sup.elself')
!     ENDIF
!     IF (me_pool == 0) then
!       WRITE(iospectral, '(/2x,a/)') '#Electronic spectral function (meV)'
!       WRITE(iospectral, '(/2x,a/)') '# jmlim: Include off-diagonal part of the self-energy'
! !      WRITE(iospectral_sup, '(/2x,a/)') '#KS eigenenergies + real and im part of electronic self-energy (meV)'
!     ENDIF
!     IF (me_pool == 0) then
!       WRITE(iospectral, '(/2x,a/)') '#K-point    Energy[meV]     A(k,w)[meV^-1]'
! !      WRITE(iospectral_sup, '(/2x,a/)') '#K-point    Band   e_nk[eV]   w[eV]   &
! !&         Real Sigma[meV]  Im Sigma[meV]'
!     ENDIF
    !
    ALLOCATE(mat_temp(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
    ALLOCATE(mat_temp2(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
    ALLOCATE(mat_temp3(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
    ALLOCATE(green_diag(ibndmax-ibndmin+1, ibndmax-ibndmin+1))
    ALLOCATE(evals_temp(ibndmax-ibndmin+1))
    ALLOCATE(ipiv(ibndmax-ibndmin+1))
    ALLOCATE(a_all_offd_resolv(ibndmax-ibndmin+1, nw_specfun, nksqtotf))
    ALLOCATE(a_all_d_resolv(ibndmax-ibndmin+1, nw_specfun, nksqtotf))
    ALLOCATE(ham_wannier(nbndsub, nbndsub))
    a_all_offd_resolv = 0.d0
    a_all_d_resolv = 0.d0
    DO ik=1, nksqtotf
      !
      ikk = 2 * ik - 1
      ikq = ikk + 1
      !
! jml: do not write in stdout
!      WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f12.7)') ik, xkf_all (:,ikk)
!      WRITE(stdout,'(5x,a)') repeat('-',67)
      !
      ! ham_wannier = cufkk_all_all(:,:,ik).H @ diag(ekk) @ cufkk_all_all(:,:,ik)
      ham_wannier = czero
      do ibnd = 1, nbndsub
        ham_wannier(ibnd, ibnd) = CMPLX(etf_all (ibnd, ikk) - ef0, zero)
      end do
      CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
          cone, cufkk_all_all(1,1,ik), nbndsub, ham_wannier, nbndsub, &
          czero, mat_temp, nbndsub)
      CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
          cone, mat_temp, nbndsub, cufkk_all_all(1,1,ik), nbndsub, &
          czero, ham_wannier, nbndsub)
      !
      DO iw=1, nw_specfun
        !
        ww = wmin_specfun + dble (iw-1) * dw
        !
        mat_temp = - ham_wannier
        DO ibnd=1, ibndmax-ibndmin+1
          mat_temp(ibnd, ibnd) = mat_temp(ibnd, ibnd) + ww
        ENDDO

        ! jml: diagonal self-energy only (Note: diagonal in energy eigenbasis)
        ! green_diag = ww - ham_wannier - cufkk_all_all.H @ diag(esigma_all) @ cufkk_all_all
        green_diag = mat_temp
        mat_temp2 = czero
        do ibnd = 1, nbndsub
          mat_temp2(ibnd, ibnd) = CMPLX(esigmar_all(ibnd,ik,iw), esigmai_all(ibnd,ik,iw))
        end do
        CALL ZGEMM('C', 'N', nbndsub, nbndsub, nbndsub, &
          cone, cufkk_all_all(1,1,ik), nbndsub, mat_temp2, nbndsub, &
          czero, mat_temp3, nbndsub)
        CALL ZGEMM('N', 'N', nbndsub, nbndsub, nbndsub, &
          -cone, mat_temp3, nbndsub, cufkk_all_all(1,1,ik), nbndsub, &
          cone, green_diag, nbndsub)


        ! jml: Note that esigma are in Wannier basis, not energy eigenbasis
        !
        ! jml: does this work?
        ! mat_temp(:,:) = mat_temp(:,:) - esigmar_all_offd(:,:,ik,iw) - ci * esigmai_all_offd(:,:,ik,iw)
        !
        !
        DO ibnd=1, ibndmax-ibndmin+1
          DO ibnd2=1, ibndmax-ibndmin+1
            ! if (iw == 1) then
            !   write(stdout, *) ibnd, ibnd2, esigmar_all_offd(ibnd, ibnd2, ik, iw), esigmai_all_offd(ibnd, ibnd2, ik, iw)
            ! end if

            mat_temp(ibnd,ibnd2) = mat_temp(ibnd,ibnd2) - CMPLX(esigmar_all_offd(ibnd,ibnd2,ik,iw), esigmai_all_offd(ibnd,ibnd2,ik,iw))
          ENDDO
        ENDDO


        CALL ZGETRF(ibndmax-ibndmin+1, ibndmax-ibndmin+1, mat_temp, ibndmax-ibndmin+1, ipiv, info)
        allocate(work(1))
        CALL ZGETRI(ibndmax-ibndmin+1, mat_temp, ibndmax-ibndmin+1, ipiv, work, -1, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        CALL ZGETRI(ibndmax-ibndmin+1, mat_temp, ibndmax-ibndmin+1, ipiv, work, lwork, info)
        deallocate(work)

        CALL ZGETRF(ibndmax-ibndmin+1, ibndmax-ibndmin+1, green_diag, ibndmax-ibndmin+1, ipiv, info)
        allocate(work(1))
        CALL ZGETRI(ibndmax-ibndmin+1, green_diag, ibndmax-ibndmin+1, ipiv, work, -1, info)
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        CALL ZGETRI(ibndmax-ibndmin+1, green_diag, ibndmax-ibndmin+1, ipiv, work, lwork, info)
        deallocate(work)

        DO ibnd=1, ibndmax-ibndmin+1
          a_all_offd(iw,ik) = a_all_offd(iw,ik) + aimag(mat_temp(ibnd,ibnd)) / pi
          ! Wannier-resolved spectral function
          a_all_offd_resolv(ibnd,iw,ik) = aimag(mat_temp(ibnd,ibnd)) / pi
          a_all_d_resolv(ibnd,iw,ik) = aimag(green_diag(ibnd,ibnd)) / pi
        END DO

        !
! jml: do not write in stdout
!        WRITE(stdout, 103) ik, ryd2ev * ww, a_all_offd(iw,ik) / ryd2mev
        !
      ENDDO ! iw
      !
! jml: do not write in stdout
!      WRITE(stdout,'(5x,a/)') repeat('-',67)
      !
    ENDDO ! ik
    DEALLOCATE(mat_temp)
    DEALLOCATE(mat_temp2)
    DEALLOCATE(mat_temp3)
    DEALLOCATE(ham_wannier)
    DEALLOCATE(cufkk_all_all)
    ! !
    ! DO ik=1, nksqtotf
    !   !
    !   ! The spectral function should integrate to 1 for each k-point
    !   specfun_sum = 0.0
    !   !
    !   DO iw=1, nw_specfun
    !     !
    !     ww = wmin_specfun + dble (iw-1) * dw
    !     fermi(iw) = wgauss(-ww/eptemp, -99)
    !     !WRITE(stdout,'(2x,i7,2x,f12.4,2x,e12.5)') ik, ryd2ev * ww, a_all_offd(iw,ik) / ryd2mev
    !     !
    !     specfun_sum = specfun_sum + a_all_offd(iw,ik)* fermi(iw) * dw !/ ryd2mev
    !     !
    !   IF (me_pool == 0) &
    !     WRITE(iospectral,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all_offd(iw,ik) / ryd2mev
    !     !
    !   ENDDO
    !   !
    !   IF (me_pool == 0) &
    !     WRITE(iospectral,'(a)') ' '
    !   IF (me_pool == 0) &
    !     WRITE(iospectral,'(2x,a,2x,e12.5)') '# Integrated spectral function ',specfun_sum
    !   !
    ! ENDDO
    ! !
    ! IF (me_pool == 0)  CLOSE(iospectral)
    ! !
    if (me_pool == 0) then
      inquire(iolength=reclen) a_all_offd(:,:)
      open(UNIT=iospectral, FILE='specfun_offdiag.elself', form='unformatted', &
           status='unknown', access='direct', recl=reclen)
      write(iospectral, rec=1) a_all_offd / ryd2mev
      close(iospectral)

      inquire(iolength=reclen) a_all_offd_resolv(:,:,:)
      open(UNIT=iospectral, FILE='specfun_offdiag_resolv.elself', form='unformatted', &
           status='unknown', access='direct', recl=reclen)
      write(iospectral, rec=1) a_all_offd_resolv / ryd2mev
      close(iospectral)

      inquire(iolength=reclen) a_all_d_resolv(:,:,:)
      open(UNIT=iospectral, FILE='specfun_diag_resolv.elself', form='unformatted', &
           status='unknown', access='direct', recl=reclen)
      write(iospectral, rec=1) a_all_d_resolv / ryd2mev
      close(iospectral)

      ! inquire(iolength=reclen) ryd2mev * (esigmar_all_offd + ci * esigmai_all_offd)
      ! open(UNIT=iospectral, FILE='esigma.bin', form='unformatted', &
      !      status='unknown', access='direct', recl=reclen)
      ! write(iospectral, rec=1) ryd2mev * (esigmar_all_offd + ci * esigmai_all_offd)
      ! close(iospectral)
    END if

!     DO ibnd=1, ibndmax-ibndmin+1
!       !
!       DO ik=1, nksqtotf
!         !
!         ikk = 2 * ik - 1
!         ikq = ikk + 1
!         !
!         !  the energy of the electron at k
!         ekk = etf_all (ibndmin-1+ibnd, ikk) - ef0
!         !
!         DO iw = 1, nw_specfun
!           !
!           ww = wmin_specfun + dble (iw-1) * dw
! ! jml: do not write in stdout
! !          WRITE(stdout,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
! !            ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all_offd(ibnd,ik,iw),&
! !            ryd2mev * esigmai_all_offd(ibnd,ik,iw)
!           !
! !          IF (me_pool == 0) &
! !          WRITE(iospectral_sup,'(2i9,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4,2x,f12.4)') ik,&
! !            ibndmin-1+ibnd, ryd2ev * ekk, ryd2ev * ww, ryd2mev * esigmar_all_offd(ibnd,ik,iw),&
! !            ryd2mev * esigmai_all_offd(ibnd,ik,iw)
!           !
!         ENDDO
!         !
!       ENDDO
!       !
!       WRITE(stdout,*) ' '
!       !
!     ENDDO
    !
!    IF (me_pool == 0)  CLOSE(iospectral_sup)
    !
    DEALLOCATE (xkf_all)
    DEALLOCATE (etf_all)
    !
  ENDIF
  !
  100 FORMAT(5x,'Gaussian Broadening: ',f10.6,' eV, ngauss=',i4)
  103 FORMAT(5x,'ik = ',i7,'  w = ',f9.4,' eV   A(k,w) = ',e12.5,' meV^-1')
  !
  RETURN
  !
  END SUBROUTINE spectral_func_q_offdiag
  !
