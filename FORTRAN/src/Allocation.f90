! Copyright (c) 2011-2013, Federal University of Rio Grande do Sul
! All rights reserved.
!
! This software is subject to a BSD License, please see the License.txt
! file for more information.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED.
!
! Authors: Luiz Felipe Kusler Possani
!          Rafael de Pelegrini Soares
!
! This is a demonstration code, for more efficient implementations
! please contact rafael@enq.ufrgs.br.
!
! PLEASE CITE AS
!   * Soares and Gerber (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie400170a;
!   * Soares et al. (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie4013979.
!

module variables
! This module contains the variables that are shared by the subroutines used to compute
! the activity coefficient according to the F-SAC model.
!
    REAL,PARAMETER:: EO = 2.395*10.0**(-4), RGAS = 0.001987, rav = 1.07, p = 0.75
    REAL,PARAMETER:: R = 66.69, Q = 50, PI = 3.14159265359

    REAL :: FPOL, ALPHA, ALPHAPRIME, SYSTEMP
    REAL :: SUMMATION
    REAL :: AEFFPRIME, DENOM

    INTEGER :: I, J, K, L, COMPSEG, COMP, n, ii, mm, nn ! Auxiliary variables

    REAL, DIMENSION(:), ALLOCATABLE :: RCOSMO, QCOSMO, RNORM, QNORM, LNGAMMACOMB
    REAL, DIMENSION(:), ALLOCATABLE :: FRAC, NUMER
    REAL, DIMENSION(:,:), ALLOCATABLE :: CONVERG, CONPR
    REAL, DIMENSION(:,:), ALLOCATABLE :: SEGGAMMAPR, SEGGAMMAOLDPR, SEGGAMMA, SEGGAMMAOLD, PROFILEMATRIX
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: SigmaComp, SIGMA
    REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: DELTAW, deltaW_HB

    CHARACTER(25), DIMENSION(:), ALLOCATABLE :: COMPOUNDS

end module variables

module deltaW_HB_data
! This module contains variables used in "deltaW_HB_reader" and "calc_deltaW" routines.
!
    implicit none
    real, dimension(:,:), allocatable :: reader

end module deltaW_HB_data

subroutine allocation(SYSCOMP,compSize)
! This routine allocates all the vectors and matrices and sets the variables contained within module "variables".
! It must be called only once in the program and before the routine "activity".
!
! DO NOT FORGET to deallocate all vectors and matrices in the end of the main program using the command "STOP".
!
    use variables
    implicit none

    character(25), dimension(compSize), intent(in) :: SYSCOMP
    integer, intent(in) :: compSize

    ! Sets constants
    COMPSEG = 51 ! NUMBER OF INTERVALS FOR THE SIGMA PROFILE
    COMP = compSize ! Number os components
    FPOL = 1

    AEFFPRIME = PI*rav**2
    ALPHA = (0.3*AEFFPRIME**(1.5))/(EO)
    ALPHAPRIME = FPOL*ALPHA

    ALLOCATE(FRAC(COMP),RCOSMO(COMP), QCOSMO(COMP), RNORM(COMP), QNORM(COMP), LNGAMMACOMB(COMP))
    ALLOCATE(PROFILEMATRIX(21,COMP), SEGGAMMA(21,COMP), SEGGAMMAOLD(21,COMP), SigmaComp(21,4,COMP))
    ALLOCATE(SIGMA(COMPSEG,4,COMP), NUMER(COMPSEG),CONVERG(21, COMP), SEGGAMMAPR(21,COMP), SEGGAMMAOLDPR(21,COMP), CONPR(21,COMP))
    ALLOCATE(DELTAW(COMP, COMP, 21, 21), deltaW_HB(COMP, COMP, 21, 21))
    ALLOCATE(COMPOUNDS(compSize))

    COMPOUNDS = SYSCOMP

    SigmaComp = 0 ! Starts SigmaComp

    !READS INDIVIDUAL SIGMA PROFILES
    DO K = 1, COMP
        CALL Comps(SYSCOMP(K), SIGMA(:,:,K), RCOSMO(K), QCOSMO(K))
    END DO

    call setSigmaComp
    call calcDeltaW

end subroutine allocation

subroutine setSigmaComp
! This routine sets the sigmaComp matrix
!
    use variables
    implicit none

    DO I = 1, COMP
        mm = 0
        DO K = 1,COMPSEG
            IF(SIGMA(K,2,I) /= 0) THEN
                mm = mm + 1
                IF(SIGMA(K,3,I) == 0) THEN
                    SigmaComp(mm,1,I) = SIGMA(K,1,I)
                    SigmaComp(mm,2,I) = SIGMA(K,2,I)
                    SigmaComp(mm,3,I) = 0
                    SigmaComp(mm,4,I) = SIGMA(K,4,I)
                    ELSE
                    SigmaComp(mm,1,I) = SIGMA(K,1,I)
                    SigmaComp(mm,2,I) = SIGMA(K,2,I) - AEFFPRIME*SIGMA(K,3,I)
                    SigmaComp(mm,3,I) = 0
                    SigmaComp(mm,4,I) = SIGMA(K,4,I)
                    mm = mm + 1
                    SigmaComp(mm,1,I) = SIGMA(K,1,I)
                    SigmaComp(mm,2,I) = AEFFPRIME*SIGMA(K,3,I)
                    SigmaComp(mm,4,I) = SIGMA(K,4,I)
                    IF(SigmaComp(mm,1,I) < 0) THEN
                        SigmaComp(mm,3,I) = 1
                        ELSE
                        SigmaComp(mm,3,I) = 2
                    ENDIF
                ENDIF
            ENDIF
        ENDDO
    ENDDO

end subroutine setSigmaComp

subroutine calcDeltaW
! This routine computes DeltaW

    use variables
    use deltaW_HB_data
    implicit none

    DELTAW = 0
    deltaW_HB = 0
    CALL deltaW_HB_reader

    DO I = 1, COMP
        DO J = 1, COMP
            DO mm = 1, 21
                DO nn = 1, 21
                    CALL calc_deltaW_HB(reader(:,:), size(reader(:,1)), SigmaComp(mm,3,I), SigmaComp(nn,3,J), &
                    SigmaComp(mm,4,I), SigmaComp(nn,4,J), COMPOUNDS(I), COMPOUNDS(J), deltaW_HB(I,J,mm,nn))
                    DELTAW(I,J,mm,nn) = (ALPHAPRIME/2.0)*(SigmaComp(mm,1,I) + SigmaComp(nn,1,J))**2.0 - deltaW_HB(I,J,mm,nn)/2.0
                ENDDO
            ENDDO
        ENDDO
    ENDDO

end subroutine calcDeltaW

subroutine deltaW_HB_reader
! This routine reads DeltaW_HB data.
!
    use deltaW_HB_data
    implicit none

    integer :: ierror, n, i
    character :: header

    open(unit = 13, file = "lib\FSAC1-wHB.csv", status = "old", action = "read")

    n = 0
    ierror = 0
    read(13,*) header

    do
        read(13,*, iostat=ierror) header
        if(ierror /= 0) exit
        n = n + 1
    end do
    rewind(13)

    allocate(reader(n,3))

    read(13,*) header
    do i = 1, n
        read(13,*) header, reader(i,1), header, reader(i,2), reader(i,3)
    enddo

    close(13)

end subroutine deltaW_HB_reader

subroutine calc_deltaW_HB(reader, n, compI_HB, compJ_HB, compI_group, compJ_group, compI, compJ, deltaW_HB)
! This routine computes DeltaW_HB.
!
    implicit none

    integer, intent(in) :: n
    real, intent(in), dimension(n,3) :: reader
    real, intent(in) :: compI_HB, compJ_HB, compI_group, compJ_group
    character(25), intent(in) :: compI, compJ
    real, intent(out) :: deltaW_HB
    integer :: i
    real :: Acc, Donn

    if ((compI_HB == 0.0).or. (compJ_HB == 0.0).or.(compI_HB == compJ_HB)) then
        Donn = 0.0
        Acc = 0.0
        deltaW_HB = 0.0
        else
        if(compI_HB == 1.0) then
            Donn = compI_group
            Acc = compJ_group
            else
            Donn = compJ_group
            Acc = compI_group
        endif
        do i = 1, n
            if((reader(i,1) == Acc).and.(reader(i,2) == Donn)) then
                deltaW_HB = reader(i,3)
                exit
            endif
            if(i == n) then
                write(*,*) "Missing HB interaction for ", compI, "and ", compJ
                stop
            endif
        enddo
    endif

end subroutine calc_deltaW_HB
