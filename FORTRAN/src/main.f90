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

program FSAC1
! This is an example program of the F-SAC method.
! The F-SAC model can be used to compute activity coefficients for a given composition and temperature.
! By using this model, the activity coefficient profile (binary mixture) and constant temperature VLE can be computed.
! This particular program illustrates these cases.
!
    implicit none

    real :: TEMPERATURE ! System temperature (KELVIN)
    character(25), dimension(:), allocatable :: Components ! Component names vector. Define the number of components in variable ncomp
    character :: yn
    integer ::  ncomp ! Number of components in the mixture
    real, dimension(:), allocatable :: FRACCOMP, LNGAMMAS, SatPressure, FRACV
    real :: Pressure
    real :: dx ! FRAC(1) increment
    integer :: i,n

    write(*,*) '---------------------------------------------------------------------------'
    write(*,*) 'Copyright (c) 2011-2012, Federal University of Rio Grande do Sul.'
    write(*,*) 'All rights reserved. This software is subject to a BSD License, please see the'
    write(*,*) 'License.txt file for more information.'
    write(*,*) ' '

    write(*,*) 'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"'
    write(*,*) 'AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,'
    write(*,*) 'THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE'
    write(*,*) 'ARE DISCLAIMED.'
    write(*,*) ' '
    write(*,*) 'Authors: Luiz Felipe Kusler Possani'
    write(*,*) '         Rafael de Pelegrini Soares'
    write(*,*) ' '

    write(*,*) 'This is a demonstration code, for more efficient implementations'
    write(*,*) 'please contact rafael@enq.ufrgs.br.'
    write(*,*) 'PLEASE CITE AS'
    write(*,*) ' * Soares and Gerber (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie400170a;'
    write(*,*) ' * Soares et al. (2013), Ind. Eng. Chem. Res. DOI:10.1021/ie4013979.'
    write(*,*) '---------------------------------------------------------------------------'
    write(*,*) ' '

    ! The number of intervals between 0.00 and 1.00 with dx = 0.01
    dx = 0.01
    n = nint((1.0d0-0.00)/dx) + 1

    do
        write(*,*) 'INSERT THE NUMBER OF COMPONENTS '
        read(*,'(I2)') ncomp
        if(ncomp > 1) exit
        write(*,*) 'WARNING: BAD ARGUMENT. TRY AGAIN. '
    enddo
    allocate(Components(ncomp), FRACCOMP(ncomp), LNGAMMAS(ncomp), SatPressure(ncomp), FRACV(ncomp))

    ! Mixture definition and memory allocation
    do i = 1, ncomp
        write(*,*) 'INSERT THE NAME OF COMPONENT ', i
        read(*,'(A25)') Components(i)
    enddo

    call allocation(Components,size(Components))
    LNGAMMAS = 0.0

    ! Example of ln gamma for a given comp. and temperature -------
    write(*,*) 'INSERT THE SYSTEM TEMPERATURE (K) '
    read(*,*) TEMPERATURE

    FRACCOMP = 0.0

    do i = 1, ncomp-1
        do
            write(*,*) 'INSERT COMPONENT ', i, ' MOLAR FRACTION (0-1)'
            read(*,*) FRACCOMP(i)
            if((FRACCOMP(i)>=0).and.(FRACCOMP(i)<=1).and.(sum(FRACCOMP)<=1.0)) exit
            write(*,*) 'WARNING: BAD ARGUMENT. TRY AGAIN.'
        enddo
    enddo
    FRACCOMP(ncomp) = 1.0 - sum(FRACCOMP)

    call calcSegGammaPure(TEMPERATURE)
    call activity(FRACCOMP,LNGAMMAS)

    do i = 1, ncomp
        write(*,*) Components(i)
        write(*,*) "GAMMA: ", exp(LNGAMMAS(i)), " LNGAMMA: ", LNGAMMAS(i)
    enddo
    ! --------------------------------------------------------------

    ! Lngamma profile
    if(ncomp == 2) then
        ! Opens output files
        OPEN(UNIT=14, FILE = "Gamma.dat", action = "write")
        WRITE(14,*) "TEMPERATURE", TEMPERATURE, "KELVIN"
        5  FORMAT (1X,A10,5X,A12,5X,A12,5X,A12,5X,A12)
        6  FORMAT (1X,A10,5X,A9,5X,A10,7X,A10,5X,A10)
        WRITE(14,6) "MOLEFRAC", "GAMMA1", "GAMMA2", "LNGAMMA1", "LNGAMMA2"
        WRITE(14,5) "X1", Components(1), Components(2), Components(1), Components(2)

        ! System's pressure profile ----------------
        write(*,*) 'COMPUTE PRESSURE PROFILE? (Y/N)'
        read(*,*) yn
        if(yn == 'Y') then
            ! Pure compound vapour pressures (for constant temperature VLE)
            write(*,*) 'INSERT COMPONENT 1 PURE VAPOR PRESSURE (BAR).'
            read(*,*) SatPressure(1)
            write(*,*) 'INSERT COMPONENT 2 PURE VAPOR PRESSURE (BAR).'
            read(*,*) SatPressure(2)

            OPEN(UNIT=15, FILE="PressureProfile.csv", ACTION="write")
            WRITE(15,*) "X1" , ";" ,"Y1" ,";" ,"P"
        endif

        FRACCOMP(1) = 0.0
        FRACCOMP(2) = 1.0 - FRACCOMP(1)
        do i = 1, n
            call activity(FRACCOMP,LNGAMMAS)
            if(yn == 'Y') call RaoultPressure(LNGAMMAS,SatPressure,Pressure,FRACV)
            ! Writes outputs
            write(14,*) FRACCOMP(1), exp(LNGAMMAS(1)), exp(LNGAMMAS(2)), LNGAMMAS(1), LNGAMMAS(2)
            write(15,*) FRACCOMP(1), ";" ,FRACV(1) ,";" ,Pressure
            FRACCOMP(1) = FRACCOMP(1) + dx
            FRACCOMP(2) = 1.0 - FRACCOMP(1)
        enddo
    endif

    stop

end program FSAC1
