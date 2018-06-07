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
subroutine RaoultPressure(LNGAMMA,SatPressure,TotalPressure,FRACV)
! This routine computes bubble pressure and composition according to the modified
! Raoult's law.
! It uses lngamma computed before and pure compounds vapour pressures for a given
! temperature.
!
    use variables
    implicit none

    real, dimension(COMP), intent(in) :: LNGAMMA, SatPressure
    real, dimension(COMP), intent(out)  :: FRACV ! Component's mass fraction on vapour fase
    real, intent(out) :: TotalPressure ! System's total pressure and
    real, dimension(COMP) :: GAMMACOMP

    do I = 1, COMP
        GAMMACOMP(I) = exp(LNGAMMA(I))
    enddo

    TotalPressure = 0.0 ! Starts TotalPressure
    FRACV = 0.0 ! Starts FRACV1

    ! Computes the system's total Pressure - modified Raoult's law
    do I = 1, COMP
        TotalPressure = TotalPressure + FRAC(I)*GAMMACOMP(I)*SatPressure(I)
    enddo

    ! Computes the compound 1 vapour mass fraction and registers it
    do I = 1, COMP
        FRACV(I) = FRAC(I)*GAMMACOMP(I)*SatPressure(I)/TotalPressure
    enddo

end subroutine RaoultPressure
