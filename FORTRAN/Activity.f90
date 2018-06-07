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

subroutine activity(CompFrac,LNGAMMA)
! Computes the activity coefficient according to the F-SAC model for a given
! composition.
! It must be called after the routines "allocation" and "calcSegGammaPure".
!
! This routine considers the temperature given in routine "calcSegGammaPure".
! Thus, in order to update the temperature, the user should call
! the "calcSegGammaPure" routine first.
!
! Please check the example programs for usage.
! This routine contains code from the COSMO-SAC implementation
! available at http://www.design.che.vt.edu.
!
    USE variables
    IMPLICIT NONE

    real, dimension(COMP), intent(in) :: CompFrac
    real, dimension(COMP), intent(out) :: LNGAMMA
    real :: LNGAMMARES

    if(sum(CompFrac) /= 1.0) then
        write(*,*) "Warning: the components molar fraction must sum 1.0"
        stop
    endif

    FRAC = CompFrac

    !CALCULATE THE MIXTURE SIGMA PROFILE
    DENOM = sum(FRAC*QCOSMO)
    DO I = 1, COMP
        DO mm =1, 21
            PROFILEMATRIX(mm,I) = FRAC(I)*SigmaComp(mm,2,I)/DENOM
        END DO
    ENDDO

    call calcSegGammaMixture
    call calcStavermanGugenheimen

    !CALCULATION OF LNGAMMAS
    DO I = 1, COMP
        LNGAMMARES = 0.0
        DO mm = 1, 21
            LNGAMMARES = LNGAMMARES + ((SigmaComp(mm,2,I)/AEFFPRIME)*(LOG(SEGGAMMA(mm,I)/(SEGGAMMAPR(mm,I)))))
        END DO
        LNGAMMA(I) = LNGAMMARES + LNGAMMACOMB(I)
    ENDDO

end subroutine activity

subroutine calcSegGammaMixture
! This routine computes the gamma values for the compounds in mixture.
!
    USE variables
    IMPLICIT NONE

    SEGGAMMA = 1.0
    SEGGAMMAOLD = 1.0
    L = 0
    DO
        SEGGAMMAOLD = SEGGAMMA
        DO I = 1, COMP
            DO mm = 1, 21
                SUMMATION = 0.0
                do J = 1, COMP
                    DO nn = 1, 21
                        SUMMATION = SUMMATION + PROFILEMATRIX(nn,J)*SEGGAMMAOLD(nn,J)* &
                        EXP(-DELTAW(I,J,mm,nn)/(RGAS*SYSTEMP))
                    END DO
                END DO
                SEGGAMMA(mm,I) = 1.0/SUMMATION
                SEGGAMMA(mm,I) = (SEGGAMMA(mm,I)+SEGGAMMAOLD(mm,I))/2.0
                CONVERG(mm,I) = ABS((SEGGAMMA(mm,I)-SEGGAMMAOLD(mm,I))/SEGGAMMAOLD(mm,I))
            END DO
        ENDDO
        L = L+1
        IF (MAXVAL(CONVERG) <=0.000001) EXIT
    ENDDO
    WRITE(*,*) "SEGGAMMA COMP niter: ", L

end subroutine calcSegGammaMixture

subroutine calcSegGammaPure(Temp)
! This routine computes the gamma values for the pure compounds in the given temperature.
! It must be called after "allocation" and before "activity".
!
! If the system temperature is constant between activity computations, this routine can
! be called only once in the program. Every time the temperature changes, this routine
! must be called.
!
! Temperature in KELVIN.
!
    use variables
    implicit none

    real, intent(in) :: Temp
    SYSTEMP = Temp

    !ITERATION FOR SEGMENT ACITIVITY COEF (PURE SPECIES)
    DO I = 1, COMP
        SEGGAMMAPR (:,I) = 1.0
        DO
            SEGGAMMAOLDPR (:,I) = SEGGAMMAPR (:,I)
            DO mm = 1, 21
                SUMMATION = 0.0
                DO nn = 1, 21
                    SUMMATION = SUMMATION + (SigmaComp(nn,2,I)/QCOSMO(I))*SEGGAMMAOLDPR(nn,I) * &
                    EXP(-DELTAW(I,I,mm,nn)/(RGAS*SYSTEMP))
                END DO
                SEGGAMMAPR(mm,I)=EXP(-LOG(SUMMATION))
                SEGGAMMAPR(mm,I)=(SEGGAMMAPR(mm,I)+SEGGAMMAOLDPR(mm,I))/2.0
                CONPR(mm,I)=ABS((SEGGAMMAPR(mm,I)-SEGGAMMAOLDPR(mm,I))/SEGGAMMAOLDPR(mm,I))
            END DO

        IF (MAXVAL(CONPR) <=0.000001) EXIT
        END DO
    END DO

end subroutine calcSegGammaPure

subroutine calcStavermanGugenheimen
! This routine uses the Staverman Gugenheimen equation to compute the combinatorial term.
!
    USE variables
    IMPLICIT NONE
    REAL :: F, V, VPRIME, BOTF, BOTV

   !THE STAVERMAN-GUGGENHEIM EQUATION
    BOTF = 0.0
    BOTV = 0.0
    F = 0.0
    V = 0.0
    VPRIME = 0.0

    DO I = 1,COMP
        RNORM(I) = RCOSMO(I)/R
        QNORM(I) = QCOSMO(I)/Q
        BOTF = BOTF + FRAC(I)*QNORM(I)
        BOTV = BOTV + FRAC(I)*RNORM(I)
    ENDDO

    DO I = 1, COMP
        F = (QNORM(I))/BOTF
        V = (RNORM(I))/BOTV
        VPRIME = (RNORM(I)**p)/sum(FRAC(:)*(RNORM(:)**p))
        LNGAMMACOMB(I) = LOG(VPRIME) + 1.0 - VPRIME - (5.0)*QNORM(I)*(LOG(V/F) + 1.0 - V/F)
    ENDDO

end subroutine calcStavermanGugenheimen
