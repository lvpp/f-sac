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

subroutine Comps(compoundName, CompoundSigmaProfile, Rtotal, Qtotal)
! This routine receives the compound's name and returns the compound's sigma profile, volume and area.
! This code is very similar to subGroups and Groups subroutines
!
    implicit none

    ! In/out variables
    character(len=25), intent(in) :: compoundName ! Compound's name
    real, dimension(51,4), intent(out) :: CompoundSigmaProfile ! sigma profile
    real, intent(out) :: Rtotal, Qtotal ! volume and area

    ! Local variables
    integer :: n, i, j, ierror ! Auxiliary variables
    real, dimension(:,:), allocatable :: xx, yy, sigmaProfile
    ! xx: 2-collum-matrix that registers: first collum - subgroups; second collum: number of times
    !that subgroup is repeated - note: for now, it only works to compounds that have a maximum of 4 diferent subgroups
    ! yy: 2-collum-matrix that registers the same imformation that matrix xx, but contains only the collums that represent
    ! a subgroups - note: this matrix does not necessarilly contains the same number of lines that xx
    ! sigmaProfile: first collum - charge density; second collum - compound's surface area

    real :: initial, inc ! Auxiliary variables

    integer :: positive_pos, negative_pos, zero_pos ! Finds the positive, zero and negative positions in sigma profile

    ! Each collum of the archive --------------------------
    character(len=25) :: Name, identify, header
    integer :: Fixed,	Counter
    real ::  MW, Rk, Qk, HBacc, HBdonn
    character(len=20) :: CAS, formula
    ! -----------------------------------------------------

    ! Compound's positve, zero and negative charge density and surface area;
    real :: Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, GroupID

    ! Allocates and initiates xx
    allocate(xx(4,2))
    xx = 0

    ! Allocates the sigma profile matrix
    allocate(sigmaProfile(51,4))

    ! Starts the volume and area
    Qtotal = 0
    Rtotal = 0

    ! Sets initial and inc
    initial = - 2.5e-2
    inc = 1e-3

    ! Starts the sigmaProfile matrix
    sigmaProfile = 0

    ! Opens the archive to read
    open(unit=10, file='pars/FSAC1-comps.csv', status='old', action='read')

    ! Starts n and name
    n = 0
    name = ' '

    ! Just to read the archive's header
    read(10,*), header

    ! Counts the lines until the compound we are looking for
    do while (name /= compoundName)
        read(10,*,iostat=ierror) identify, name
        if(ierror /= 0) then
            print*, 'Warning: compound name not found.'
            stop
        endif
        n = n + 1
    enddo
    ! Rewinds the input archive in order to read the lines again
    rewind(10)

    ! Reads the lines just above the compound that we are looking for
    do i = 1, n
        read(10,*) identify
    enddo

    ! Reads the line of the compound, registering each collum
    read(10,*) identify, name, CAS, formula, MW, Fixed, Counter &
    , xx(1,1),  xx(1,2),  xx(2,1),  xx(2,2), xx(3,1),  xx(3,2),  xx(4,1),  xx(4,2)!, xx(5,1)

    ! Counts the number of lines of xx that really exist (not zero)
    n = 0
    do i = 1, 4
        if(xx(i,1) /= 0) then
        n = n + 1
        endif
    enddo

    ! Allocates yy
    allocate(yy(n,2))

    ! Builds the yy matrix
    do j = 1, 2
        do i = 1, n
            yy(i,j) = xx(i,j)
        enddo
    enddo

    zero_pos = 26 ! Zero position - always the element 26

    ! Computes the sigma positions and charges to each line of yy
    do i = 1, n
        call subGroups(yy(i,1), Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, Rk, Qk, HBacc, HBdonn, GroupID)
        positive_pos = nint((sigma_pos - initial)/inc)
        negative_pos = nint((sigma_neg - initial)/inc)

        sigmaProfile(zero_pos,2) = sigmaProfile(zero_pos,2) + yy(i,2)*Q_zero
        sigmaProfile(positive_pos + 1,2) = sigmaProfile(positive_pos + 1,2) + yy(i,2)*Q_pos
        sigmaProfile(negative_pos + 1,2) = sigmaProfile(negative_pos + 1,2) + yy(i,2)*Q_neg

        sigmaProfile(zero_pos,1) = 0.0
        sigmaProfile(positive_pos + 1,1) = sigma_pos
        sigmaProfile(negative_pos + 1,1) = sigma_neg

        sigmaProfile(zero_pos,3) = 0.0
        sigmaProfile(positive_pos + 1,3) = HBacc
        sigmaProfile(negative_pos + 1,3) = HBdonn

        sigmaProfile(zero_pos,4) = GroupID
        sigmaProfile(positive_pos + 1,4) = GroupID
        sigmaProfile(negative_pos + 1,4) = GroupID

        Qtotal = Qtotal + yy(i,2)*Qk
        Rtotal = Rtotal + yy(i,2)*Rk
    enddo

    CompoundSigmaProfile(:,1) = sigmaProfile(:,1)
    CompoundSigmaProfile(:,2) = sigmaProfile(:,2)
    CompoundSigmaProfile(:,3) = sigmaProfile(:,3)
    CompoundSigmaProfile(:,4) = sigmaProfile(:,4)

    deallocate(sigmaProfile, xx, yy)    ! Deallocates all vectors and matrices
    close(10) ! Closes the archive

end subroutine Comps

subroutine subGroups(SubGroupName,Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, Rk, Qk, HBacc, HBdonn, GroupID)
! This routine receives the subgroup's name and returns the subgroups's positve, zero
! and negative charge density and surface area; also returns subgroup's volume and area
! This code is very similar to Comps and Groups subroutines
!
    implicit none

    integer :: n, i ! Auxiliary variables

    ! In/out variables
    real, intent(in) :: SubGroupName ! Subgroup name - enters from subGroup subroutine
    real, intent(out) :: Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, Rk, Qk, HBacc, HBdonn, GroupID ! Subgroup's positve, zero
    ! and negative charge density and surface area; subgroup's volume and area

    ! Each collum of the archive --------------------------
    character(len=15) :: Group, SubGroup
    real :: SubGroupID, MW
    ! -----------------------------------------------------

    ! Opens the archive to read
    open(unit=11, file='pars/FSAC1-subgroups.csv', status='old', action='read')

    ! Initiates n and raeds Group's name
    n = 0
    read(11,*) Group

    ! Counts the lines until the subgroup we are looking for
    do while (SubGroupID /= SubGroupName)
        read(11,*) Group, GroupID, SubGroup, SubGroupID
        n = n + 1
    enddo
    ! Rewinds the input archive in order to read the lines again
    rewind(11)

    ! Reads the lines just above the subgroup that we are looking for
    do i = 1, n
        read(11,*) Group
    enddo

    ! Reads the line of the subgroup, registering each collum
    read(11,*) Group, GroupID, SubGroup, SubGroupID,  MW, Rk, Qk

    ! Computes the group's positve, zero and negative charge density and surface area
    call Groups(GroupID, Qk, Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, HBacc, HBdonn)

    ! Closes the archive
    close(11)

end subroutine subGroups

subroutine Groups(GroupName, Qk, Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, HBacc, HBdonn)
! This routine receives the group's name and returns the groups's positve, zero
! and negative charge density and surface area
! This code is very similar to Comps and subGroups subroutines
!
    implicit none

    integer :: n, i ! Auxiliary variables

    ! In/out variables
    real, intent(in) :: GroupName, Qk ! Group name and total area - enter from Group subroutine
    real, intent(out) :: Q_pos, Q_neg, Q_zero, sigma_pos, sigma_neg, HBacc, HBdonn

    ! Each collum of the archive --------------------------
    character(len=15) :: Group
    real :: GroupID, Charge
    ! -----------------------------------------------------

    ! Opens the archive to read
    open(unit=12, file='pars/FSAC1-groups.csv', status='old', action='read')

    ! Starts n and reads Group's name
    n = 0
    read(12,*) Group


    ! Counts the lines until the group we are looking for
    do while (GroupID /= GroupName)
        read(12,*) Group, GroupID, Charge, Q_pos, Q_neg, sigma_pos, HBacc, HBdonn
        n = n + 1
    enddo
    ! Rewinds the input archive in order to read the lines again
    rewind(12)

    ! Reads the lines just above the subgroup that we are looking for
    do i = 1, n
        read(12,*) Group
    enddo

    ! Reads the line of the group, registering each collum
    read(12,*) Group, GroupID, Charge, Q_pos, Q_neg, sigma_pos, HBacc, HBdonn

    ! Computes the zero surface area
    Q_zero = Qk - (Q_pos + Q_neg)

    ! Computes the negative chrge density - if the positive charge density is zero, so negative is
    if(Q_neg == 0) then
        sigma_neg = 0
        else
        sigma_neg = - sigma_pos*Q_pos/Q_neg
    endif

    ! Closes the archive
    close(12)

end subroutine Groups
