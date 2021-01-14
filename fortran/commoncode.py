#!/usr/bin/env python

import numpy as np

SourceBase = """

subroutine ran2(r)
    implicit none
    real(8), intent(out) :: r
    integer, parameter :: NTAB = 32
    integer, parameter :: IM1=2147483563, IM2=2147483399, IMM1=2147483562
    integer, parameter :: IA1=40014, IA2=40692, IQ1=53668, IQ2=52774
    integer, parameter :: IR1=12211, IR2=3791, NDIV=1+IMM1/NTAB  
    real(8), parameter :: AM=1.d0/IM1, EPS=1.2d-7, RNMX=1.d0-EPS
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    integer :: j, k

    if (state(1).le.0) then
        state(1)=max(-state(1),1)
        state(2)=state(1)
        do j=NTAB+8,1,-1
             k=state(1)/IQ1
             state(1)=IA1*(state(1)-k*IQ1)-k*IR1
             if (state(1).lt.0) state(1)=state(1)+IM1
             if (j.le.NTAB) state(3+j)=state(1)
        enddo
        state(3)=state(4)
    endif
    k=state(1)/IQ1
    state(1)=IA1*(state(1)-k*IQ1)-k*IR1
    if (state(1).lt.0) state(1)=state(1)+IM1
    k=state(2)/IQ2
    state(2)=IA2*(state(2)-k*IQ2)-k*IR2
    if (state(2).lt.0) state(2)=state(2)+IM2
    j=1+state(3)/NDIV
    state(3)=state(3+j)-state(2)
    state(3+j)=state(1)
    if(state(3).lt.1)state(3)=state(3)+IMM1
    r=min(AM*state(3),RNMX)
end subroutine

subroutine ran2seed(seedval)
    implicit none
    integer, intent(in) :: seedval
    integer, dimension(3+32) :: state 
    COMMON /ran2data/ state
    state(1) = abs(seedval)
end subroutine

subroutine ran2array(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    integer :: i
    external :: ran2
    real(8) :: r
    do i = 1, n
        call ran2(r)
        rarray(i) = r
    enddo
end subroutine

subroutine ran2int(n, i)
    !returns an integer on [0,n)
    implicit none
    integer, intent(in) :: n
    integer, intent(out) :: i
    external :: ran2
    real(8) :: r
    call ran2(r)
    i = int(real(n) * r)
    i = min(n-1, max(0, i))
end subroutine

subroutine ran2norm(r)
    implicit none
    real(8), intent(out) :: r
    external :: ran2
    real(8) :: r1, r2, rsq
    real(8), save :: rsaved = 0.1
    logical, save :: hassaved = .false.
    if (hassaved) then
        r = rsaved
        hassaved = .false.
    else
        rsq = 2.d0
        do while (rsq == 0.d0 .or. rsq >= 1.d0)
            call ran2(r1)
            call ran2(r2)
            r1 = 2.d0 * r1 - 1.d0
            r2 = 2.d0 * r2 - 1.d0
            rsq = r1*r1 + r2*r2
        enddo
        rsq = sqrt(-2.d0 * log(rsq) / rsq)
        r = r1 * rsq
        rsaved = r2 * rsq
        hassaved = .true.
    endif
end subroutine

subroutine ran2normarray(n, rarray)
    implicit none
    integer, intent(in) :: n
    real(8), dimension(n), intent(out) :: rarray
    external :: ran2norm
    integer :: i
    real(8) :: r
    do i = 1, n
        call ran2norm(r)
        rarray(i) = r
    enddo
end subroutine


subroutine erfc(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the complementary error function erfc(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
end subroutine erfc

subroutine erf(x, value)
    real(8), intent(in) :: x
    real(8), intent(out) :: value
    !Returns the error function erf(x) with fractional 
    !error everywhere less than 1:2^10-7
    real(8) :: t, z
    z = abs(x)
    t = 1./(1.+0.5*z)
    value = t*exp(-z*z - 1.26551223 + t*(1.00002368 + t*(.37409196 + &
          & t*(.09678418 + t*(-.18628806 + t*(.27886807 + t*(-1.13520398 + &
          & t*(1.48851587 + t*(-.82215223 + t*.17087277)))))))))
    if (x < 0.) value = 2. - value
    value = 1. - value
end subroutine erf


integer function GetPairIndFromij(i, j)
    integer, intent(in) :: i, j
    if (i > j) then
        GetPairIndFromij = i * (i + 1) / 2 + j
    else
        GetPairIndFromij = j * (j + 1) / 2 + i
    endif      
end function

subroutine GetijFromPairInd(ind, i, j)
    integer, intent(in) :: ind
    integer, intent(out) :: i, j
    real(8) :: r
    r = dble(ind) + 0.001d0
    r = sqrt(1.d0 + 8.d0 * r) * 0.5d0 - 0.5d0
    i = int(r)
    j = ind - i*(i+1)/2
end subroutine

subroutine RandomRotMat3D(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, sphi, cphi
    real(8) :: Ang, q0, q1, q2, q3
    real(8), dimension(3) :: Vec
    real(8), parameter :: pi = 3.1415926535897931D0
    external :: ran2
    !get a random rotation angle
    call ran2(Ang)
    Ang = (2.*Ang - 1) * MaxAng
    !get a random vector about which to rotate
    call ran2(cphi)
    cphi = 2.*cphi-1.
    cphi = max(-1.,min(1.,cphi))
    sphi = sqrt(1.-cphi*cphi)
    call ran2(theta)
    theta = 2.d0 * pi * theta
    Vec = (/ cos(theta)*sphi, sin(theta)*sphi, cphi /) 
    !make intermediate variables
    q0 = cos(0.5*Ang)
    Vec = sin(0.5*Ang) * Vec
    q1 = Vec(1)
    q2 = Vec(2)
    q3 = Vec(3)
    !assemble the rotation matrix
    RotMat(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
    RotMat(1,2) = 2.*(q1*q2 + q0*q3)
    RotMat(1,3) = 2.*(q1*q3 - q0*q2)
    RotMat(2,1) = 2.*(q1*q2 - q0*q3)
    RotMat(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
    RotMat(2,3) = 2.*(q2*q3 + q0*q1)
    RotMat(3,1) = 2.*(q1*q3 + q0*q2)
    RotMat(3,2) = 2.*(q2*q3 - q0*q1)
    RotMat(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
end subroutine

subroutine RandomRotMatxyz(MaxAng, RotMat)
    real(8), intent(in) :: MaxAng
    real(8), dimension(3,3) :: RotMat
    real(8) :: theta, costheta, sintheta, axis
    external :: ran2
    !get a random rotation angle
    call ran2(theta)
    theta = (2.d0 * theta - 1.d0) * MaxAng
    sintheta = sin(theta)
    costheta = cos(theta)
    !get a random axis about which to rotate
    call ran2(axis)
    RotMat = 0.d0
    if (axis < 0.33333333333d0) then
        RotMat(1,1) = costheta
        RotMat(1,2) = -sintheta
        RotMat(2,1) = sintheta
        RotMat(2,2) = costheta
        RotMat(3,3) = 1.d0
    elseif (axis < 0.66666666667d0) then
        RotMat(2,2) = costheta
        RotMat(2,3) = -sintheta
        RotMat(3,2) = sintheta
        RotMat(3,3) = costheta
        RotMat(1,1) = 1.d0    
    else
        RotMat(1,1) = costheta
        RotMat(1,3) = -sintheta
        RotMat(3,1) = sintheta
        RotMat(3,3) = costheta
        RotMat(2,2) = 1.d0
    endif
end subroutine
"""

def PostLoad(Lib):
    #set the fortran random number generator with a seed based on numpy RNG
    seedval = np.random.get_state()[1][0]
    Lib.Module.ran2seed(seedval)

def AddCommonFortran(Lib):
    """Adds common / base fortran code to Lib."""
    Lib.Code += "\n" + SourceBase + "\n"
    Lib.PostLoadList.append(PostLoad)
