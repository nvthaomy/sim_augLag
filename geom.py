#!/usr/bin/env python

#LAST MODIFIED: 12-09-08

#Conventions:
#all angles are in radians


import numpy as np
import time
import sys

import sim.fortran as fortran

#globals
pi = np.pi
pi2 = 2 * pi
DegPerRad = 180./pi
RadPerDeg = pi/180.


#======== FORTRAN CODE ========
Lib = fortran.Module("geomlib")
Lib.Code = """
function ClipTrig(x)
    implicit none
    real(8) :: ClipTrig
    real(8), intent(in) :: x
    ClipTrig = min(1.D0, max(-1.D0, x))
end function

function NormRad(Rad)
    implicit none
    real(8) :: NormRad
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: pi2 = 2.D0*pi
    real(8), intent(in) :: Rad
    NormRad = mod(Rad + pi, pi2) - pi
    if (NormRad < -pi) NormRad = NormRad + pi2
end function

function NormDeg(Deg)
    implicit none
    real(8) :: NormDeg
    real(8), intent(in) :: Deg
    NormDeg = mod(Deg + 180.D0, 360.D0) - 180.D0
    if (NormDeg < -180.D0) NormDeg = NormDeg + 360.D0
end function

subroutine Centroid(Pos, Ret, N, Dim)
    implicit none
    integer, intent(in) :: N, Dim
    real(8), dimension(N,Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = sum(Pos, 1) / real(N)
end subroutine

subroutine CentroidInd(Pos, AtomInd, Ret, N, Dim, NAtomInd)
    implicit none
    integer, intent(in) :: N, Dim, NAtomInd
    real(8), dimension(N,Dim), intent(in) :: Pos
    integer, dimension(NAtomInd), intent(in) :: AtomInd
    real(8), dimension(Dim), intent(out) :: Ret
    integer :: i
    Ret = 0.
    do i = 1, NAtomInd
        Ret = Ret + Pos(AtomInd(i) + 1, :)
    enddo
    Ret = Ret / dble(NAtomInd)
end subroutine

real(8) function Length(Vec, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Vec
    Length = sqrt(sum(Vec*Vec))
end function

subroutine UnitVec(Vec, Ret, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Vec
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = Vec / sqrt(sum(Vec*Vec))
end subroutine

subroutine CrossProd3(r1, r2, Ret, Dim)
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: r1, r2
    real(8), dimension(Dim), intent(out) :: Ret
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Ret(1) = r1(2)*r2(3)-r1(3)*r2(2)
    Ret(2) = r1(3)*r2(1)-r1(1)*r2(3)
    Ret(3) = r1(1)*r2(2)-r1(2)*r2(1)
end subroutine

real(8) function RadFromTrig(SinVal, CosVal)
    implicit none
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), intent(in) :: SinVal, CosVal
    real(8), external :: ClipTrig
    RadFromTrig = acos(ClipTrig(CosVal))
    if (SinVal < 0.) RadFromTrig = 2.*pi - RadFromTrig
end function

real(8) function DegFromTrig(SinVal, CosVal)
    implicit none
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: DegPerRad = 180.D0/pi
    real(8), intent(in) :: SinVal, CosVal
    real(8), external :: ClipTrig, RadFromTrig
    DegFromTrig = DegPerRad * RadFromTrig(SinVal, CosVal)
end function

real(8) function Angle3Rad(Pos1, Pos2, Pos3, Dim)
    implicit none
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: DegPerRad = 180.D0/pi
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Pos1, Pos2, Pos3
    real(8), dimension(Dim) :: Vec21, Vec23
    real(8) :: Norm, Phi
    real(8), external :: ClipTrig, NormRad
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    if (all(Pos1 == Pos2) .or. all(Pos2 == Pos3)) then
        Angle3Rad = 0.
        return
    endif
    Vec21 = Pos1 - Pos2
    Vec23 = Pos3 - Pos2
    Norm = sqrt(sum(Vec21*Vec21)*sum(Vec23*Vec23))
    Phi = ClipTrig(dot_product(Vec21, Vec23) / Norm)
    Angle3Rad = acos(Phi)
end function


real(8) function Dihedral3Rad(Pos1, Pos2, Pos3, Pos4, Dim)
    implicit none
    real(8), parameter :: pi = 3.1415926535897931D0
    real(8), parameter :: DegPerRad = 180.D0/pi
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Pos1, Pos2, Pos3, Pos4
    real(8), dimension(Dim) :: Vec12, Vec23, Vec34
    real(8), dimension(Dim) :: Norm12, Norm34
    real(8) :: Norm, Phi
    real(8), external :: ClipTrig, NormRad
    external :: CrossProd3
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Vec12 = Pos2 - Pos1
    Vec23 = Pos3 - Pos2
    Vec34 = Pos4 - Pos3
    call CrossProd3(Vec12, Vec23, Norm12, Dim)
    call CrossProd3(Vec23, Vec34, Norm34, Dim)
    Norm = sqrt( sum(Norm12*Norm12) * sum(Norm34*Norm34) )
    Phi = ClipTrig(dot_product(Norm12, Norm34) / Norm)
    Phi = acos(Phi)
    if (dot_product(Vec12, Norm34) < 0.) Phi = -Phi
    Dihedral3Rad = NormRad(Phi)
end function


subroutine GetVecMappingRad(Vec1, Vec2, Vec, Ang, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Vec1, Vec2
    real(8), dimension(Dim), intent(out) :: Vec
    real(8), intent(out) :: Ang
    real(8), dimension(Dim) :: v1, v2
    real(8) :: CosAng
    real(8), external :: ClipTrig
    external :: CrossProd3, UnitVec
    call UnitVec(Vec1, v1, Dim)
    call UnitVec(Vec2, v2, Dim)
    CosAng = ClipTrig(dot_product(v1, v2))
    Ang = acos(CosAng)
    if (CosAng == 1.) then
        Vec = v1
    elseif (CosAng == -1.) then
        Vec(2:Dim) = v1(1:Dim-1)
        Vec(1) = v1(Dim)
        Vec = Vec - dot_product(Vec, v1)
    else
        call CrossProd3(v1, v2, Vec, Dim)
    endif
end subroutine


subroutine RotateArrayAboutPoint(Pos, RotMat, Point, Ret, N, Dim)
    integer, intent(in) :: N, Dim
    real(8), dimension(N, Dim), intent(in) :: Pos
    real(8), dimension(Dim, Dim), intent(in) :: RotMat
    real(8), dimension(Dim), intent(in) :: Point
    real(8), dimension(N, Dim), intent(out) :: Ret
    integer :: i
    do i = 1, N
        Ret(i,:) = Pos(i,:) - Point
    enddo
    Ret = matmul(Ret, RotMat)
    do i = 1, N
        Ret(i,:) = Ret(i,:) + Point
    enddo
end subroutine

subroutine RotatePointAboutPoint(Pos, RotMat, Point, Ret, Dim)
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Pos
    real(8), dimension(Dim, Dim), intent(in) :: RotMat
    real(8), dimension(Dim), intent(in) :: Point
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = matmul(Pos - Point, RotMat) + Point
end subroutine


subroutine RotMat3Rad(Vec, Ang, Ret, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: Vec
    real(8), intent(in) :: Ang
    real(8), dimension(Dim,Dim), intent(out) :: Ret
    real(8), dimension(Dim) :: v
    real(8) :: rcos, rsin
    external :: UnitVec
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Ret = 0.
    rcos = cos(Ang)
    rsin = sin(Ang)
    call UnitVec(Vec, v, 3)
    Ret(1,1) =         rcos + v(1)*v(1)*(1-rcos)
    Ret(2,1) =  v(3) * rsin + v(2)*v(1)*(1-rcos)
    Ret(3,1) = -v(2) * rsin + v(3)*v(1)*(1-rcos)
    Ret(1,2) = -v(3) * rsin + v(1)*v(2)*(1-rcos)
    Ret(2,2) =         rcos + v(2)*v(2)*(1-rcos)
    Ret(3,2) =  v(1) * rsin + v(3)*v(2)*(1-rcos)
    Ret(1,3) =  v(2) * rsin + v(1)*v(3)*(1-rcos)
    Ret(2,3) = -v(1) * rsin + v(2)*v(3)*(1-rcos)
    Ret(3,3) =         rcos + v(3)*v(3)*(1-rcos)
end subroutine


subroutine RotMat3EulerRad(Phi, Theta, Psi, Ret)
    real(8), dimension(3,3), intent(out) :: Ret
    real(8), intent(in) :: Phi, Theta, Psi
    real(8) :: Sin1, Cos1, Sin2, Cos2, Sin3, Cos3
    Sin1 = sin(Phi)
    Cos1 = cos(Phi)
    Sin2 = sin(Theta)
    Cos2 = cos(Theta)
    Sin3 = sin(Psi)
    Cos3 = cos(Psi)
    Ret(1,1) = Cos1*Cos3 - Sin1*Cos2*Sin3
    Ret(1,2) = Sin1*Cos3 + Cos1*Cos2*Sin3
    Ret(1,3) = Sin2*Sin3
    Ret(2,1) = -Cos1*Sin3 - Sin1*Cos2*Cos3
    Ret(2,2) = -Sin1*Sin3 + Cos1*Cos2*Cos3
    Ret(2,3) = Sin2*Cos3
    Ret(3,1) = Sin1*Sin2
    Ret(3,2) = -Cos1*Sin2
    Ret(3,3) = Cos2
end subroutine


subroutine RotMat3Q(q0, q1, q2, q3, Ret)
    implicit none
    real(8), dimension(3,3), intent(out) :: Ret
    real(8), intent(in) :: q0, q1, q2, q3
    Ret(1,1) = q0*q0 + q1*q1 - q2*q2 - q3*q3
    Ret(1,2) = 2.*(q1*q2 + q0*q3)
    Ret(1,3) = 2.*(q1*q3 - q0*q2)
    Ret(2,1) = 2.*(q1*q2 - q0*q3)
    Ret(2,2) = q0*q0 - q1*q1 + q2*q2 - q3*q3
    Ret(2,3) = 2.*(q2*q3 + q0*q1)
    Ret(3,1) = 2.*(q1*q3 + q0*q2)
    Ret(3,2) = 2.*(q2*q3 - q0*q1)
    Ret(3,3) = q0*q0 - q1*q1 - q2*q2 + q3*q3
end subroutine


subroutine EulerFromRotMat3Rad(RotMat, Phi, Theta, Psi, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(Dim,Dim), intent(in) :: RotMat
    real(8), intent(out) :: Phi, Theta, Psi
    real(8) :: sphi
    real(8), external :: RadFromTrig
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Theta = acos(min(1., max(-1., RotMat(3,3))))
    if (RotMat(3,3)==1. .or. (RotMat(3,1)==0. .and. RotMat(3,2)==0.) &
        & .or. (RotMat(1,3)==0. .and. RotMat(2,3)==0.)) then
        Psi = 0.
        Phi = RadFromTrig(RotMat(1,2), RotMat(1,1))
    else
        sphi = sqrt(1.-RotMat(3,3)*RotMat(3,3))
        Phi = RadFromTrig(RotMat(3,1)/sphi, -RotMat(3,2)/sphi)
        Psi = RadFromTrig(RotMat(1,3)/sphi, RotMat(2,3)/sphi)
    endif
end subroutine


subroutine SpherePoints(N, Points)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N,3), intent(out) :: Points
    real(8) :: off, y, phi, r
    real(8) :: inc 
    real(8), parameter :: pi = 3.1415926535897931D0 
    integer :: k
    inc = pi * (3. - sqrt(5.))
    Points = 0.
    off = 2. / real(N)
    do k = 1, N
        y = real(k-1) * off - 1. + (off * 0.5)
        r = sqrt(max(1. - y*y, 0.))
        phi = real(k-1) * inc
        Points(k,1) = cos(phi)*r
        Points(k,2) = y
        Points(k,3) = sin(phi)*r
    enddo
end subroutine


subroutine SphereSurfaceAreas(Pos, Radii, Points, BoxL, Areas, NSphere, NPoints, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8), dimension(NPoints, Dim), intent(in) :: Points
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NSphere), intent(out) :: Areas
    real(8), parameter :: pi = 3.141592653589D0 
    integer, intent(in) :: NSphere, NPoints
    integer :: i, j, k
    real(8), dimension(NPoints,Dim) :: ThisPoints
    real(8) :: AreaPerPoint
    logical, dimension(NPoints) :: Exposed
    real(8), dimension(NSphere) :: RadiiSq
    real(8), dimension(Dim) :: iPos, jPos
    real(8), dimension(Dim) :: iBoxL, Dist
    
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    Areas = 0.
    RadiiSq = Radii*Radii
    do i = 1, NSphere
        iPos = Pos(i,:)
        AreaPerPoint = 4.*pi*Radii(i)**2 / real(NPoints)
        Exposed = .true.
        do k = 1, NPoints
            ThisPoints(k,:) = Points(k,:) * Radii(i) + iPos
        enddo
        do j = 1, NSphere
            if (i == j) cycle
            jPos = Pos(j,:)
            Dist = jPos - iPos
            Dist = Dist - BoxL * anint(Dist * iBoxL)
            jPos = Dist + iPos
            if (.not. any(Exposed)) exit
            !first check if spheres are far from each other
            if (sum(Dist**2) > (Radii(i) + Radii(j))**2) cycle
            do k = 1, NPoints
                if (.not. Exposed(k)) cycle
                if (sum((ThisPoints(k,:) - jPos)**2) < RadiiSq(j)) Exposed(k) = .false.
            enddo
        enddo
        Areas(i) = AreaPerPoint * real(count(Exposed))
    enddo
end subroutine


subroutine SphereVolumes(Pos, Radii, dx, Volumes, NSphere, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8),  intent(in) :: dx
    real(8), dimension(NSphere), intent(out) :: Volumes
    integer, intent(in) :: NSphere
    real(8), dimension(NSphere) :: RadiiSq
    real(8) :: minDistSq, DistSq, dV
    integer :: i,j
    real(8), dimension(Dim) :: Pos2, minPos, maxPos
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    RadiiSq = Radii*Radii
    Volumes = 0.
    dV = dx*dx*dx
    minPos = (/minval(Pos(:,1) - Radii), minval(Pos(:,2) - Radii), minval(Pos(:,3) - Radii)/)
    maxPos = (/maxval(Pos(:,1) + Radii), maxval(Pos(:,2) + Radii), maxval(Pos(:,3) + Radii)/)
    maxPos = maxPos + dx * 0.5    
    !first do a coarse grid check to see which spheres are where
    Pos2 = minPos
    do while (all(Pos2 < maxPos))
        j = 0
        minDistSq = huge(1.d0)
        do i = 1, NSphere
            DistSq = sum((Pos(i,:) - Pos2)**2)
            if (DistSq < minDistSq .and. DistSq < RadiiSq(i)) then
                minDistSq = DistSq
                j = i
            endif
        enddo
        if (j > 0) Volumes(j) = Volumes(j) + dV
        Pos2(1) = Pos2(1) + dx
        do i = 1, 2
            if (Pos2(i) >= maxPos(i)) then
                Pos2(i) = minPos(i)
                Pos2(i+1) = Pos2(i+1) + dx
            endif
        enddo
    enddo   
end subroutine


!======== POSITIONS ========


subroutine cubiclattice(N, Dim, Pos)
    integer, intent(in) :: N, Dim
    real(8), dimension(N, Dim), intent(out) :: Pos
    integer :: i, j
    integer, dimension(Dim) :: Counts, MaxCounts
    real(8), dimension(Dim) :: Spc
    real(8) :: Distrib
    OffSet = 0
    Distrib = dble(N)**(1.d0/dble(Dim))
    if (abs(Distrib-int(Distrib))>1.d-12) then
        Distrib = dble(int(Distrib)) + 1.d0
    endif
    MaxCounts = dnint(Distrib)
    Spc = 1.d0 / Distrib
    Counts = 0
    do i = 1, N
        Pos(i,:) = dble(Counts)*Spc - 0.5d0
        Counts(1) = Counts(1) + 1
        do j = 1, (Dim-1), 1
            if (Counts(j) == MaxCounts(j)) then
                Counts(j) = 0
                Counts(j+1) = Counts(j+1) + 1
            endif
        enddo
    enddo
end subroutine 


subroutine fcclattice(N, Dim, Pos)
    integer, intent(in) :: N, Dim
    real(8), dimension(N, Dim), intent(out) :: Pos
    integer :: i, j
    integer :: OffSet
    integer, dimension(Dim) :: Counts, MaxCounts
    real(8), dimension(Dim) :: Spc
    real(8) :: Distrib
    OffSet = 0
    Distrib = (2.d0*dble(N))**(1.d0/dble(Dim))
    if (abs(Distrib-int(Distrib))>1.d-12) then
        Distrib = dble(int(Distrib)) + 1.d0
    endif
    MaxCounts = dnint(Distrib)
    MaxCounts(1) = dnint(Distrib/2.d0)
    Spc = 1.d0 / Distrib
    Counts = 0
    do i = 1, N
        Pos(i,:) = dble(Counts)*Spc - 0.5d0
        Pos(i,1) = dble(2*Counts(1))*Spc(1) - 0.5d0 + dble(OffSet)*Spc(1)
        Counts(1) = Counts(1) + 1
        do j = 1, (Dim-1), 1
            if (Counts(j) == MaxCounts(j)) then
                Counts(j) = 0
                Counts(j+1) = Counts(j+1) + 1
                Offset = 1 - OffSet
            endif
        enddo
    enddo
end subroutine 


subroutine bcclattice(N, Dim, Pos)
    !only works for 3D
    integer, intent(in) :: N, Dim
    real(8), dimension(N, Dim), intent(out) :: Pos
    integer :: i, j
    integer, dimension(Dim) :: OffSet
    integer, dimension(Dim) :: Counts, MaxCounts
    real(8), dimension(Dim) :: Spc
    real(8) :: Distrib
    OffSet = 0
    Distrib = (4.d0 * dble(N))**(1.d0/dble(Dim))
    if (abs(Distrib-int(Distrib))>1.d-12) then
        Distrib = dble(int(Distrib)) + 1.d0
    endif
    MaxCounts = anint(Distrib) 
    MaxCounts(1:2) = anint(Distrib/2.d0)
    Spc = 1.d0 / Distrib
    Counts = 0
    do i = 1, N
        Pos(i,:) = dble(2.d0*Counts)*Spc - 0.5d0 + dble(OffSet)*Spc
        Pos(i,Dim) = dble(Counts(Dim))*Spc(Dim) - 0.5d0
        Counts(1) = Counts(1) + 1
        do j = 1, (Dim-1), 1
            if (Counts(j) == MaxCounts(j)) then
                Counts(j) = 0
                Counts(j+1) = Counts(j+1) + 1
                if (j==Dim-1) Offset(1:Dim-1) = 1. - Offset(1:Dim-1)
            endif
        enddo
    enddo
end subroutine


subroutine minimizepos(Pos, NSteps, ds, N, Dim)
    !only works for 3D
    integer, intent(in) :: N, Dim, NSteps
    real(8), intent(in) :: ds
    real(8), dimension(N, Dim), intent(inout) :: Pos
!f2py intent(in,out,inplace) :: Pos
    real(8), dimension(N, Dim) :: Force
    real(8), dimension(Dim) :: Fij, rij, Posmin, Posmax, Poscen
    real(8) :: dij2
    integer :: i, j, k
    external :: Centroid
    !normalize
    call Centroid(Pos, Poscen, N, Dim)
    do i = 1, N
        Pos(i,:) = Pos(i,:) - PosCen
    enddo
    Posmin = minval(Pos, 1)
    Posmax = maxval(Pos, 1)
    do i = 1, N
        Pos(i,:) = (Pos(i,:) - PosMin) / (PosMax - PosMin) - 0.5d0
    enddo
    !minimize
    do k = 1, NSteps
        !compute forces
        Forces = 0.d0
        do i = 1, N
            !pair forces, as -1/r
            do j = i+1, N
                rij = Pos(j,:) - Pos(i,:)
                rij = rij - dnint(rij)
                dij2 = dot_product(rij, rij)
                Fij = -rij / (dij2 + 0.001d0)
                Force(i,:) = Force(i,:) + Fij
                Force(j,:) = Force(j,:) - Fij
            enddo
        enddo
        !update
        Pos = Pos + ds * Force
    enddo
    !check bounds
    Pos = Pos - dnint(Pos)
end subroutine



!======== RMSD ========
function drms(Pos1, Pos2, NPos, Dim)
    implicit none
    real(8) :: dRMS
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos1, Pos2
    integer :: i,j
    real(8) :: d1, d2
    dRMS = 0.
    do i = 1, NPos
        do j = i+1, NPos
            d1 = sqrt(sum((Pos1(i,:)-Pos1(j,:))**2))
            d2 = sqrt(sum((Pos2(i,:)-Pos2(j,:))**2))
            dRMS = dRMS + (d1-d2)**2
         enddo
    enddo
    dRMS = sqrt(dRMS * 2.d0 / (NPos * (NPos - 1)))
end function




!======== CONTACT MAPS AND CLUSTERING ========

subroutine contactmap(Pos, BoxL, Cutoff, CMap, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), intent(in) :: Cutoff
    integer, dimension(NPos, NPos), intent(out) :: CMap
    integer :: i, j
    real(8), dimension(Dim) :: distvec, iBoxL
    real(8) :: distsq, cutsq 
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    cutsq = Cutoff * Cutoff
    CMap = 0
    do i = 1, NPos
        do j = i+1, NPos
            distvec = Pos(i,:) - Pos(j,:)
            distvec = distvec - BoxL * anint(distvec * iBoxL)
            distsq = sum(distvec * distvec)
            if (distsq <= cutsq) then
                CMap(i,j) = 1.d0
                CMap(j,i) = 1.d0
            endif
        enddo
    enddo            
end subroutine



!======== MINIMUM IMAGING ========

subroutine minimage(Pos, BoxL, ReimagedPos, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: ReimagedPos
    real(8), dimension(Dim) :: iBoxL
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    do i = 1, NPos
        ReimagedPos(i,:) = Pos(i,:) - BoxL * anint(Pos(i,:) * iBoxL)
    enddo            
end subroutine

subroutine minimagen(Pos, BoxL, n, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: n
    real(8), dimension(Dim) :: iBoxL
    integer :: i
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    do i = 1, NPos
        n(i,:) = anint(Pos(i,:) * iBoxL)
    enddo            
end subroutine

subroutine reimage(Pos, RefPos, BoxL, ReimagedPos, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: RefPos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: ReimagedPos
    integer :: i
    real(8), dimension(Dim) :: distvec, iBoxL
    real(8) :: distsq, cutsq 
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    do i = 1, NPos
        distvec = Pos(i,:) - RefPos
        distvec = distvec - BoxL * anint(distvec * iBoxL)
        ReimagedPos(i,:) = RefPos + distvec
    enddo            
end subroutine

"""


#======== VECTOR FUNCTIONS ========

def Length(Vec):
  "Returns the length of a vector."
  return Lib.Module.length(Vec)
  
def UnitVec(Vec):
  "Returns a vector of unit length."
  return Lib.Module.unitvec(Vec)

def RandVec(Dim = 3):
  "Returns a random unit vector on a sphere."
  while True:
    Vec = 2.*np.random.rand(Dim) - 1.
    if (Vec*Vec).sum() <= 1: break
  return UnitVec(Vec)


#======== ANGLE FUNCTIONS ========

def NormRad(Rad):
  "Normalizes Rad to be within (-pi,pi)."
  return Lib.Module.normrad(Rad)

def NormDeg(Deg):
  "Normalizes Deg to be within (-180.,180.)."
  return Lib.Module.normdeg(Deg)

def NearestAngle(Rad, RefRad):
  "Returns the nearest Rad +- n*2pi to RefRad."
  return ((Rad - RefRad + pi) % pi2) + RefRad - pi

def RadFromTrig(SinVal, CosVal):
  "Determines the angle in radians from sin and cosine values."
  return Lib.Module.radfromtrig(SinVal, CosVal)

def DegFromTrig(SinVal, CosVal):
  "Determines the angle in degrees from sin and cosine values."
  return Lib.Module.degfromtrig(SinVal, CosVal)

def Angle(Pos1, Pos2, Pos3):
  "Calculates the angle formed by three positions."
  return Lib.Module.angle3rad(Pos1, Pos2, Pos3)

def Dihedral(Pos1, Pos2, Pos3, Pos4):
  "Calculates the dihedral angle formed by four positions."
  return Lib.Module.dihedral3rad(Pos1, Pos2, Pos3, Pos4)

def GetVecMapping(Vec1, Vec2):
  """Returns the axis and angle between two vectors.""
Returns Vec, Ang such that Vec1 = dot(Vec2, RotMat(Vec, Ang))"""
  return Lib.Module.getvecmappingrad(Vec1, Vec2)


#======== SPHERE ROUTINES =======

def SpherePoints(NPoints):
  """Returns NPoints points equally spaced on a unit sphere, using the
golden spiral algorithm.  Return array is dimensions (3,NPoints)"""
  return Lib.Module.spherepoints(NPoints)

def SphereSurfaceAreas(Pos, Radii, NPoints = 1000, BoxL = [-1, -1, -1]):
  """Uses the Shrake-Rupley algorithm to compute the surface area
of a set of spheres with center positions in Pos and radii in Radii.
NPoints is the number of points to use on each sphere.
BoxL is optional for minimum imaging.
Returns an array of the surface area per sphere."""
  #make a unit sphere
  Points = SpherePoints(NPoints)
  return Lib.Module.spheresurfaceareas(Pos, Radii, Points, BoxL)

def SphereVolumes(Pos, Radii, dx = 1.0):
  """Uses a grid-based algorithm to compute the volume of a set of
overlapping spheres with center positions in Pos and radii in Radii.
dx is grid bin width for dividing up space along each axis.
FineScale is used to do a faster calculation by first coarse-
graining space into boxes of dimension dx and then dx/FineScale.
FineScale only works with the library.
Returns an array of the volume per sphere."""
  return Lib.Module.spherevolumes(Pos, Radii, dx)


#======== ROTATION ROUTINES ========

def RotMat(Vec, Ang):
  "Returns the rotation matrix about a vector Ang degrees."
  return Lib.Module.rotmat3rad(Vec, Ang)

def RotMatEuler(Phi, Theta, Psi):
  "Returns the rotation matrix for Euler angles Phi, Theta, and Psi."
  return Lib.Module.rotmat3eulerrad(Phi, Theta, Psi)

def RotMatQ(q0, q1, q2, q3):
  "Returns the rotation matrix for quarternions q0, q1, q2, and q3."
  return Lib.Module.rotmat3q(q0, q1, q2, q3)

def RotateAboutPoint(Pos, RotMatrix, Point):
    "Rotates a position matrix through a point."
    #translate, rotate the positions, and translate back
    if len(Pos.shape) == 1:
      return Lib.Module.rotatepointaboutpoint(Pos, RotMatrix, Point)
    else:
      return Lib.Module.rotatearrayaboutpoint(Pos, RotMatrix, Point)

def EulerFromRotMat(RotMatrix):
  "Gets the Euler angles implied by a rotation matrix."
  return Lib.Module.eulerfromrotmat3rad(RotMatrix)

def RandRotMat(dAng):
  "Returns a random rotation matrix."
  Vec = RandVec(3)
  Ang =  dAng * (2*np.random.random() - 1)
  return RotMat(Vec, Ang)	

def RandRotMatEuler(dAng, Phi, Theta, Psi):
  "Returns a random rotation matrix and new Euler angles."
  r = np.dot(RotMatEuler(Phi, Theta, Psi), RandRotMat(dAng))
  Phi, Theta, Psi = EulerFromRotMat(r)
  return r, Phi, Theta, Psi

def RotMatMapped(Pos1, Pos2):
  """Extracts the rotation matrix from two arrays of three points each.
Returns RotMat such that Pos1 is aligned to dot(Pos2, RotMat)"""
  #check for identity
  if np.all(Pos1 == Pos2):
    return np.identity(Pos1.shape[1], Pos1.dtype)
  #align first vector
  Vec, Ang = GetVecMapping(Pos1[1] - Pos1[0], Pos2[1] - Pos2[0])
  rm = RotMat(Vec, Ang)
  #align pos2 to pos1
  Pos2 = np.dot(Pos2, rm)
  #from the remaining vec, get the perpindicular components to the aligned
  vperp = UnitVec(Pos1[1] - Pos1[0])
  v1 = Pos1[2] - Pos1[1]
  v2 = Pos2[2] - Pos2[1]
  v1 = v1 - np.dot(v1, vperp) * vperp
  v2 = v2 - np.dot(v2, vperp) * vperp
  #check to make sure points are not colinear
  if sum(v1*v1) == 0 or sum(v2*v2) == 0:
    return rm
  else:
    Vec, Ang = GetVecMapping(v1, v2)
    rm = np.dot(rm, RotMat(Vec, Ang))
    return rm

def RotMatRMSD(Pos1, Pos2, Center = True):
  """Extracts the rotation matrix from two arrays of position coordinates,
that would minimize the RMSD when both arrays are aligned.
Returns RotMat such that Pos1 is aligned to dot(Pos2, RotMat)"""
  #check for identity
  if np.all(Pos1 == Pos2):
    return np.identity(Pos1.shape[1], Pos1.dtype)
  #make copies
  p1 = Pos1.copy() 
  p2 = Pos2.copy()
  if Center:
    p1 = p1 - p1.mean(axis=0)
    p2 = p2 - p2.mean(axis=0)
  d1, d2 = np.shape(p1)
  #calculate correlation matrix
  C = np.dot(p2.transpose(), p1)
  #get singular value decomp
  V, S, Wt = np.linalg.svd(C)
  #check determinants
  Vdet = np.linalg.det(V)
  Wtdet = np.linalg.det(Wt)
  #correct for possible reflection
  Wt[-1,:] = Wt[-1,:] * (Wtdet * Vdet)
  #return rotation matrix
  return np.dot(V, Wt)


def AlignmentRMSD(Pos1, Pos2, Center = True):
  """Returns the translation vectors, rotation matrix, and sum of
residuals squared (Pos1Vec, Pos2Vec, RotMat, Resid), for aligning
Pos1 to Pos2, such that Pos1 + Pos1Vec is aligned to
dot(Pos2 + Pos2Vec, RotMat)."""
  d1, d2 = Pos1.shape
  #get centers
  if Center:
    Pos1Vec = -np.average(Pos1, axis=0)
    Pos2Vec = -np.average(Pos2, axis=0)
    p1 = Pos1 + Pos1Vec
    p2 = Pos2 + Pos2Vec
  else:
    Pos1Vec, Pos2Vec = np.zeros(d2, Pos1.dtype), np.zeros(d2, Pos2.dtype)
    p1, p2 = Pos1, Pos2
  #check for identity
  if np.all(p1 == p2):
    return Pos1Vec, Pos2Vec, np.identity(d2, Pos1.dtype), 0.
  #calculate E0
  E0 = np.sum(p1*p1, axis=None) + np.sum(p2*p2, axis=None)
  #calculate correlation matrix
  C = np.dot(np.transpose(p2), p1)
  #get singular value decomp
  V, S, Wt = np.linalg.svd(C)
  #check determinants
  Vdet = np.linalg.det(V)
  Wtdet = np.linalg.det(Wt)
  #if it's a reflection, reflect along lowest eigenvalue
  if Vdet*Wtdet < 0.: S[-1] = -S[-1]
  #compute risiduals
  Residuals = np.max(E0 - 2. * np.sum(S), 0.)
  #correct for possible reflection
  Wt[-1,:] = Wt[-1,:] * (Wtdet * Vdet)
  #calculate rotation matrix
  U = np.dot(V, Wt)
  return Pos1Vec, Pos2Vec, U, Residuals


def RMSD(Pos1, Pos2, Center = True):
  """Returns the mean squared displacement for aligning
Pos1 to Pos2, such that Pos1 + Pos1Vec is aligned to
dot(Pos2 + Pos2Vec, RotMat)."""
  d1, d2 = Pos1.shape
  #get centers
  if Center:
    p1 = Pos1 - np.average(Pos1, axis=0)
    p2 = Pos2 - np.average(Pos2, axis=0)
  else:
    p1, p2 = Pos1, Pos2
  #check for identity
  if np.all(p1 == p2) or d1 == 0: return 0.
  #calculate E0
  E0 = np.sum(p1*p1, axis=None) + np.sum(p2*p2, axis=None)
  #calculate correlation matrix
  C = np.dot(np.transpose(p2), p1)
  #get singular value decomp
  V, S, Wt = np.linalg.svd(C)
  #check determinants
  Vdet = np.linalg.det(V)
  Wtdet = np.linalg.det(Wt)
  #if it's a reflection, reflect along lowest eigenvalue
  if Vdet*Wtdet < 0.: S[-1] = -S[-1]
  #compute residuals
  Residuals = np.max(E0 - 2. * np.sum(S), 0.)
  return np.sqrt(Residuals / d1)

def dRMSD(Pos1, Pos2):
  """Returns the distance-based RMSD."""
  return Lib.Module.drms(Pos1, Pos2)


#======== LATTICES ========
def CubicLattice(N, Dim = 3):
    "Returns a cubic lattice on [-0.5, 0.5) in each dimension."
    return Lib.Module.cubiclattice(N, Dim)

def BCCLattice(N, Dim = 3):
    "Returns an BCC lattice on [-0.5, 0.5) in each dimension."
    return Lib.Module.bcclattice(N, Dim)

def FCCLattice(N, Dim = 3):
    "Returns an FCC lattice on [-0.5, 0.5) in each dimension."
    return Lib.Module.fcclattice(N, Dim)

def RandomMinPos(N, Dim = 3, NSteps = 100, ds = 0.001):
    """Returns a random array of positions on [-0.5, 0.5);
minimization is performed to place positions apart."""
    Pos = np.random.rand(N, Dim)
    return Lib.Module.minimizepos(Pos, NSteps, ds)

def MinimizePos(Pos, NSteps = 100, ds = 0.001):
    """Removes overlaps of positions in Pos so as to minimize overlaps;
returned array has positions in [-0.5, 0.5)."""
    return Lib.Module.minimizepos(Pos, NSteps, ds)
    

#======== CONTACT MAPS AND CLUSTERING ========
def ContactMap(Pos, BoxL, Cutoff):
    """Returns a contact map for sites in Pos.  BoxL will be used to minimum
image pair distances, and should be an array of same dimensions as Pos[i,:].
If the system is not periodic in one or more dimensions, use BoxL[i] = 0"""
    return Lib.Module.contactmap(Pos, BoxL, Cutoff)
    
def ClusterStats(Pos, BoxL, Cutoff):
    """Performs clustering analysis of positions in Pos. Will group positions
in contact into a single group.  Returns (ClustDist, ClustGroups) where 
ClustDist gives the number of clusters of a given size (i.g., ClustDist[2]
tells how many clusters have just two positions) and ClustGroups is a list of 
lists each which contains position indices in clusters. BoxL will be used to 
minimum image pair distances, and should be an array of same dimensions as Pos[i,:].
If the system is not periodic in one or more dimensions, use BoxL[i] = 0"""
    cm = Lib.Module.contactmap(Pos, BoxL, Cutoff)
    N = len(Pos)
    ClustGroups = []
    for i in range(N):
        Neighs = [j for (j,cmj) in enumerate(cm[i,:])
                  if cmj > 0 and i != j]
        ClustGroups.append(set([i] + Neighs))
    Cont = True
    while Cont:
        Cont = False
        for (i, Group1) in enumerate(ClustGroups):
            for (j, Group2) in enumerate(ClustGroups[i+1:]):
                if len(Group1.intersection(Group2)):
                    ClustGroups[i] = Group1.union(Group2)
                    del ClustGroups[j + i + 1]
                    Cont = True
                    break
            if Cont:
                break
    #sort groups
    ClustGroups = [sorted(list(Group)) for Group in ClustGroups]
    ClustGroups.sort()
    #make distribution
    ClustDist = np.zeros(N+1, dtype=int)
    for Group in ClustGroups:
        ClustDist[len(Group)] += 1
    return ClustDist, ClustGroups
                    
            
#======== REIMAGING ========
def Minimage(Pos, BoxL):
    """Returns a reimaged version of Pos to the central box.
If the system is not periodic in one or more dimensions, use BoxL[i] = 0"""
    if len(Pos.shape) == 2:
        return Lib.Module.minimage(Pos, BoxL)
    else:
        return Lib.Module.minimage(Pos[np.newaxis,:], BoxL)
    
def MinimageN(Pos, BoxL):
    """Returns a the number of BoxL tranlations to reimage each position in Pos.
If the system is not periodic in one or more dimensions, use BoxL[i] = 0"""
    if len(Pos.shape) == 2:
        return Lib.Module.minimagen(Pos, BoxL)
    else:
        return Lib.Module.minimagen(Pos[np.newaxis,:], BoxL)
            
def Reimage(Pos, RefPos, BoxL):
    """Returns a reimaged version of Pos with images taken neares to RefPos.
If the system is not periodic in one or more dimensions, use BoxL[i] = 0"""
    if len(Pos.shape) == 2:
        return Lib.Module.reimage(Pos, RefPos, BoxL)
    else:
        return Lib.Module.reimage(Pos[np.newaxis,:], RefPos, BoxL)
        

#======== COMPILING ========
Lib.KeepSource = True
Lib.Load()


#======== TEST ========
def Test():
    "Runs tests of the compiled routines."
    NLoop = 1000
    v1 = np.array([1,-3,5], float)
    v2 = np.array([2,-6,-2], float)
    v3 = np.array([-1,-1,5], float)
    v4 = np.array([3,2,1], float)
    M = RotMat(v1, 30.)
    Pos = np.random.rand(15).reshape((5,3))
    Radii = np.random.rand(5)
    for j in range(NLoop):
      l = []
      l.extend([Length(v1), UnitVec(v1)])
      l.extend([NormRad(1.2*pi), NormRad(-1.2*pi), NormDeg(-185.), NormDeg(185.)])
      l.extend([RadFromTrig(1,0), RadFromTrig(0,-1), RadFromTrig(np.sin(1.), np.cos(1.))])
      l.extend([DegFromTrig(1,0), DegFromTrig(0,-1), DegFromTrig(np.sin(1.), np.cos(1.))])
      l.extend([Dihedral(v1, v2, v3, v4), Dihedral(v2, v4, v1, v3)])
      l.extend([GetVecMapping(v1, v2)])
      l.extend([np.abs(RotMat(v1, 51.2)).sum(), np.abs(RotMatEuler(32., 24., 80.)).sum(), \
                np.abs(RotMatQ(1., .2, .4, .894)).sum()])
      l.extend([np.abs(RotateAboutPoint(Pos, M, v2)).sum()])
      l.extend([EulerFromRotMat(M)])
      l.extend([SpherePoints(5)])
    l.extend([SphereSurfaceAreas(Pos, Radii, 1000)])
    l.extend([SphereVolumes(Pos, Radii, 1.)])
    l.extend([SphereVolumes(Pos, Radii, 0.1)])
    for x in l: print x
    print CubicLattice(9, 2)
    print BCCLattice(9, 2)
    print FCCLattice(9, 2)
    
def TestCluster():
    NPos = 1000
    Dim = 3
    BoxL = np.array([10.]*Dim, dtype=float)
    Cutoff = 0.7
    Pos = np.random.rand(NPos * Dim).reshape((NPos, Dim)) * BoxL - 0.5 * BoxL
    ClustDist, ClustGroups = ClusterStats(Pos, BoxL, Cutoff)
    print "Cluster distribution:"
    for (n, ct) in enumerate(ClustDist):
        print "   %7imers : %7i" % (n, ct)
    print "Cluster groups:"
    for Group in ClustGroups:
        print "  ", Group




