! Program to solve the heat equation using the FTCS finite-difference method
! to compile: gfortran -O3 -ffree-line-length-none -freal-4-real-8 -o HeatEquation heateqn.f90
program HeatEquation
implicit none
integer :: iNX,iNY,i,j,iFinished,iNG,iFileReference,iCounter
real(8) :: rL,rVN,rDX,rInverseDX2,rDY,rInverseDY2,rTime,rTotalTime,rAlpha,rDT,rPi,rL1,rL2,rLInf
real(8),allocatable,dimension(:,:) :: rX,rY,rSolution,rInitialCondition,rSolutionNew,rExact
character(32) :: sFileName,sTimeStep

rPi = 4.d0*atan(1.d0)

! @todo: parallelise with MPI

! Simulation parameters
iNX = 256 ! Number of nodes in x direction
iNY = iNX ! Same number of nodes in y direction
rL = 1.d0 ! Length of domain
rTotalTime = 10.d0 ! End time
rAlpha = 1.d-2 ! Diffusivity
rVN = 0.25d0 ! VN number
iNG = 1 ! Number of ghost cells

rDX = rL/dble(iNX-1)
rInverseDX2 = 1.d0/rDX**2
rDY = rL/dble(iNY-1)
rInverseDY2 = 1.d0/rDY**2
rDT = rVN*min(rDX,rDY)**2/abs(rAlpha)

allocate(rX(1:iNX,1:iNY),rY(1:iNX,1:iNY),rSolution(1-iNG:iNX+iNG,1-iNG:iNY+iNG),rInitialCondition(1:iNX,1:iNY),rSolutionNew(1-iNG:iNX+iNG,1-iNG:iNY+iNG),rExact(1:iNX,1:iNY))

do j=1,iNY 
    do i=1,iNX 
        rX(i,j)=dble(i-1)*rDX
        rY(i,j)=dble(j-1)*rDY
    end do
end do

! Initial condition
rInitialCondition = 0.d0
do j=1,iNY 
    do i=1,iNX 
        rInitialCondition(i,j) = sin(rPi*rX(i,j)/rL)*sin(rPi*rY(i,j)/rL)
    end do
end do
rSolution(1:iNX,1:iNY)=rInitialCondition

! Boundary conditions 
do i=1,iNG
   rSolution(1-i,:) = 0.d0 ! Dirichlet BC
   rSolution(iNX+i,:) = 0.d0 ! Dirichlet BC
end do
do i=j,iNG
    rSolution(:,1-j) = 0.d0 ! Dirichlet BC
    rSolution(:,iNY+j) = 0.d0 ! Dirichlet BC
end do

! Solution computation
rTime = 0.d0
iFinished = 0
rSolutionNew = 0.d0
iCounter = 1
do while(iFinished.eq.0)
   rTime = rTime + rDT
   iCounter = iCounter + 1
   print*,'Timestep: ',iCounter,'Time: ',rTime
   ! Check if finished
   if(rTime.gt.rTotalTime) then
      rDT = rTotalTime-(rTime-rDT)
      rTime = rTotalTime
      iFinished = 1
   elseif(rTime.eq.rTotalTime) then
      iFinished = 1
   end if

   do j=1,iNY 
      do i=1,iNX 
         rSolutionNew(i,j) = rSolution(i,j)+rAlpha*rDT*rInverseDX2*(rSolution(i+1,j)-2.d0*rSolution(i,j)+rSolution(i-1,j)) &
                           & +rAlpha*rDT*rInverseDY2*(rSolution(i,j+1)-2.d0*rSolution(i,j)+rSolution(i,j-1))
      end do
   end do
   rSolution = rSolutionNew
end do

! Error norm calc
rExact = rInitialCondition*exp(-rAlpha*2*rPi**2/rL**2*rTotalTime)
rL1 = 0.d0
rL2 = 0.d0
rLInf = 0.d0
do j=1,iNY 
    do i=1,iNX 
        rL1 = rL1 + abs(rExact(i,j)-rSolution(i,j))
        rL2 = rL2 + (rExact(i,j)-rSolution(i,j))**2
        rLInf = max(abs(rExact(i,j)-rSolution(i,j)),rLInf)
    end do
end do

! rDummyArray = 0.d0
! CALL MPI_ALLREDUCE(rL1,rDummyArray,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IER)
! rL1 = rDummyArray
! rDummyArray = 0.d0
! CALL MPI_ALLREDUCE(rL2,rDummyArray,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IER)
! rL2 = rDummyArray
! rDummyArray = 0.d0
! CALL MPI_ALLREDUCE(rLInf,rDummyArray,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IER)
! rLInf = rDummyArray

! Scale by dx
rL1=rL1*rDX*rDY
rL2=sqrt(rL2*rDX*rDY)
! Write error norm data
sFileName="ErrorNorms.dat"
iFileReference=32
open(iFileReference,file=sFilename,position="append",status="unknown",action="write")
write(iFileReference,FMT=194)rDX,rL1,rL2,rLInf
close(iFileReference)

! Write initial condition data
iFileReference = 12
! write(sTimeStep,FMT='(F8.5)') rTime
! sFileName=trim(adjustl(sTimeStep))//".tec.dat"
sFileName = "0.00000.tec.dat"
print*, 'Writing file: ', sFileName
open (unit=iFileReference,file=sFileName,action='write')
write(iFileReference,*)'VARIABLES= "X","Y","U0"'
write(iFileReference,*) 'ZONE F=POINT',' I=',iNX,' J=',iNY
write(iFileReference,*)'DT=(DOUBLE DOUBLE DOUBLE)'
do j=1,iNY
    do i=1,iNX
    write(iFileReference,FMT=193) rX(i,j),rY(i,j),rInitialCondition(i,j)
    end do
end do
close(iFileReference)

! Write solution data
iFileReference = 13
write(sTimeStep,FMT='(F8.5)') rTime
sFileName=trim(adjustl(sTimeStep))//".tec.dat"
print*, 'Writing file: ', sFileName
open (unit=iFileReference,file=sFileName,action='write')
write(iFileReference,*)'VARIABLES= "X","Y","U","U0"'
write(iFileReference,*) 'ZONE F=POINT',' I=',iNX,' J=',iNY
write(iFileReference,*)'DT=(DOUBLE DOUBLE DOUBLE DOUBLE)'
do j=1,iNY
    do i=1,iNX
    write(iFileReference,FMT=194) rX(i,j),rY(i,j),rSolution(i,j),rExact(i,j)
    end do
end do
close(iFileReference)

193 FORMAT(3E25.15)
194 FORMAT(4E25.15)

end program HeatEquation