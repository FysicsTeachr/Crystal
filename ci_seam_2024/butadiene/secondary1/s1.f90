program main
implicit none
integer*1, dimension(24):: qq,qq2, mcPoint
real, dimension(24):: qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real*8,dimension(:),allocatable:: pot, pot2
real::start,finish,time,dumdu,now,addIt,qShrink, acceptance
real*8:: Ecut0, vbar, ein, ein2, vbarTMP,se

integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y,xx,same,energyChecked, rand, oldRand
integer:: gradLine, termCount, keepStoring, rand2, oldRand2, yyy, tsign, same2, initi, signed

character (len=11) :: numnum(3)
character (len=10) ::output1
character (len=10) ::set1, set2, set1b, set2b
real, Dimension(123)::grad, coup
logical::there
character (len=1024) :: text
character (len=1024) :: text2, text3, text4, text5, text6
character (len=1024) :: xxx,xxy
character(len=1024):: reader

real:: CX, CY, CZ, NX,NY,NZ, H1X,H1Y,H1Z, H2X, H2Y,H2Z, H3X,H3Y,H3Z, nscal, seamCut
real::CXp,CYp,CZp,NXp,NYp,NZp,H1Xp,H1Yp,H1Zp,H2Xp,H2Yp,H2Zp,H3Xp,H3Yp,H3Zp,oldH1xp,oldH1yp
real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup

real::C1X,C1Y,C1Z,C2X,C2Y,C2Z,C3X,C3Y,C3Z,C4X,C4Y,C4Z
real::H7X,H7Y,H7Z,H8X,H8Y,H8Z,H9X,H9Y,H9Z,H10X,H10Y,H10Z,H5X,H5Y,H5Z,H6X,H6Y,H6Z

real, dimension(8,123):: G
real, dimension(8):: gradAlong, coupAlong, qStretch
real, dimension(:,:), allocatable:: gradMatrix, coupMatrix

integer:: current, shellCounter, waiter,unitt,ax, dun
character (len=13) ::measure, foundp, shell,foundling
character::from

waiter=2
oldRand=0
shell="../shelll.txt"
measure="../measur.txt"
foundp="../foundp.txt"
foundling="foundling.txt"

  do
    open(10, file=measure)
    read(10,*,iostat=irr) current, shellCounter,from, dun
    if(irr/=0.or.current==0)then
      close(10)
      cycle
    end if
    rewind(10)
    if(current==shellCounter) then
      write(10,*) 0, shellCounter, "s",dun
    else
      write(10,*) current+1, shellCounter, "s",dun
    end if
    close(10)
    exit
  end do

!  do
!   if(current ==0) then
!     call sleepqq(100)
!
!     do
!       open(10, file=measure)
!       read(10,*,iostat=irr) current, shellCounter,from, dun
!       if(irr/=0)then
!         close(10)
!         cycle
!       end if
!       close(10)
!       exit
!     end do
!
!
!     cycle
!   else
!     exit
!   end if
!  end do  

!  if(current==shellCounter) then
!    open(12, file=measure)
!    write(12,*) 0, shellCounter, "s",dun+1
!    close(12)
!  else
!    open(12, file=measure)
!    write(12,*) current+1, shellCounter, "s",dun+1
!    close(12)
!  end if

  open(132,file="coord.txt")
  write(132,*) current, shellCounter

  open(10, file= shell)
!  read(10,*)
  do i=1,current-1
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
    read(10,*)
  end do

  read(10,*,iostat=irr) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
  write(132,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
  read(10,*,iostat=irr)  qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
  write(132,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
  read(10,*,iostat=irr) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
  write(132,*) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
  read(10,*,iostat=irr)  qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
  write(132,*) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
  read(10,*,iostat=irr) qq2(21), qq2(22), qq2(23), qq2(24)
  write(132,*) qq2(21), qq2(22), qq2(23), qq2(24)
  close(10)

write(*,*) qq2

  qt = 1     
  do i=7,24
    qt(i)=1  
  end do
    qt(1) = qt(1)/2.0
    qt(2) = qt(2)/2.0
  qt=qq2*1.0d0*qt

!write(*,*) qt

! C1X=1.50313 + qt(1) + qt(3)
! C1Y=  0.921457 + qt(5)
! C1Z=    -0.637969 + qt(2)
! C2X=  0.757837 + qt(1)
! C2Y=  0
! C2Z=  0
! C3X=  -0.757837 - qt(1)
! C3Y=  0
! C3Z=    0
! C4X=  -1.36276 - qt(1) - qt(4)
! C4Y=  -1.16938 - qt(6)
! C4Z=   -0.655378 + qt(2)
! H5X=  2.57549 + qt(1) + qt(3) + qt(7)
! H5Y=  0.900799 + qt(5) + qt(8)
! H5Z=   -0.611677 + qt(2) + qt(9)
! H6X=  1.05318 + qt(1) + qt(10)
! H6Y=  1.71906 + qt(11)
! H6Z=    -1.20076 + qt(12)
! H7X=  1.23105 + qt(1) + qt(13)
! H7Y=  -0.79691 + qt(14)
! H7Z=  0.544652 + qt(15)
! H8X=  -1.12147 - qt(1) + qt(16)
! H8Y=  0.920091 + qt(17)
! H8Z=  -0.450208 + qt(18)
! H9X=  -1.86258 - qt(1) + qt(19) - qt(4)
! H9Y=  -2.05995 + qt(20) -    qt(6)
! H9Z=  -0.915735 + qt(2) + qt(21)
! H10X= -1.10957 - qt(1) + qt(22)
! H10Y=  0.0229366 + qt(23)
! H10Z=  1.02401 + qt(24)

 C1X= 1.48227 + qt(1) + qt(3)
 C1Y=  1.11308 + qt(5)
 C1Z=   -0.00486039 + qt(2)
 C2X=  0.681794  + qt(1)
 C2Y=  0
 C2Z=  0
 C3X=  -0.681794 - qt(1)
 C3Y=  0
 C3Z=    0
 C4X=  -1.76433 - qt(1) - qt(4)
 C4Y=  -1.00594- qt(6)
 C4Z=   -0.00486039 + qt(2)
 H5X=  2.5508+ qt(1) + qt(3) + qt(7)
 H5Y=  1.03324+ qt(5) + qt(8)
 H5Z=   -0.00963474 + qt(2) + qt(9)
 H6X=  1.05751 + qt(1) + qt(10)
 H6Y=  2.09595 + qt(11)
 H6Z=    -0.014923  + qt(12)
 H7X=  1.16363  + qt(1) + qt(13)
 H7Y=  -0.970061 + qt(14)
 H7Z=   0.00259615 + qt(15)
 H8X=  -2.52962 - qt(1) + qt(16)
 H8Y=  -0.757754 + qt(17)
 H8Z=  -0.7326 + qt(18)
 H9X=  -1.38611 - qt(1) + qt(19) - qt(4)
 H9Y=  -1.99998 + qt(20) -    qt(6)
 H9Z=  -0.245768 + qt(2) + qt(21)
 H10X= -2.25317 - qt(1) + qt(22)
 H10Y=  -1.06462 + qt(23)
 H10Z=  0.962647 + qt(24)

if(                                                 &
                ((c1x-h5x)**2+(c1y-h5y)**2+(c1z-h5z)**2 > 2.11**2 .and. &
                 (c2x-h5x)**2+(c2y-h5y)**2+(c2z-h5z)**2 > 2.11**2 .and. &
                 (c3x-h5x)**2+(c3y-h5y)**2+(c3z-h5z)**2 > 2.11**2 .and. &
                 (c4x-h5x)**2+(c4y-h5y)**2+(c4z-h5z)**2 > 2.11**2) &
                 .or.                                             &
                ((c1x-h6x)**2+(c1y-h6y)**2+(c1z-h6z)**2 > 2.11**2 .and. &
                 (c2x-h6x)**2+(c2y-h6y)**2+(c2z-h6z)**2 > 2.11**2 .and. &
                 (c3x-h6x)**2+(c3y-h6y)**2+(c3z-h6z)**2 > 2.11**2 .and. &
                 (c4x-h6x)**2+(c4y-h6y)**2+(c4z-h6z)**2 > 2.11**2) &
                 .or.                                             &
                ((c1x-h7x)**2+(c1y-h7y)**2+(c1z-h7z)**2 > 2.11**2 .and. &
                 (c2x-h7x)**2+(c2y-h7y)**2+(c2z-h7z)**2 > 2.11**2 .and. &
                 (c3x-h7x)**2+(c3y-h7y)**2+(c3z-h7z)**2 > 2.11**2 .and. &
                 (c4x-h7x)**2+(c4y-h7y)**2+(c4z-h7z)**2 > 2.11**2) &
                 .or.                                             &
                ((c1x-h8x)**2+(c1y-h8y)**2+(c1z-h8z)**2 > 2.11**2 .and. &
                 (c2x-h8x)**2+(c2y-h8y)**2+(c2z-h8z)**2 > 2.11**2 .and. &
                 (c3x-h8x)**2+(c3y-h8y)**2+(c3z-h8z)**2 > 2.11**2 .and. &
                 (c4x-h8x)**2+(c4y-h8y)**2+(c4z-h8z)**2 > 2.11**2) &
                 .or.                                             &
                ((c1x-h9x)**2+(c1y-h9y)**2+(c1z-h9z)**2 > 2.11**2 .and. &
                 (c2x-h9x)**2+(c2y-h9y)**2+(c2z-h9z)**2 > 2.11**2 .and. &
                 (c3x-h9x)**2+(c3y-h9y)**2+(c3z-h9z)**2 > 2.11**2 .and. &
                 (c4x-h9x)**2+(c4y-h9y)**2+(c4z-h9z)**2 > 2.11**2) &
                 .or.                                             &
                ((c1x-h10x)**2+(c1y-h10y)**2+(c1z-h10z)**2 > 2.11**2 .and. &
                 (c2x-h10x)**2+(c2y-h10y)**2+(c2z-h10z)**2 > 2.11**2 .and. &
                 (c3x-h10x)**2+(c3y-h10y)**2+(c3z-h10z)**2 > 2.11**2 .and. &
                 (c4x-h10x)**2+(c4y-h10y)**2+(c4z-h10z)**2 > 2.11**2) &
                  .or.                                          &
                ((c2x-c1x)**2+(c2y-c1y)**2+(c2z-c1z)**2 > 2.5**2 .and. &
                 (c3x-c1x)**2+(c3y-c1y)**2+(c3z-c1z)**2 > 2.5**2 .and. &
                 (c4x-c1x)**2+(c4y-c1y)**2+(c4z-c1z)**2 > 2.5**2) &
                 .or.                                             &
                ((c3x-c2x)**2+(c3y-c2y)**2+(c3z-c2z)**2 > 2.5**2 .and. &
                 (c4x-c2x)**2+(c4y-c2y)**2+(c4z-c2z)**2 > 2.5**2) &
                 .or.                                            &
                ((c4x-c3x)**2+(c4y-c3y)**2+(c4z-c3z)**2 > 2.5**2) &
                  ) then

!(write something to a file for submit.sh to read)
open(124,file="cancel.txt")
!write(128,*) "cancelling, current:" , current
write(124,*) "y"
close(124)
else

         open(unit=24, file="tgeom.xyz")
          write(24,*) "10"
          write(24,*) ""
          write(24,*)"C" , " ", C1X, " ", C1Y, " ", C1Z
          write(24,*)"C" , " ", C2X, " ", C2Y, " ", C2Z
          write(24,*)"C" , " ", C3X, " ", C3Y, " ", C3Z
          WRITE(24,*)"C" , " ", C4X, " ", C4Y, " ", C4Z
          WRITE(24,*)"H" , " ", H5X, " ", H5Y, " ", H5Z
          write(24,*)"H" , " ", H6X, " ", H6Y, " ", H6Z
          write(24,*)"H" , " ", H7X, " ", H7Y, " ", H7Z
          write(24,*)"H" , " ", H8X, " ", H8Y, " ", H8Z
          WRITE(24,*)"H" , " ", H9X, " ", H9Y, " ", H9Z
          WRITE(24,*)"H" , " ", H10X, " ", H10Y, " ", H10Z
         close(24)

!(write something else)          
open(124,file="cancel.txt")
write(124,*) "n"
close(124)

end if
end program
