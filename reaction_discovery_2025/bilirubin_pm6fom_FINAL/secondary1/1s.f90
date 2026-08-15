program main
implicit none
integer*1, dimension(14):: qq,qq2, mcPoint
real, dimension(14):: qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real*8,dimension(:),allocatable:: pot, pot2
real::start,finish,time,dumdu,now,addIt,qShrink, acceptance,dum, midx,midy,midz
real*8:: Ecut0, vbar, ein, ein2, vbarTMP,se
real, dimension(79,3)::arr
real::lo, low1,low2,low3,lowd,low1i,low2i,low3i,lowdi
real,dimension(3)::vec1,vec2,vec3,vec

integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij
integer:: msize,p,s,chki,ibar2,dump,grth,leth,el,t,xx,same,energyChecked, rand, oldRand
integer:: gradLine, termCount, keepStoring, rand2, oldRand2, yyy, tsign, same2, initi, signed, qq2d
integer:: pre
character (len=11) :: numnum(3)
character (len=10) ::output1
character (len=10) ::set1, set2, set1b, set2b
real, Dimension(123)::grad, coup
logical::there
character (len=1024) :: text, command
character (len=1024) :: text2, text3, text4, text5, text6
character (len=1024) :: xxx,xxy
character(len=3):: reader

real:: CX, CY, CZ, NX,NY,NZ, H1X,H1Y,H1Z, H2X, H2Y,H2Z, H3X,H3Y,H3Z, nscal, seamCut
real::CXp,CYp,CZp,NXp,NYp,NZp,H1Xp,H1Yp,H1Zp,H2Xp,H2Yp,H2Zp,H3Xp,H3Yp,H3Zp,oldH1xp,oldH1yp
real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup
real*4:: x,y,z

real::C1X,C1Y,C1Z,C2X,C2Y,C2Z,C3X,C3Y,C3Z,C4X,C4Y,C4Z,C5X,C5Y,C5Z,C6X,C6Y,C6Z
real::C7X,C7Y,C7Z,C8X,C8Y,C8Z,C9X,C9Y,C9Z,C10X,C10Y,C10Z,C11X,C11Y,C11Z,C12X,C12Y,C12Z
real::C13X,C13Y,C13Z,C14X,C14Y,C14Z,C15X,C15Y,C15Z,N16X,N16Y,N16Z,O17X,O17Y,O17Z,O18X,O18Y,O18Z
real::O19X,O19Y,O19Z,O20X,O20Y,O20Z,O21X,O21Y,O21Z,H22X,H22Y,H22Z,H23X,H23Y,H23Z,H24X,H24Y,H24Z
real::H25X,H25Y,H25Z,H26X,H26Y,H26Z,H27X,H27Y,H27Z,H28X,H28Y,H28Z,H29X,H29Y,H29Z,H30X,H30Y,H30Z
real::H31X,H31Y,H31Z,H32X,H32Y,H32Z,H33X,H33Y,H33Z,H34X,H34Y,H34Z,H35X,H35Y,H35Z,H36X,H36Y,H36Z
real::H37X,H37Y,H37Z,H38X,H38Y,H38Z,H39X,H39Y,H39Z,H40X,H40Y,H40Z,H41X,H41Y,H41Z,H42X,H42Y,H42Z


real, dimension(8,123):: G
real, dimension(8):: gradAlong, coupAlong, qStretch
real, dimension(:,:), allocatable:: gradMatrix, coupMatrix

integer:: current, shellCounter, waiter,unitt,ax, dun, ios
character (len=13) ::measure, foundp, shell,foundling,preshell, pre_str
character::from

waiter=2
oldRand=0
shell="../shelll.txt"
preshell="../preShl.txt"
measure="../measur.txt"
foundp="../foundp.txt"
foundling="foundling.txt"

  do
    open(10, file=measure)
    read(10,*,iostat=irr) current, shellCounter,from, dun
    if(irr/=0.or.current==0)then
      close(10)
      call sleep(1)
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

  open(10,file="status.txt")
  write(10,*) "r"
  close(10)

  open(10, file= shell)
  do i=1,current-1
    read(10,*)
    read(10,*)
    read(10,*)
  end do

  read(10,*,iostat=irr) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
  read(10,*,iostat=irr)  qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
  read(10,*,iostat=irr) qq2(11), qq2(12), qq2(13), qq2(14)
  close(10)

  open(132,file="coord.txt")
  write(132,*) current, shellCounter
  write(132,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
  write(132,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
  write(132,*) qq2(11), qq2(12), qq2(13), qq2(14)
  close(132)

  open(110,file=preshell)
  do i=1,current-1
      read(110,*)
  end do
  read(110,*) pre, qq2d
  close(110)

 !copy pre.xyz from ../prefolder and save as pregeom.xyz
 write(command, '("cp ../prefolder/", i0, ".xyz pregeom.xyz")') pre
 call execute_command_line(command)

   qq2=0
  if(qq2d>0) then
    qq2(Abs(qq2d))=1
  else
    qq2(Abs(qq2d))=-1
  end if

  qt = 90.0
  qt(11) =20.0
  qt(12)= 0.6
  qt(13) = 0.6
  qt(14) = 0.6
  qt=qq2*1.0d0*qt
  open(132,file="coord1.txt")
  write(132,*) current, shellCounter
  write(132,*) qt(1), qt(2), qt(3), qt(4), qt(5)
  write(132,*) qt(6), qt(7), qt(8), qt(9), qt(10)
  write(132,*) qt(11), qt(12), qt(13), qt(14)
  close(132)

!write(*,*) qq2

open(10,file="cancel.txt")
write(10,*) "n"
close(10)

open(10,file="proceed.txt")
write(10,*) "y"
close(10)

end program
