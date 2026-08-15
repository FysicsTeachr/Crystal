program main
implicit none
integer*1, dimension(24):: qq,qq2, mcPoint
real, dimension(24):: qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real*8,dimension(:),allocatable:: pot, pot2
real::start,finish,time,dumdu,now,addIt,qShrink, acceptance
real*8:: Ecut0, vbar, ein, ein2, vbarTMP,se
integer, dimension(:), allocatable ::shellArr

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
character(len=1):: cancl

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

open(124, file="cancel.txt")
read(124,*) cancl
if(cancl=="y") then
  vbarTMP=99999
  se=99998
else
            open(unit=14, file="ttc.out")

            do

              read (14,"(a)",iostat=irr) text ! Read line into character variable

              if(index(text, "1  singlet").ne.0) then
                first= index(text, "-")
                xxx= text(first:first+13)
              end if

              if(index(text, "2  singlet").ne.0) then 
                energyChecked=1
                first= index(text, "-")
                xxy= text(first:first+13)
                close(14)
                exit
              end if

!              if(index(text, "Job finished").ne.0 .and. same==1) then
!                 energyChecked=1
!                 close(14)
!                 exit
!              end if

              if(index(text, "Job terminated").ne.0) then
                energyChecked=3
                exit
              end if

            end do

            if (energyChecked==3) then
              ein=99999.0
              ein2=99999.0+10.0
              vbarTMP=(ein2+ein)/2.0
            end if

            if (energyChecked==1) then
              read(xxx,*) ein                   ! Convert xxx into integer
              read(xxy,*)  ein2
              vbarTMP=(ein2+ein)/2.0
              oldrand=rand
            end if
            close(14)
            vbarTMP=(ein2+ein)/2.0
            se= Abs(ein2-ein)

end if

inquire(file=foundling, exist=there)
if (there) then
            open(12, file=foundling, status="old", position="append", action="write")
else
            open(12, file=foundling)
end if

  open(132,file="coord.txt")
  read(132,*) current, shellCounter
  read(132,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
  read(132,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
  read(132,*) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
  read(132,*) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
  read(132,*) qq2(21), qq2(22), qq2(23), qq2(24)
write(132,*) "writing to foundling"
write(12,*) current, vbarTMP, se
write(12,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
write(12,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
write(12,*) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
write(12,*) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
write(12,*) qq2(21), qq2(22), qq2(23), qq2(24)
close(12)

!open(12,file="../shellCo.txt")
!read(12,*) shellCounter
!allocate(shellArr(shellCounter))
!do i=1,shellCounter
!read(12,*) shellArr(i) 
!end do
!
!rewind(12)
!
!shellArr(current)=1
!write(12,*) shellCounter
!do i=1,shellCounter
!write(12,*) shellArr(i)
!end do
!close(12)

!do
!            open(10, file=measure)
!            read(10,*,iostat=irr) current, shellCounter, from, dun
!            if(irr/=0) then
!              close(10)
!              cycle
!            end if
!            rewind(10)    
!            write(10,*) current, shellCounter, from, dun-1
!            close(10)
!            exit
!end do  
          

end program
