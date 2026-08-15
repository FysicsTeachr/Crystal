program main
implicit none
integer*1, dimension(14):: qq,qq2, mcPoint
real, dimension(14):: qt
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

real::  nscal, seamCut
real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup

real, dimension(8,123):: G
real, dimension(8):: gradAlong, coupAlong, qStretch
real, dimension(:,:), allocatable:: gradMatrix, coupMatrix

integer:: current, shellCounter, waiter,unitt,ax, dun
character (len=13) ::measure, foundp, shell,foundling,preFoundp
character::from

waiter=2
oldRand=0
shell="../shelll.txt"
measure="../measur.txt"
foundp="../foundp.txt"
preFoundp="prefoundp.txt"
foundling="foundling.txt"

            open(unit=14, file="ttc.out")

            do

              read (14,"(a)",iostat=irr) text 
              if (irr /= 0) then
                energyChecked = 3
                exit
              end if


              if(index(text, "FINAL ENERGY:").ne.0) then
                first= index(text, "-")
                xxx= text(first:first+12)
              end if

              !FINAL ENERGY: -66.1162495623 a.u.
              if(index(text, "Job finished").ne.0) then
                 energyChecked=1
                 exit
              end if

              if(index(text, "Job terminated").ne.0) then
                energyChecked=3
                exit
              end if

              if(index(text, "DL-FIND ERROR:").ne.0) then
                energyChecked=3
                exit
              end if

            end do

            if (energyChecked==3) then
              open(134,file="cancel.txt")
              write(134,*)"y"
              close(134)
              ein=99999.0
            end if

            if (energyChecked==1) then
              read(xxx,*) ein       ! Convert xxx into integer
            end if
            close(14)
            vbarTMP=ein

!end if

open(132,file="coord.txt")
read(132,*) current, shellCounter
read(132,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
read(132,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
read(132,*) qq2(11), qq2(12), qq2(13), qq2(14)

inquire(file=foundling, exist=there)
if (there) then
            open(12, file=foundling, status="old", position="append", action="write")
else
            open(12, file=foundling)
end if
  
write(12,*) current, vbarTMP
write(12,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
write(12,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
write(12,*) qq2(11), qq2(12), qq2(13), qq2(14)
close(12)

end program
