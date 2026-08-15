program main
implicit none
integer*1, dimension(14):: qq,qq2, mcPoint
real, dimension(14):: qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real::start,finish,time,dumdu,now,addIt,qShrink, acceptance
real*8:: Ecut0, vbar,ein, ein2, vbarTMP,se,pot, max_value
integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y,xx,same,energyChecked, rand, oldRand
integer:: gradLine, termCount, keepStoring, rand2, oldRand2, yyy, tsign, same2, initi, signed
logical::there
character (len=1024) :: text, text2, text3, text4, text5, text6
character (len=1024) :: xxx,xxy
real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup
integer:: current, shellCounter, waiter,unitt,ax, dun, countn
character (len=13) ::measure, foundp, shell,foundling,preFoundp

shell="../shelll.txt"
measure="../measur.txt"
foundp="../foundp.txt"
foundling="foundling.txt"

vbartmp = 9898.0
countn = 0
max_value = -999999.0
open(unit=167, file="values.txt", status='unknown', action='read')
do
  read(167,*,iostat=irr) pot
  if (irr /= 0) exit
  countn = countn + 1
  if (pot > max_value) max_value = pot
end do
close(167)
if (countn == 3) then  !!change here
  vbartmp = max_value
endif

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
