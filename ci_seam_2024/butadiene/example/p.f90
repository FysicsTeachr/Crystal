!I am setting a tolerance in the sorting because there seem to be duplicates otherwise. Can be 
!fixed by setting gog's and qq's as integers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
implicit none
integer*1, dimension (:,:),allocatable:: gog,nogo,arrae,maygo, nogoold
integer*1, dimension(:),allocatable:: qq,qq2,qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric,leMay,riMay
real*8,dimension(:),allocatable:: pot, pot2, pot3, seampot, seampot2,seampot3, pot2old, seampot2old
real::start,finish,time,dumdu,now,addIt, acceptance
real*8:: Ecut0, vbar, ein, ein2, vbarTMP, se, ecutMax, incr
logical:: there
integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij,current,unique2
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y
integer:: xx,same,energyChecked, rand, oldRand,numbr, ibar2old, sCount, b
integer:: termCount,keepStoring,rand2,oldRand2,yyy,tsign,same2,initi,shellCounter,ibar3,rae,iba,dun

character::from
character (len=11) :: numnum(3)
character (len=10) ::set1, set2, set1b, set2b, shell, measure, foundp
character(len=24):: foundling1, foundling2, foundling3, foundling4, foundling5, foundling6, foundling7, foundling8, foundling9,foundling0
character(len=24):: foundlingq, foundlingw, foundlinge, foundlingr, foundlingt, foundlingy, foundlingu, foundlingi,foundlingo,foundlingp
character(len=24):: foundlinga, foundlings, foundlingd, foundlingf, foundlingg, foundlingh, foundlingj,foundlingk,foundlingl,foundlingz
character (len=1024) :: text
character (len=1024) :: text2, text3, text4, text5, text6
character (len=1024) :: xxx,xxy
character(len=1024):: reader, filename, format_string

real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup, nscal, seamCut

real, dimension(8,123):: G
!real, dimension(8):: qStretch

write(*,*) "started"

open(unit=128,file="converted.txt")

scal=1.0
d=24  !!--------> dimensions 

allocate(symmetric(d))
symmetric=0
!symmetric(3) =1  !!--------> Set = 1 for dimensions that are symmetric. NOTE: Fix this

set1="set002.txt"

set1b="set002.txt"
shell="shelll.txt"
measure="measur.txt"
foundp="foundp.txt"
foundling1="secondary1/foundling.txt"
foundling2="secondary2/foundling.txt"
foundling3="secondary3/foundling.txt"
foundling4="secondary4/foundling.txt"
foundling5="secondary5/foundling.txt"
foundling6="secondary6/foundling.txt"
foundling7="secondary7/foundling.txt"
foundling8="secondary8/foundling.txt"
foundling9="secondary9/foundling.txt"
foundling0="secondary0/foundling.txt"
foundlingq="secondaryq/foundling.txt"
foundlingw="secondaryw/foundling.txt"
foundlinge="secondarye/foundling.txt"
foundlingr="secondaryr/foundling.txt"
foundlingt="secondaryt/foundling.txt"
foundlingy="secondaryy/foundling.txt"
foundlingu="secondaryu/foundling.txt"
foundlingi="secondaryi/foundling.txt"
foundlingo="secondaryo/foundling.txt"
foundlingp="secondaryp/foundling.txt"
foundlinga="secondarya/foundling.txt"
foundlings="secondarys/foundling.txt"
foundlingd="secondaryd/foundling.txt"
foundlingf="secondaryf/foundling.txt"
foundlingg="secondaryg/foundling.txt"
foundlingh="secondaryh/foundling.txt"
foundlingj="secondaryj/foundling.txt"
foundlingk="secondaryk/foundling.txt"
foundlingl="secondaryl/foundling.txt"
foundlingz="secondaryz/foundling.txt"



shellcounter=0

seamCut=0.0165

Open(UNIT = 12, file=measure)      
write(12,*) 0
Close(12)

dump=0

msize=9999999  ! max basis size (for allocating memory).
Ecut0= -154.706     !!! cutoff energy for accepting lattices.
incr=0.004
EcutMax=-153
nscal=1.0      ! lattice spacing."

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

open(unit=180, file="radiux.txt") !    

Allocate(rad(d))                  !                      

ct=d*2+1

rad=0
write(180,*)rad
do iii=1,d
    rad=0
    do i=1,d
      if (i==iii) then
      rad(iii)=1
      end if
    end do
    write(180,*)rad
end do

do iii=1,d
  if(symmetric(iii)==0) then         !!!!symmetry code needs correction
    rad=0
    do i=1,d
      if (i==iii) then
      rad(iii)=1
      end if
    end do
      write(180,*)rad*(-1)
  else
    ct=ct-1
  end if
end do
rewind(180)
!write(*,*) "a"

Allocate(arrae(ct,d))     !Reading radiux.txt
do i=1,ct
read(180,*) (arrae(i,j),j=1,d)
end do

allocate (gog(d,msize))
allocate (nogo(d,msize))
!allocate (maygo(d,msize))
allocate (le(msize))
allocate (ri(msize))
allocate (leNo(msize))
allocate (riNo(msize))
allocate(qq(d))
allocate(qq2(d))
allocate(qt(d))
allocate(pot(msize))
allocate(seampot(msize))
allocate(seampot2(msize))
allocate(pot2(msize))
!write(*,*) "aa"
ibar=1
ibar2=1
le=-1
ri=-1
leNo=-1
riNo=-1
same=2
energyChecked=2
text="abc"
initi=1   !2 if center needs to be added
oldRand=0
oldRand2=0

!Defining the first point
do i=1,d
 gog(i,1)=0
 nogo(i,1)=125
end do


INQUIRE( FILE=set1, EXIST=THERE )
IF ( THERE ) then
open(unit=12, file=set1)
read(12,*)

do 
  read(12,*,iostat=irr) rae, iba
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if(irr/=0 .or. iba>1000) then
  if(irr/=0 ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    close(12)
    exit
  end if

  read(12,*) qq2(1), qq2(2),qq2(3), qq2(4), qq2(5), qq2(6), qq2(7), qq2(8),qq2(9),qq2(10), qq2(11), qq2(12), qq2(13),qq2(14),qq2(15), qq2(16)
  read(12,*) qq2(17), qq2(18),qq2(19),qq2(20), qq2(21), qq2(22), qq2(23),qq2(24)
  read(12,*) vbarTMP, se

  if(rae==2) then
    call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)    
    if (unique.ne.0) then
          if (unique==1) then
            call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbarTMP,msize,seampot2,se)
          else if (unique==-1) then
            call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbarTMP,msize,seampot2,se)
          end if    
    end if

  else if(rae==1)then
    call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
    if (unique.ne.0) then
          if (unique==1) then
            call addTo(ibar,ri,gog,iii,qq2,d,pot,vbarTMP,msize,seampot,se)
          else if (unique==-1) then
            call addTo(ibar,le,gog,iii,qq2,d,pot,vbarTMP,msize,seampot,se)
          end if
    end if

  end if
end do
end if

INQUIRE( FILE=foundp, EXIST=THERE ) 
IF ( THERE ) then

open(unit=12, file=foundp)
read(12,*)
do
  read(12,*,iostat=irr) rae ,vbarTMP, se
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(irr/=0 ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    close(12)
    exit
  end if
write(*,*) "foundp"
  read(12,*) qq2(1), qq2(2),qq2(3), qq2(4), qq2(5)
  read(12,*) qq2(6), qq2(7), qq2(8),qq2(9),qq2(10)
  read(12,*) qq2(11), qq2(12), qq2(13),qq2(14),qq2(15)
  read(12,*) qq2(16), qq2(17), qq2(18),qq2(19),qq2(20)
  read(12,*) qq2(21), qq2(22), qq2(23),qq2(24)


  if(vbarTMP>Ecut0 .or. Abs(se)>seamCut) then
    call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)
    if (unique.ne.0) then
          if (unique==1) then
            call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbarTMP,msize,seampot2,se)
!write(*,*) "f rejected"
          else if (unique==-1) then
            call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbarTMP,msize,seampot2,se)
!write(*,*) "f rejected"
          end if
    end if

  else if(vbarTMP<Ecut0 .and. Abs(se)<seamCut)then
    call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
    if (unique.ne.0) then
          if (unique==1) then
            call addTo(ibar,ri,gog,iii,qq2,d,pot,vbarTMP,msize,seampot,se)
!write(*,*) "f accected"
          else if (unique==-1) then
            call addTo(ibar,le,gog,iii,qq2,d,pot,vbarTMP,msize,seampot,se)
!write(*,*) "f accected"
          end if
    end if
  end if
  
end do

END IF





write(*,*) "ibar, ibar2", ibar, ibar2

 Open(UNIT = 12, file=set1)
 write(12,*)
 Close(12)


 open(12, file=set1b, status="old", position="append", action="write")
 do i=2,ibar2
   write(12,*) 2,i
   write(12,*) nogo(1,i), nogo(2,i),nogo(3,i),nogo(4,i), nogo(5,i), nogo(6,i), nogo(7,i),nogo(8,i),nogo(9,i), nogo(10,i), nogo(11,i), nogo(12,i),nogo(13,i),nogo(14,i), nogo(15,i),nogo(16,i)
   write(12,*)  nogo(17,i),nogo(18,i),nogo(19,i), nogo(20,i), nogo(21,i), nogo(22,i),nogo(23,i),nogo(24,i)
   write(12,*)  pot2(i), seampot2(i)
 end do
 do i=2,ibar
   write(12,*) 1,i
   write(12,*) gog(1,i), gog(2,i), gog(3,i), gog(4,i), gog(5,i), gog(6,i), gog(7,i), gog(8,i), gog(9,i), gog(10,i), gog(11,i), gog(12,i), gog(13,i), gog(14,i), gog(15,i), gog(16,i)
   write(12,*) gog(17,i), gog(18,i), gog(19,i), gog(20,i), gog(21,i), gog(22,i), gog(23,i), gog(24,i)
   write(12,*)  pot(i), seampot(i)
 end do
 close(12)

write(*,*) ibar, ibar2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do
  
  ecut0=ecut0+incr
 write(*,*) "ecut", ecut0
  if (ecut0>ecutMax)then
    write(*,*) "increase ecutMax"
    exit
  end if
  
  allocate(nogoOld(d,msize))
  allocate(pot2old(msize))
  allocate(seampot2old(msize))
  nogoOld=nogo
  ibar2old=ibar2
  pot2old=pot2
  seampot2old=seampot2

  deallocate(nogo)
  deallocate(pot2)
  deallocate(seampot2)
  ibar2=1
  allocate(nogo(d,msize))
  allocate(pot2(msize))
  allocate(seampot2(msize))
  do i=1,d
    nogo(i,1)=125
  end do
  leno=-1
  rino=-1

  do i=2, ibar2old
  vbarTMP=pot2old(i)
  se=seampot2old(i)

    do chki=1,d
      qq2(chki)=nogoOld(chki,i)
    end do
  
    if(vbarTMP<=Ecut0 .AND. Abs(se)<=seamcut) then

        call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
      if (unique==1) then
        call addTo(ibar,ri,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
      else if (unique==-1) then
        call addTo(ibar,le,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
      end if
    write(*,*) "transfer ibar2 to ibar", i, ibar
    write(128,*) i, ibar
    write(128,*) qq2
    write(128,*) vbarTMP, se

    else if(vbarTMP>eCUT0 .or. Abs(se)>seamcut) then
        call sort(qq2,d,nogo,leno,rino,msize,i,unique,iii)
      if (unique==1) then
        call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
      else if (unique==-1) then
        call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
      end if
    end if
  end do

  deallocate(nogoOld)
  deallocate(pot2old)
  deallocate(seampot2old)

  write(*,*) "shuffle completed, ibar , ibar2, ibar2old"
  write(*,*) ibar, ibar2, ibar2old


ci=1
cf=ibar

!write(*,*) "aaa"

do while (cf>=ci)

 if (ibar>msize/8) then !To avoid unwanted symmetry breaks
  write(4,*) "killing code due to the msize/5 allocation", ibar
  write(*,*) "killing code due to the msize/5 allocation", ibar
  dump=1
  exit
 end if


 ibar3=1
 shellCounter=0
 allocate(maygo(d,msize))
 allocate(pot3(msize))
 allocate(seampot3(msize))
 allocate (leMay(msize))
 allocate (riMay(msize))
 leMay=-1
 riMay=-1


 do i=1,d
  mayGo(i,1)=125
 end do

! write(*,*) "bbb"

 do i=ci,cf
 
  do chki=1,d
   qq(chki)=gog(chki,i)
  end do
  
  do t=1,d
   do tsign = -1,1,2
   
    qq2=qq
    qq2(t)=qq2(t)+tsign*1.0


    !write(*,*) qq2
    !if(qq2(4)>=0 .and. qq2(4)<4.6 .and.Abs(qq2(5))<4.6.and.Abs(qq2(6))<4.6.and.Abs(qq2(7))<4.6.and.Abs(qq2(8))<4.6 .and. qq2(1)>-1.5   .and.   qq2(2)>-1.9  .and. qq2(3) >-1.9 )then
            
      call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)
      if (unique.ne.0) then
      call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
      if (unique.ne.0) then
          call sort(qq2,d,maygo,leMay,riMay,msize,i,unique,iii)
          if (unique==1) then
            call addTo(ibar3,riMay,maygo,iii,qq2,d,pot3,0,msize,seampot3,0)
          else if (unique==-1) then
            call addTo(ibar3,leMay,maygo,iii,qq2,d,pot3,0,msize,seampot3,0)
          end if
      end if
      end if
    !end if
   end do
  end do
 end do
 
! write(*,*) "ccc", ibar3

 shellCounter=ibar3

write(*,*) "ibar3", ibar3

if(ibar3==1) then
 deallocate(maygo)
 deallocate(pot3)
 deallocate(seampot3)
 deallocate (leMay)
 deallocate (riMay)
exit
end if

if(ibar3>1) then

 open(12, file=shell)
 do i=2,ibar3 
     write(12,*) mayGo(1,i),mayGo(2,i),mayGo(3,i),mayGo(4,i),mayGo(5,i)
     write(12,*) mayGo(6,i),mayGo(7,i),mayGo(8,i),maygo(9,i),maygo(10,i)
     write(12,*) mayGo(11,i),mayGo(12,i),mayGo(13,i),mayGo(14,i),mayGo(15,i)
     write(12,*) mayGo(16,i),mayGo(17,i),mayGo(18,i),maygo(19,i),maygo(20,i)
     write(12,*) mayGo(21,i),mayGo(22,i),mayGo(23,i),mayGo(24,i)

 end do
 close(12) 

 open(12, file=measure)
 write(12,*) 1, ibar3-1, "p",0
 close(12)

! open(12,file="shellCo.txt")
! write(12,*) ibar3-1
! do i=2, ibar3
!   write(12,*) 0
! end do
! close(12)
 
 do
  do
    open(12, file=measure)
    read(12,*,iostat=irr) current, yyy, from,dun
      if(irr/=0) then
        close(12)
        cycle
      end if
    close(12)
    exit
  end do
  
  !scount=0
  !open(12,file="shellCo.txt")
  !read(12,*)
  !do i=1, ibar3-1
  !  read(12,*) dun
  !  scount=scount+dun
  !end do

  if(current==0 .and. from.eq."s") then! .and. scount==ibar3-1) then
    exit
  else
    call sleep(7)
  end if
 end do

! write(*,*) "waiting 100 ms for subjob completion"
! call sleepqq(100)          !time for individual nodes to finish calculations.
 write(*,*) ibar, ibar2

 open(14, file=foundp)
 do b=1,25
 if (b < 10) then
      format_string = "(A9,I1,A14)"
    else if (b<100) then
      format_string = "(A9,I2,A14)"
 end if
 write(filename,format_string) "secondary",b,"/foundling.txt"
 open(13, file=filename)
 do
    read(13,*,iostat=irr) numbr, vbarTMP, se
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) numbr, vbarTMP, se
    read(13,*,iostat=irr) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
    read(13,*,iostat=irr) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
    read(13,*,iostat=irr) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
    read(13,*,iostat=irr) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
    read(13,*,iostat=irr) qq2(21), qq2(22), qq2(23), qq2(24)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(21), qq2(22), qq2(23), qq2(24)
 end do
 end do

 close(14)
 open(14, file=foundp)
 do
   read(14,*,iostat=irr) numbr, vbarTMP, se
   IF (irr/=0) then
     close(14)
     EXIT
   end if

   read(14,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
   read(14,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
   read(14,*) qq2(11), qq2(12), qq2(13), qq2(14), qq2(15)
   read(14,*) qq2(16), qq2(17), qq2(18), qq2(19), qq2(20)
   read(14,*) qq2(21), qq2(22), qq2(23), qq2(24)

   
   if(vbarTMP.le.Ecut0 .and. Abs(se).le.seamCut) then
     call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
     if (unique==1) then
       call addTo(ibar,ri,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
     else if (unique==-1) then
       call addTo(ibar,le,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
     end if
     open(12, file=set1b, status="old", position="append", action="write")
     write(12,*) 1,ibar
     write(12,*) qq2
     write(12,*)  vbarTMP, se
     close(12)

   else if (vbarTMP>Ecut0 .or. Abs(se)>seamCut) then
     call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)
     if (unique==1) then
       call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
     else if (unique==-1) then
       call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
     end if
     open(12, file=set1b, status="old", position="append", action="write")
     write(12,*) 2,ibar2
     write(12,*) qq2
     write(12,*)  vbarTMP, se
     close(12)
   end if
 end do

 write(*,*) "shell finished, ibar=" ,ibar

 Open(UNIT = 12, file=shell)      
 write(12,*)
 Close(12)
 Open(UNIT = 12, file=measure)
 write(12,*) 0, ibar3-1, "p",0
 Close(12)
 Open(UNIT = 12, file=foundp)
 write(12,*)
 Close(12)

end if

 deallocate (leMay)
 deallocate (riMay)

 deallocate(maygo)
 deallocate(pot3)
 deallocate(seampot3)
 ci=cf+1
 cf=ibar

end do

end do

write(*,*) "job completed"

end program


function grth(d,qq,i,iii,msize,gog)                     !returns 1 if greater, else 0
integer:: d,i,iii,a,grth,z,xx,b,msize
integer*1, dimension(d,msize):: gog
integer*1, dimension(d):: qq
xx=0
do z=1,d
  if(qq(z)>gog(z,iii))then
    xx=1
    exit
  else if(qq(z)<gog(z,iii))then
    exit
  end if
end do
grth=xx
return
end function


function leth(d,qq,i,iii,msize,gog)                     !returns 1 if lesser, else 0
integer:: d,i,iii,a,leth,z,xx,b,msize
integer*1, dimension(d,msize):: gog
integer*1, dimension(d):: qq
xx=0
do z=1,d
  if(qq(z)<gog(z,iii))then
    xx=1
    exit
  else if(qq(z)>gog(z,iii))then
    exit
  end if
end do
leth=xx
return
end function


subroutine sort(qq,d,gogo,le,ri,msize,i,unique,iii) !unique is (0)/(1)/(-1) if the vector is
!(already in list)/(on right side of lowest element in the binary-sort tree)/(on left side)
integer:: i,d,msize,iii,unique,grth,leth
integer*1,dimension(d)::q,qq
integer*1,dimension(d,msize)::gogo
integer, dimension(msize):: le,ri
         unique=0
         iii=1
         do
             if ((grth(d,qq,i,iii,msize,gogo))==1) then !If (the vector being sorted is greater
                                                      !than the iii-th in the binary tree) then
                 if (ri(iii)>0) then  !If (there is an element to the right of iii)
                   iii = ri(iii) !Then, the next time when this the loop repeats,
                     ! compare the vector with the element on the right side of iii
                 else
                   unique=1   !Else report to the main code that the vector is a new element. (iii and unique are the outputs of this section)
                 exit
                 end if
             else if (leth(d,qq,i,iii,msize,gogo)==1) then
                 if (le(iii)>0) then
                   iii = le(iii)
                 else
                   unique=-1
                 exit
                 end if
             else
                 exit    !Unique remains 0 if the vector is same as iii
             end if
         end do
end subroutine

subroutine addTo(iba,ri,gogo,iii,qq,d,pot,vbartmp,msize,seampot,se)
integer:: iba,iii,el,d
integer, dimension(msize):: ri
integer*1,dimension(d,msize)::gogo
integer*1,dimension(d)::qq
real*8,dimension(msize)::pot,seampot
real*8:: vbarTMP, se
iba=iba+1
ri(iii)=iba                                             !
do el=1,d
  gogo(el,iba)=qq(el)
end do
pot(iba)=vbartmp
seampot(iba)=se
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
