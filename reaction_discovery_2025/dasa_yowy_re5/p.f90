!coordinates 3, 11 (terminal angles) being additionally minimized in additio to 1
program main
implicit none
integer*1, dimension (:,:),allocatable:: gog,nogo,arrae,maygo,nogoold,may
integer*1, dimension(:),allocatable:: qq,qq2,qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric,leMay,riMay, mia, qq2dir
real*8,dimension(:),allocatable:: pot,pot2,pot3,seampot,seampot2,seampot3,pot2old, seampot2old 
real::start,finish,time,dumdu,now,addIt, acceptance
real*8:: Ecut0, vbar, ein, ein2, vbarTMP, se, ecutMax, incr
logical:: there
integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij,current,unique2
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y, ia, qq2d
integer:: xx,same,energyChecked, rand, oldRand,numbr, ibar2old, sCount, b
integer:: termCount,keepStoring,rand2,oldRand2,yyy,tsign,same2,initi,shellCounter,ibar3,rae,iba,dun

character::from,stat
character (len=11) :: numnum(3)
character (len=10) ::set1, set2, set1b, set2b, shell, measure, foundp
character (len=1024) :: text
character (len=1024) :: text2, text3, text4, text5, text6
character (len=1024) :: xxx,xxy
character (len=1024)::  filename, format_string, namefile

real:: scal,  sumGrad, sumCoup, dotGrad, dotCoup, nscal, seamCut

real, dimension(8,123):: G

character(len=100)::i_str, ibar_str, numbr_str, ibar2_str,command
!real, dimension(8):: qStretch
logical :: file_exists

real*4, dimension(:,:,:), allocatable::allfiles, yesfiles, nofiles, nofilesTMP
character(len=1) ::atom
integer::xi, bbb

bbb=399   !number of secondary jobs

allocate(yesfiles(100000,42,3))
allocate(nofiles(1000000,42,3))
allocate(nofilesTMP(1000000,42,3))
open(113,file="rejected_trace.txt")
open(unit=131, file="new_lattice.xyz")
open(unit=132, file="climbing.xyz")

!open(unit=128,file="1.txt")
!read(128,*)
!read(128,*)
!do i=1,42
!  read(128,*) atom, yesfiles(1,i,1), yesfiles(1,i,2), yesfiles(1,i,3)
!end do
!close(128)

write(*,*) "started"



open(unit=128,file="converted.txt")

scal=1.0
d=14  !!--------> dimensions 
msize=9999999

allocate(symmetric(d))
symmetric=0
!symmetric(3) =1  !!--------> Set = 1 for dimensions that are symmetric. NOTE: Fix this

set1="set002.txt"

set1b="set002.txt"
shell="shelll.txt"
measure="measur.txt"
foundp="foundp.txt"


shellcounter=0

seamCut=0.0165
se=0

Open(UNIT = 12, file=measure)      
write(12,*) 0
Close(12)

dump=0

Ecut0= -71.24     !!! cutoff energy for accepting lattices.
incr=0.008
EcutMax=-71.20
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

  read(12,*) qq2(1), qq2(2),qq2(3), qq2(4), qq2(5), qq2(6), qq2(7), qq2(8),qq2(9),qq2(10), qq2(11), qq2(12), qq2(13),qq2(14)
  read(12,*) vbarTMP

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
  read(12,*,iostat=irr) rae ,vbarTMP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(irr/=0 ) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    close(12)
    exit
  end if
write(*,*) "foundp"
  read(12,*) qq2(1), qq2(2),qq2(3), qq2(4), qq2(5)
  read(12,*) qq2(6), qq2(7), qq2(8),qq2(9),qq2(10)
  read(12,*) qq2(11), qq2(12), qq2(13),qq2(14)


  if(vbarTMP>Ecut0) then
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

  else if(vbarTMP<Ecut0)then
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
   write(12,*) nogo(1,i), nogo(2,i),nogo(3,i),nogo(4,i), nogo(5,i), nogo(6,i), nogo(7,i),nogo(8,i),nogo(9,i), nogo(10,i), nogo(11,i), nogo(12,i),nogo(13,i),nogo(14,i)
   write(12,*)  pot2(i)
 end do
 do i=2,ibar
   write(12,*) 1,i
   write(12,*) gog(1,i), gog(2,i), gog(3,i), gog(4,i), gog(5,i), gog(6,i), gog(7,i), gog(8,i), gog(9,i), gog(10,i), gog(11,i), gog(12,i), gog(13,i), gog(14,i)
   write(12,*)  pot(i)
 end do
 close(12)

write(*,*) ibar, ibar2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!open(197,file="filesinfo.txt")

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

    do chki=1,d
      qq2(chki)=nogoOld(chki,i)
    end do
  
    if(vbarTMP<=Ecut0) then

      call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
      if (unique==1) then
       call addTo(ibar,ri,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)

      else if (unique==-1) then
       call addTo(ibar,le,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
      end if

!      write(*,*) "transfer ibar2 to ibar", i, ibar
!      write(128,*) i, ibar
!      write(128,*) qq2
!      write(128,*) vbarTMP
!      do b=1,42
!        yesfiles(ibar,b,1) = nofiles(i-1,b,1)
!        yesfiles(ibar,b,2) = nofiles(i-1,b,2)
!        yesfiles(ibar,b,3) = nofiles(i-1,b,3)
!      end do

  ! 132: if vbarTMP<=Ecut0 then transfer file i-1.pun from nofolder to yesfolder and rename it as ibar.pun
  !Also copy it to climbing lattices

      write(command, '("mv nofolder/", i0, ".pun yesfolder/", i0, ".pun")') i-1, ibar
      call execute_command_line(command)

      write(command, '("cp yesfolder/", i0, ".pun climbinglat/")') ibar
      call execute_command_line(command)

!      write(132,*) 42
!      write(132,*) vbarTMP, ibar
!      do b=1,42
!        write(132,*) yesfiles(ibar,b,1), yesfiles(ibar,b,2), yesfiles(ibar,b,3)
!      end do


    else if(vbarTMP>eCUT0) then
      call sort(qq2,d,nogo,leno,rino,msize,i,unique,iii)
      if (unique==1) then
        call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
      else if (unique==-1) then
        call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
      end if

       write(command, '("cp nofolder/", i0, " nofoldertmp/", i0)') i-1, ibar2-1
       call execute_command_line(command)


!       do b=1,42
!         nofilesTMP(ibar2-1,b,1) = nofiles(i-1,b,1)
!         nofilesTMP(ibar2-1,b,2) = nofiles(i-1,b,2)
!         nofilesTMP(ibar2-1,b,3) = nofiles(i-1,b,3)
!       end do

    end if
  end do

  write(command, '("rm nofolder")')
  call execute_command_line(command)
  ! Rename nofoldertmp to nofolder
  write(command, '("mv nofoldertmp nofolder")')
  call execute_command_line(command)

  ! Create a new folder named nofoldertmp
  write(command, '("mkdir nofoldertmp")')
  call execute_command_line(command)


!  nofiles = nofilesTMP

!  deallocate(nofiles)
!  allocate(nofiles(ibar2-1,42,3))
!  do i=1,ibar2-1
!   do t=1, 42
!    do b=1,3
!    nofiles(i,t,b)=nofilesTMP(i,t,b)
!    end do
!   end do
!  end do
!  deallocate(nofilesTMP)

!  write(197,*) 'nofiles updated'
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
 allocate(may(d,msize))
 allocate(mia(msize))
 allocate(qq2dir(msize))
 allocate(pot3(msize))
 allocate(seampot3(msize))
 allocate (leMay(msize))
 allocate (riMay(msize))
 leMay=-1
 riMay=-1

 do i=1,d
  mayGo(i,1)=125
  may(i,1) = 125
 end do

 mia(1) =999999
 qq2dir(1)=999997

 do i=ci,cf
 
  do chki=1,d
   qq(chki)=gog(chki,i)
  end do
  
  do t=1,d
   do tsign = -1,1,2
    qq2d=tsign*t
    qq2=qq
    qq2(t)=qq2(t)+tsign*1.0

    if(qq2(2).le.4 .and. qq2(4).le.4 .and. qq2(6).le.4 .and. qq2(8).le.4 .and. qq2(10).le.0  &
    .and. qq2(2)>-4 .and. qq2(4)>-4 .and. qq2(6)>-4 .and. qq2(8)>-4 .and. qq2(10)>-1 &
    .and. qq2(1)<1 .and. qq2(3)<1 .and. qq2(5)<1 .and. qq2(7)<1 &
    .and. Abs(qq2(1))<2 .and. Abs(qq2(3))<2 .and. Abs(qq2(5))<2.and. Abs(qq2(7))<2.and. Abs(qq2(9))<1  &
    .and. Abs(qq2(11))<1 .and.  Abs(qq2(12))<1 .and. Abs(qq2(13))<1 .and. Abs(qq2(14))<1  )then

      call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)
      if (unique.ne.0) then
        call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
        if (unique.ne.0) then
          call sort(qq2,d,maygo,leMay,riMay,msize,i,unique,iii)
          if (unique==1) then
            call addTo(ibar3,riMay,maygo,iii,qq2,d,pot3,0,msize,seampot3,0)
            do ia=1,d
              may(ia,ibar3) = qq(ia)
            end do
              mia(ibar3) = i
              qq2dir(ibar3) = qq2d
          else if (unique==-1) then
            call addTo(ibar3,leMay,maygo,iii,qq2,d,pot3,0,msize,seampot3,0)
            do ia=1,d
              may(ia,ibar3) = qq(ia)
            end do
              mia(ibar3) = i
              qq2dir(ibar3) = qq2d
          end if
        end if
      end if
    end if
   end do
  end do
 end do
 
! write(*,*) "ccc", ibar3

 shellCounter=ibar3

write(*,*) "ibar3", ibar3

if(ibar3==1) then
 deallocate(maygo)
 deallocate(mia)
 deallocate(qq2dir)
 deallocate(may)
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
     write(12,*) mayGo(11,i),mayGo(12,i),mayGo(13,i),mayGo(14,i)
 end do
 close(12) 


 !remove all files from prefolder
 write(command, '("rm -f prefolder/*")')
 call execute_command_line(command)

 open(12, file="preShl.txt")
 do i=2,ibar3
   write(12,*) mia(i), qq2dir(i)

   !copy mia(i).pun from yesfolder to prefolder

   !speedup idea: can make an array of mia(i) values and only use the below command when there is a new mia(i) to avoid repeatition

   write(command, '("cp yesfolder/", i0, ".pun prefolder/")') mia(i)
   call execute_command_line(command)

 end do
 close(12)

 open(12, file=measure)
 write(12,*) 1, ibar3-1, "p",0
 close(12)
 
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

  if(current==0 .and. from.eq."s") then! .and. scount==ibar3-1) then
    exit
  else
    call sleep(1)
  end if
 end do

! call sleep(20)


 do b=1,bbb              
   if (b < 10) then
     format_string = "(A9,I1,A11)"
   else if (b<100) then
     format_string = "(A9,I2,A11)"
   else
     format_string = "(A9,I3,A11)"
   end if
     write(namefile,format_string) "secondary",b,"/status.txt"
     do
       open(13, file=namefile)
       read(13,*,iostat=irr) stat
       if(irr/=0) then
        close(13)
        cycle
       end if
       close(13)
       !tcount=tcount+1
       if(trim(adjustl(stat))=="d") then
         exit
       end if
       call sleep(1)
     end do
 end do


 write(*,*) ibar, ibar2

 open(14, file=foundp)
 do b=1,bbb
 if (b < 10) then
      format_string = "(A9,I1,A14)"
    else if (b<100) then
      format_string = "(A9,I2,A14)"
    else
      format_string = "(A9,I3,A14)"
 end if
 write(filename,format_string) "secondary",b,"/foundling.txt"
 open(13, file=filename)
 do
    read(13,*,iostat=irr) numbr, vbarTMP
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) numbr, vbarTMP
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
    read(13,*,iostat=irr) qq2(11), qq2(12), qq2(13), qq2(14)
     IF (irr/=0) then
       close(13, status="delete")
       EXIT
     end if
    write(14,*) qq2(11), qq2(12), qq2(13), qq2(14)

 end do
 end do

 close(14)

write(*,*) "check1- foundp made"


! write(*,*)"opening allfiles"
! open(14,file="allfiles.txt")
! write(*,*) "opened allfiles"
! allocate(allfiles(ibar3-1,42,3))  !!!Might need to be ibar if bugs
! allfiles=9999.0
! do
!   read(14,*,iostat=irr) xi
!   if(irr/=0)then
!     close(14)
!     exit
!   end if
!   read(14,*) !to ignore line with 42
!   read(14,*) !ignore energy line
!   do b=1,42
!     read(14,*) atom,allfiles(xi,b,1), allfiles(xi,b,2), allfiles(xi,b,3)
!   end do
! end do
! call execute_command_line('rm allfiles.txt')

! write(*,*) "allfiles read complete"

 open(14, file=foundp)
 do
   read(14,*,iostat=irr) numbr, vbarTMP
   IF (irr/=0) then
     close(14)
     EXIT
   end if

! Additional debug statements to track file position

   read(14,*) qq2(1), qq2(2), qq2(3), qq2(4), qq2(5)
   read(14,*) qq2(6), qq2(7), qq2(8), qq2(9), qq2(10)
   read(14,*) qq2(11), qq2(12), qq2(13), qq2(14)

   if(vbarTMP.le.Ecut0) then
     call sort(qq2,d,gog,le,ri,msize,i,unique,iii)
     if (unique==1) then
       call addTo(ibar,ri,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)

     else if (unique==-1) then
       call addTo(ibar,le,gog,iii,qq2,d,pot,vbartmp,msize,seampot,se)
     end if
     if(unique.ne.0) then
       open(12, file=set1b, status="old", position="append", action="write")
       write(12,*) 1,ibar
       write(12,*) qq2
       write(12,*)  vbarTMP
       close(12)

       !transfer numbr.pun from allfolder to yesfolder and rename it as ibar.pun
       write(command, '("mv allfolder/", i0, ".pun yesfolder/", i0, ".pun")') numbr, ibar
       call execute_command_line(command)

!       if(allfiles(numbr,1,1)<9998.0) then
!         do b=1,42
!           yesfiles(ibar,b,1) = allfiles(numbr,b,1)
!           yesfiles(ibar,b,2) = allfiles(numbr,b,2)
!           yesfiles(ibar,b,3) = allfiles(numbr,b,3)
!         end do
!       end if

      write(command, '("cp yesfolder/", i0, ".pun newlat/")') ibar
      call execute_command_line(command)

       write(131,*) 42
       write(131,*) vbarTMP, ibar, mia(numbr+1)    !!!! jun 30: mia used ibar3
       do b=1,42
         write(131,*) yesfiles(ibar,b,1), yesfiles(ibar,b,2), yesfiles(ibar,b,3)
       end do

     end if
   else if (vbarTMP>Ecut0) then
     call sort(qq2,d,nogo,leNo,riNo,msize,i,unique,iii)
     if (unique==1) then
       call addTo(ibar2,rino,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)

     else if (unique==-1) then
       call addTo(ibar2,leno,nogo,iii,qq2,d,pot2,vbartmp,msize,seampot2,se)
     end if
     if(unique.ne.0) then
       open(12, file=set1b, status="old", position="append", action="write")
       write(12,*) 2,ibar2
       write(12,*) qq2
       write(12,*)  vbarTMP
       close(12)
!       if(allfiles(numbr,1,1)<9998.0) then
!         do b=1,42
!           nofiles(ibar2-1,b,1) = allfiles(numbr,b,1)
!           nofiles(ibar2-1,b,2) = allfiles(numbr,b,2)
!           nofiles(ibar2-1,b,3) = allfiles(numbr,b,3)
!         end do
!       end if

       write(command, '("mv -f allfolder/", i0, ".pun nofolder/", i0, ".pun")') numbr, ibar2-1 
       call execute_command_line(command)

       write(113,*) vbarTMP, ibar2, mia(numbr+1)
     end if
   end if
 end do

 write(*,*) "allfiles distributed"
! deallocate(allfiles)

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
 deallocate(mia)
 deallocate(qq2dir)
 deallocate(may)
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
