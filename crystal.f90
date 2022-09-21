!Updates in this version:
!1. Cleanup
!2. Can now declare symmetric coordinates to skip negative values in that dimension
!
!Readme notes:
!1. The gout.txt file initially should contain the energy at the equilibrium geometry,
!with whatever the convention specified within the crystal code.
!Currently, the gout.txt can be replaced by something like the following line*:
! SCF Done:  E(RTPSSh) =  -4444.4
!
!*Note that this file will change on execution,
!so it must be re-newed each time before submitting the code
!
!2. To run the file in Frontera after step 1, 
!enter the following 2 lines (from the same folder):
!ifort crystal.f90
!sbatch run.sh
!
!3. On running, the output file that will contain cordinates and energies will be "lattices.txt"
!
!
!Comments: 
!For now, the change in energy is used to determine if Gaussian calculations have completed. 
!Please make sure every symmetric cordinate q is set as symmetric(q)=1 inside the code. 


program main
implicit none
real, dimension (:,:),allocatable:: gog,nogo,arrae
real, dimension(:),allocatable:: qq,qq2, mcPoint
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real*8,dimension(:),allocatable:: pot, pot2
real*8::Ecut0,vbar,start,finish,time,step,dumdu,now,addIt, vbarTMP, ein, einOld
integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,first,last,ci,cf,i,ibar,iii,d
integer:: msize,p,s,chki,ibar2,z,dump,grth,leth,el,t,y,xx

character (len=1024) :: text
character (len=1024) :: xxx

real:: FeX, FeY, H1X,H1Y, H2X, H2Y, nscal

d=3  !!--------> This should be the dimensions of normal mode coordinates 

allocate(symmetric(d)) 
symmetric=0
symmetric(3) =1  !!--------> Set = 1 for dimensions that are symmetric. Others are 0.
                 !!Crystal skips the negative direction of these coordinates


Open(UNIT = 12, file="lattices.txt")      !Output lattices go here
write(12,*) 
Close(12)

dump=0

msize=5000000  ! max basis size (for allocating memory).
Ecut0= 10000       !!! cutoff energy for accepting lattices.
nscal=1.0        ! lattice spacing.

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

open(unit=180, file="radiux.txt") !        <<<<<<<<<<<<<<<
                                  
Allocate(rad(d))                  !                      |
                                  !                      |
!defined below for simplicial growth. otherwise, this array can be externally replaced by the desired replacement of growth type. In that case, comment out this part.

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
  if(symmetric(iii)==0) then
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


Allocate(arrae(ct,d))     !Now reading the above file
do i=1,ct
read(180,*) (arrae(i,j),j=1,d)
end do

allocate (gog(d,msize))              
allocate (nogo(d,msize))               
allocate (le(msize))
allocate (ri(msize))
allocate (leNo(msize))
allocate (riNo(msize))
allocate(qq(d))
allocate(qq2(d))
allocate(mcPoint(d))
allocate(pot(msize))
allocate(pot2(msize))

ci=1
cf=1
ibar=1
ibar2=1
le=-1
ri=-1
leNo=-1
riNo=-1
step=0


!			Defining the first seed point			
do i=1,d										
gog(i,1)=0	
nogo(i,1)=99999										
end do											


!do while (step<=0) !-->this loop is for modified-crystal termination algorithms.
ci=1
cf=ibar


open(unit=14, file="gout.txt")
do
   read (14,"(a)",iostat=irr) text ! Read line into character variable
   
   IF (irr/=0) EXIT
   
   if(index(text, "SCF Done").ne.0) then ! Check each line for substring
   
     first= index(text, "=");
     last= index(text,"A.U.");
     
     xxx= text(first+1:last-1);
     
   end if  
end do
read(xxx,*) ein                    ! Convert xxx into integer
einOld=ein
vbarTMP=ein
close(14)







do while (cf>=ci)
!!
if (ibar>msize/8) then !To avoid unwanted symmetry breaks, stop run and increase mzise next time
write(4,*) "killing code due to the msize/5 allocation", ibar
write(*,*) "killing code due to the msize/5 allocation", ibar
dump=1
exit
end if
!!

dumdu=0

do i=ci,cf	



  do chki=1,d
  qq(chki)=gog(chki,i)
  qq2(chki)=gog(chki,i)
  end do

  do t=2,ct
      do chki=1,d
        qq(chki)=gog(chki,i)+arrae(t,chki)*1.0
      end do

      do chki=1,d
       if (qq(chki) .eq. qq2(chki)) then
         mcPoint(chki)=qq(chki)
       else
         mcPoint(chki)=qq(chki) !!! Choose either qq(chki) or (qq(chki)+qq2(chki))/2.0
       end if
      end do


      call sort(mcPoint,d,nogo,leNo,riNo,msize,i,unique,iii)
      if(unique.ne.0) then
        dumdu=unique
        call sort(qq,d,gog,le,ri,msize,i,unique,iii)    
        if(unique.ne.0) then

          !!!!!!!!!!! Transformation Matrix goes here!!!!!!!!!!!!
          
          FeX=-1.55+0.01*qq(1)+0.02*qq(2)+0.015*qq(3)
          FeY=-0.05+0.02*qq(1)+0.01*qq(2)+0.01*qq(3)
          H1X=0.43+0.015*qq(1)
          H1Y=0.012*qq(1)+0.011*qq(2)
          H2X=-0.43-0.015*qq(1)+0.005*qq(3)
          H2Y= -0.011*qq(1)-0.012*qq(2)
          
          open(unit=24, file="gin.txt")
          write(24,*) "%NProcShared=56"
          write(24,*) "%mem=160GB"
          write(24,*) "#P TPSSh def2TZVPP Nosym opt=tight"
          write(24,*)
          write(24,*) "DFT_17 bill min"
          write(24,*)
          write(24,*) "+1, 1"
          write(24,*) "Fe" , " ", "-1", " ", FeX, " ", FeY, " ", 0.0
          write(24,*) "H" , " ", "-1", " ", H1X, " ", H1Y, " ", 0.0
          write(24,*) "H" , " ", "-1", " ", H2X, " ", H2Y, " ", 0.0
          WRITE(24,*) "H  0    -1.393088       -1.5650434      0.0"
          WRITE(24,*)"P  0    -3.6520628      -0.66765275     0.0"
          WRITE(24,*)"H  0    -4.1016344      -1.4668488      1.06641"
          WRITE(24,*)"H  0    -4.6950598      0.2823582       0.0"
          WRITE(24,*)"H  0    -4.1016344      -1.4668488      -1.06641"
          WRITE(24,*)"P  0    -1.4171272      -0.30480754     2.197"
          WRITE(24,*)"H  0    -1.6083279      0.76914505      3.09046"
          WRITE(24,*)"H  0    -2.2705079      -1.2382428      2.81076"
          WRITE(24,*)"H  0    -0.18798692     -0.77062316     2.69433"
          WRITE(24,*)"P  0    -1.9096042      2.1473263       0.0"
          WRITE(24,*)"H  0    -3.2399509      2.6126614       0.0"
          WRITE(24,*)"H  0    -1.4301754      2.9387032       -1.06171"
          WRITE(24,*)"H  0    -1.4301754      2.9387032       1.06171"
          WRITE(24,*)"P  0    -1.4171272      -0.30480754     -2.197"
          WRITE(24,*)"H  0    -0.18798692     -0.77062316     -2.69433"
          WRITE(24,*)"H  0    -2.2705079      -1.2382428      -2.81076"
          WRITE(24,*)"H  0    -1.6083279      0.76914505      -3.09046"
          write(24,*)
          
          close(24)


          do zzz=1,90000000

            open(unit=14, file="gout.txt")

            do
              read (14,"(a)",iostat=irr) text ! Read line into character variable
              IF (irr/=0) EXIT
              if(index(text, "SCF Done").ne.0) then ! Check each line for substring
   
                first= index(text, "=");
                last= index(text,"A.U.");
     
                xxx= text(first+1:last-1);

              end if  
            end do
            read(xxx,*) ein                   ! Convert xxx into integer
            close(14)


            
            if(ein.eq.9999999) then           !  Option to force exit
              write(*,*) "ibar is", " ", ibar
              write(*,*) "exiting. Modify here to save lattices to file"
              dump=1
              exit
            end if
            
            if(ein.eq.einOld) then
              call sleep(1)
            else
              vbarTMP=ein
              einOld=ein
              write(*,*) ein
              exit
            end if
          end do
          
          
          
          if(vbarTMP<=Ecut0 .and. step==0 .and. ein.ne.9999999) then
            if (unique==1) then
            !unique is 0 if the vector is already in the list
            !Else it is 1 or -1. Depends on if it goes to left or right branch at the end.
              call addTo(ibar,ri,gog,iii,qq,d,pot,vbartmp,msize)
              !adds vector to gog, and adds vector's label to ri(iii)
              
              open(12, file="lattices.txt", status="old", position="append", action="write")
              write(12,*) qq, vbartmp
              close(12)
              
            else if (unique==-1) then
              call addTo(ibar,le,gog,iii,qq,d,pot,vbartmp,msize)
              
              open(12, file="lattices.txt", status="old", position="append", action="write")
              write(12,*) qq, vbartmp
              close(12)
              
            end if
          else if(vbarTMP>Ecut0 .and. step==0)then
            if (dumdu==1) then
              call addTo2(ibar2,riNo,nogo,iii,mcPoint,d,pot2,vbartmp,msize)
            else if (dumdu==-1) then
              call addTo2(ibar2,leNo,nogo,iii,mcPoint,d,pot2,vbartmp,msize)
            end if
          end if
          
          
          
        end if
      end if
      
    if (dump==1) then   
      exit
    end if
  
  end do
  
  if (dump==1) then   
    exit
  end if
  
end do

ci=cf+1											!
cf=ibar	
call cpu_time(now)
!write(4,*) ibar, now-start
!write(4,*) "ibar,ci,cf: ", ibar, ci, cf

!!if (step>=1) then  !steps where you want to terminate after one iteration
!!exit
!!end if

end do   										!
!_______________________________________________!

!write(4,*) ecut0, "   ", ibar ,"   ",now-start

!do i=1,ibar
!  write(4,*) (gog(chki,i),chki=1,d)
!end do
!do i=1,ibar2
!  write(4,*) (nogo(chki,i),chki=1,d)
!end do


call cpu_time(finish)
time=finish-start
!write(*,*) time

deallocate (gog)              
deallocate (nogo)               
deallocate (le)
deallocate (ri)
deallocate (leNo)
deallocate (riNo)
deallocate(qq)
deallocate(qq2)
deallocate(mcPoint)
deallocate(pot)
deallocate(pot2)
deallocate(rad)                    !                   |
deallocate(arrae)

rewind(180)

!write(*,*) "ibar, ibar2" , ibar, ibar2

!end do
end program


function grth(d,qq,i,iii,msize,gog)			!returns 1 if greater, else 0
integer:: d,i,iii,a,grth,z,xx,b,msize
real, dimension(d,msize):: gog
real, dimension(d):: qq
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


function leth(d,qq,i,iii,msize,gog)			!returns 1 if lesser, else 0
integer:: d,i,iii,a,leth,z,xx,b,msize
real, dimension(d,msize):: gog
real, dimension(d):: qq
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
real,dimension(d)::q,qq
real,dimension(d,msize)::gogo
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


subroutine addTo(iba,ri,gogo,iii,qq,d,pot,vbartmp,msize)
integer:: iba,iii,el,d
integer, dimension(msize):: ri
real,dimension(d,msize)::gogo
real,dimension(d)::qq
real*8,dimension(msize)::pot
iba=iba+1	
ri(iii)=iba						!
do el=1,d
  gogo(el,iba)=qq(el)
end do
pot(iba)=vbartmp
end subroutine

subroutine addTo2(iba,ri,gogo,iii,qq,d,pot2,dimen,msize)
integer:: iba,iii,el,d
integer, dimension(msize):: ri
real,dimension(d,msize)::gogo
real,dimension(d)::qq
real*8,dimension(msize)::pot2
iba=iba+1	
ri(iii)=iba						!
do el=1,d
  gogo(el,iba)=qq(el)
end do
pot2(iba)=dimen
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
