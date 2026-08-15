program main
      implicit none
      real*8:: x,y,z
      integer:: irr
      character:: atom

      open(unit=11, file="pregeom.xyz")
      open(unit=12, file="gradient.dat")
read(11,*)
read(11,*)
write(12,'(A)') "PM7  1SCF CI=(4,2) OPEN=(4,4) GRADIENTS GEO-OK ROOT=1 FROCC FLOCC=0.10 PREC VECTORS MXROOT=5  MECI=5  +"
!write(12,'(A)') "AM1  1SCF GRADIENTS GEO-OK PREC VECTORS +"
write(12,'(A)') "SINGLET S2 MEMORY=2000"
write(12,'(A)')
write(12,'(A)') "geometry"

Do
read(11,*,iostat=irr) atom, x, y, z
if(irr/=0) then
        exit
end if
!write(12,*'(A,F15.10," 1")') atom, x, "1", y, "1" ,z,"1"
write(12,'(A,1X,F14.10," 1 ",F14.10," 1 ",F14.10," 1")') atom,x,y,z
end do

write(12,*)
write(12,'(A)') "OCCUP"
write(12,'(A)') "2.0 2.0 0.0 0.0"

end program
