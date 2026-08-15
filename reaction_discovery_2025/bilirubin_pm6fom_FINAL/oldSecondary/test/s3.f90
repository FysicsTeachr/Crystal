program main
implicit none
integer*1, dimension(14):: qq,qq2, mcPoint
real, dimension(14):: qt
integer, dimension (:), allocatable :: le,ri,rad,leNo,riNo,symmetric
real*8,dimension(:),allocatable:: pot, pot2
real::start,finish,time,dumdu,now,addIt,qShrink, acceptance
real*8:: Ecut0, vbar,ein, ein2, vbarTMP,se, emax
integer, dimension(:), allocatable ::shellArr

integer::j,k,l,ct,qc,qd,qe,qf,v,unique,zzz,irr,irr2,first,last,ci,cf,i,ibar,iii,d,iij,evaluations
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

 evaluations=0
 vbarTMP=55555.0

            open(unit=14, file="gradient.log")

            do
              read (14,"(a)",iostat=irr) text
              if (irr /= 0) then
                if (evaluations>10) then
                        energyChecked=1
                else        
                        energyChecked = 3
                end if
                exit
              end if

              if(index(text, "E (change) =").ne.0) then
                first= index(text, "E (change) =") +13
                xxx= text(first:first+12)
                xxx= trim(adjustl(xxx))
                evaluations=evaluations+1
              end if

              if(index(text, "Converged! =D").ne.0) then
                 energyChecked=1
                 exit
              end if
              if(index(text, "Converged").ne.0) then
                 energyChecked=1
                 exit
              end if
              if(index(text, "converged").ne.0) then
                 energyChecked=1
                 exit
              end if
              if(index(text, "CONVERGED").ne.0) then
                 energyChecked=1
                 exit
              end if
              if(index(text, "linalg").ne.0) then
                 energyChecked=3
                 exit
              end if
              if(index(text, "Linalg").ne.0) then
                 energyChecked=3
                 exit
              end if
              if(index(text, "linAlg").ne.0) then
                 energyChecked=3
                 exit
              end if
              if(index(text, "LinAlg").ne.0) then
                 energyChecked=3
                 exit
              end if
            end do

            close(14)

            emax=-272.365

            if (energyChecked==1) then
              read(xxx,*) ein     ! Convert xxx into real named ein
              vbarTMP=ein
            end if
            if (energyChecked==3) then
              open(134,file="../cancel.txt")
              write(134,*)"y"
              close(134)
              ein=99999.0
            end if

            open(14, file="../values.txt", status="unknown", position="append", action="write") 
            write(14,*) vbartmp
            close(14) 

            if (energyChecked==3 .or. vbartmp>emax) then
              open(134,file="proceed.txt")
              write(134,*)"n"
              close(134)
            end if
            if (energyChecked==1 .and. vbartmp<emax) then
              open(134,file="proceed.txt")
              write(134,*)"y"
              close(134)
            end if

end program
