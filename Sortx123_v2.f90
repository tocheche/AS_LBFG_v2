!***********************
!
Subroutine Sortx123(x,lbfg)
!Collects distances and indexes (ordering label(OL)) of principal atom (PA) and their first 4 neighbor atoms (NAs).
!The initial input config, x(3*(Nat1+Nat2+N4OUT)), is contained in file 202 (or 201).
!Sortx123() generates T1,T2,T3,T23 for a given x configuration; T1,T2,T3,T23 are used by Sorting4X to obtain files necessary to find the neighbors of: At_1(Ga) in temp240,temp241; At_2(Sb) in temp260,temp261; At_3(As) in temp280,temp281.
!Files temp240,temp260,temp280 are used in EnConfig() to obtain the data necessary to compute Uel(n,x) & DUel(LL,:).
!File 111 provides the indexes of each atom and its first 4 neighbors for the x configuration.
!
use DataType !for N13,N2IN,N4OUT,N24
integer, parameter :: dp = kind(1.0d0)
Real(dp):: T23c(8),T21c(8),T31c(8)
integer:: Nb1(8),Nb2(8),Nb3(8)
Real(dp):: T1(N13,3), T2(N2IN,3), T3(N4OUT,3), T23(N24,3)
integer:: i, i1, lbfg !,stgy
Real(dp):: x(3*(N13+N2IN+N4OUT))
!-----------------------------------------------
!Generating T1,T2,T3,T23,T23c,T21c,T31c,Nb1,Nb2,Nb3 from x:
do i=1,N13
    T1(i,1)=x(3*(i-1)+1);T1(i,2)=x(3*(i-1)+2);T1(i,3)=x(3*(i-1)+3)
    !write(451,*) T1(i,1),T1(i,2),T1(i,3)
    !write(900,*) T1(i,1),T1(i,2),T1(i,3)
end do
continue
!
do i=N13+1,N13+N2IN
    T2(i-N13,1)=x(3*(i-1)+1); T2(i-N13,2)=x(3*(i-1)+2); T2(i-N13,3)=x(3*(i-1)+3);
    T23(i-N13,1)=x(3*(i-1)+1);T23(i-N13,2)=x(3*(i-1)+2);
    T23(i-N13,3)=x(3*(i-1)+3);
    !write(452,*) T2(i-N13,1),T2(i-N13,2),T2(i-N13,3)
    !write(423,*) T23(i-N13,1),T23(i-N13,2),T23(i-N13,3)
   ! write(900,*) T2(i-N13,1),T2(i-N13,2),T2(i-N13,3)
end do
continue
!
do i=N13+N2IN+1,N13+N2IN+N4OUT
    T3(i-(N13+N2IN),1)=x(3*(i-1)+1);T3(i-(N13+N2IN),2)=x(3*(i-1)+2);
    T3(i-(N13+N2IN),3)=x(3*(i-1)+3); 
    T23(i-N13,1)=x(3*(i-1)+1);T23(i-N13,2)=x(3*(i-1)+2);
    T23(i-N13,3)=x(3*(i-1)+3); 
    !write(453,*) T3(i-N13,1),T3(i-N13,2),T3(i-N13,3)
    !write(423,*) T23(i-N13,1),T23(i-N13,2),T23(i-N13,3)
   ! write(900,*) T3(i-N13,1),T3(i-N13,2),T3(i-N13,3)
end do
continue
!-----------------------------------------------
!
! For Ga as a PA:
open(unit=240,file='temp240.dat',status='replace') !240 (241) is going to be used in EnConfig(..).
!open(unit=241,file='temp241.dat',status='replace')
Do i1=1,N13
Call Sorting4X(i1,T1,T23,T23c,N13,N24,Nb1) !There are 1,..j (j.ge.4) neighbors At_i2 (NA) of At_i1 (PA); the sequence 1,..j of At_i2 is rearranged in Sorting4X() and the first 4 numbers (indexes of the closest 4 NA) are kept.
write(240,*)T23c !T23c contains the first 4 increasing ordered PA-NA distances (given PA) & the first 4 corresponding NA indexes (4NAI) which are obtained by ordering the sequence NAI 1,2,..j (as function of PA-NA distance).
!write(241,*)Nb1  !4NAI and the 4 initial indexes (from T23) of the closest 4 NA. 
    if (lbfg==1)then
    write(111,*)i1,N13+Nb1(5),N13+Nb1(6),N13+Nb1(7),N13+Nb1(8) ! 111 provides the indexes of each PA (Ga atom) and its first 4 NA, all from the config.   
    end if   
continue
end Do
close(240) !;close(241) 
!-----------------------------------------------------------------
! For Sb as a PA:
open(unit=260,file='temp260.dat',status='replace') !260 (261) is going to be used in EnConfig(..).
!open(unit=261,file='temp261.dat',status='replace')
Do i1=1,N2IN
Call Sorting4X(i1,T2,T1,T21c,N2IN,N13,Nb2)
write(260,*)T21c !similar to the above T23c but NA are Ga atoms (that's why we put an index 1).
!write(261,*)Nb2  !4NAI and the 4 initial indexes (from T1) of the closest 4 NA (Nb2(:).le.Nat1). 
    if (lbfg==1)then
    write(111,*)i1+N13,Nb2(5),Nb2(6),Nb2(7),Nb2(8) !see above about 111, with PA Sb.
    end if
end Do
close(260) !;close(261) 
!-----------------------------------------------------------------
! For As as a PA:
open(unit=280,file='temp280.dat',status='replace') !280 (281) is going to be used in EnConfig(..).
!open(unit=281,file='temp281.dat',status='replace')
Do i1=1,N4OUT
Call Sorting4X(i1,T3,T1,T31c,N4OUT,N13,Nb3)
write(280,*)T31c !similar to the above T21c but PA are As atoms (that's why we put an index 3).
!write(281,*)Nb3  !!4NAI and the 4 initial indexes (from T1) of the closest 4 NA (Nb3(:).le.Nat1). 
    if (lbfg==1)then
    write(111,*)i1+N13+N2IN,Nb3(5),Nb3(6),Nb3(7),Nb3(8) !see  above about 111, with PA As.
    end if
end Do
close(280) !;close(281)
    end Subroutine Sortx123
!**************************************************************************************************