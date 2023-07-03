   !TO BE DELETED/REPLACED:
!v07
!(visualize with ZincBlende-structure-planes.nb)    
!=========================================
    
Subroutine FCCbox(nx,ny,nzmin,nzmax)
    ! Generates FCC atom positions in a box with latt_ct=1; for each set (i,j,k)v07 8 atom locations are generated; the box has 8*(2*nx+1)*(2*ny+1)*(nzmin+nzmax+1) total # atoms! In x and y directions the coords vary between -nx,nx and -ny,ny, respectively; in z direction the coords vary between -nzmin,nzmax.
    ! Atoms are placed in 4 types of xy planes that we name layer types (LTs) as follows. Let LT1 be the plane z=0. Then for  k>0:
    !       LT3 is upwardly shifted by k*a/2 with respect to LT1; LT1 and LT3 are filled with At1.
    !       LT2 is upwardly shifted by k*3*a/4 with respect to LT1.
    !       LT4 is upwardly shifted by k*a/4 with respect to LT1; LT2 and LT4 are filled with At2.
    !       LT2, LT3, LT4 are horizontally shifthed with respect to LT1 (visualize with ZincBlende-structure-planes.nb).
    ! Any other orientation of the BULK atoms can EASILY be obtained with this subroutine by rotations.   
    ! 11 collects: type At1(Ga), index(ordering label), FCC bulk atom coords; z layers at 0,-0.5,0.5,...
    ! 12 collects: type At2(Sb or As), index(ordering label), FCC bulk atom coords; z layers at 0,-0.25,0.25,-0.75, 0.75,......
Integer:: nx,ny,nzmin,nzmax,i,j,k, n, m ! n is index   
n=0; m=0;
    Do 1 i=-nx,nx
        Do 1 j=-ny,ny
                Do 1 k=-nzmin,nzmax
            n=n+1
             m=n*4-3        
    Write(11,*) 1, m, i, j, k                        !At1 in LT1
    Write(11,*) 1, m+1, (2*i+1)/2d0, (2*j-1)/2d0, k  !At1 in LT1
    Write(11,*) 1, m+2, i, (2*j+1)/2d0, (2*k+1)/2d0  !At1 in LT3
    Write(11,*) 1, m+3, (2*i+1)/2d0, j, (2*k+1)/2d0  !At1 in LT3
   ! 
    Write(21,*) 2, m, (4*i-1)/4d0, (4*j+1)/4d0, (4*k+3)/4d0   !At2 in LT2
    Write(21,*) 2, m+1, (4*i+1)/4d0, (4*j-1)/4d0, (4*k+3)/4d0 !At2 in LT2
    Write(21,*) 2, m+2, (4*i+1)/4d0, (4*j+1)/4d0, (4*k+1)/4d0 !At2 in LT4
    Write(21,*) 2, m+3, (4*i+3)/4d0, (4*j-1)/4d0, (4*k+1)/4d0 !At2 in LT4
    
1   continue 
End Subroutine FCCbox
!*******************************************************************
!  
Subroutine QDHTF(x,y,z,log)
    !Generates hemi-torus (HT) function in Cartesian coords.
    !Rq-torus radius; Rt-radius of HT; Rt1-internal radius of HT; Rt2-outer radius of HT.
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp):: Rt1, Rt2, x, y, z, rho
Logical:: log !, log1, log2
    Rt1=Rt-Rq
    Rt2=Rt+Rq  
	rho=dsqrt(x**2+y**2)        
    log=((rho.ge.Rt1).and.(rho.le.Rt2).and.((rho-Rt)**2+z**2.le.Rq**2).and.(z.gt.0d0))! log=true for (x,y,z) inside+boundary HT    
!continue        
End Subroutine QDHTF
!*********************************************     
!    
Subroutine QDConeF(x,y,z,log)
    !Generates Cone function in Cartesian
    !Rc-cone radius; h-cone height
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp)::  x, y, z, rho
Logical:: log
    rho=dsqrt(x**2+y**2)   
    log=(rho.le.(h-z)*Rc/h) .and.(z.le.h).and.(z.gt.0d0)! log=true for (x,y,z) inside the Cone and its surface
!continue
    End Subroutine QDConeF
!*********************************************     
!    
Subroutine QDConeTruncF(x,y,z,log)
    !Generates truncated cone function in Cartesian
    !Rct-cone radius; height=ht (height of cone); height truncated cone=hct (distance between the centers of the two bases) 
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp)::  x, y, z, rho
Logical:: log
    rho=dsqrt(x**2+y**2)   
    log=(rho.le.(ht-z)*Rct/ht) .and.(z.le.hct).and.(z.gt.0d0)! log=true for (x,y,z) inside the Cone and its surface
    if (hct.ge.ht) then
    write (6,61)
61  format ('The height of truncated cone hct can not be larger than its height. Change hct in VFF_v2.inp.')
    STOP
    end if
!continue
    End Subroutine QDConeTruncF
!*********************************************
!    
Subroutine QDPyramidF(x,y,z,log)
    !Generates rectangular pyramid function in Cartesian
    !b=Rc-half of the baze; h-pyramid height
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp)::  x, y, z, fpyr
Logical:: log
    fpyr=(h-z)*Rc/h   
   log=(x.ge.-fpyr).and.(x.le.fpyr).and.(y.ge.-fpyr).and.(y.le.fpyr).and.(z.gt.0d0).and.(z.le.h) ! log=true for (x,y,z) inside the pyramid and its surface
!continue
    End Subroutine QDPyramidF
!*********************************************
!    
Subroutine CSQDF(x,y,z,log)
    !Generates core-shell QD
    !RCO-core radius
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp)::  x, y, z
Logical:: log
log=(x**2+y**2+z**2).le.RCO**2 !(x,y,z) inside the core and its surface
!continue
End Subroutine CSQDF
!*********************************************      
!
Subroutine SelectQDT(x1,y1,z1,x2,y2,z2,log1,log2,QD)
! One selects the QD type; QD=1 for hemi-torus (HT); QD=2 for cone; QD=21 for truncated cone; QD=3 for pyramid; QD=4 for core-shell
integer, parameter :: dp = kind(1.0d0)
integer:: QD
Logical:: log1 , log2
Real(dp):: x1,y1,z1 ,x2,y2,z2
!
If (QD==1) then ! for HT
    Call QDHTF(x1,y1,z1,log1)  !interogates if x1,y1,z1 are IN or OUT HT
    Call QDHTF(x2,y2,z2,log2)  !interogates if x2,y2,z2 are IN or OUT HT
elseif (QD==2) then ! for Cone
    Call QDConeF(x1,y1,z1,log1) !interogates if x1,y1,z1 are IN or OUT Cone  
    Call QDConeF(x2,y2,z2,log2) !interogates if x2,y2,z2 are IN or OUT Cone
elseif (QD==21) then ! for Cone
    Call QDConeTruncF(x1,y1,z1,log1) !interogates if x1,y1,z1 are IN or OUT Cone  
    Call QDConeTruncF(x2,y2,z2,log2) !interogates if x2,y2,z2 are IN or OUT Cone
elseif (QD==3) then ! for Pyramid
    Call QDPyramidF(x1,y1,z1,log1) !interogates if x1,y1,z1 are IN or OUT Pyramid  
    Call QDPyramidF(x2,y2,z2,log2) !interogates if x2,y2,z2 are IN or OUT Pyramid
elseif (QD==4) then ! for core-shell
    Call CSQDF(x1,y1,z1,log1) !interogates if x1,y1,z1 are in core 
    Call CSQDF(x2,y2,z2,log2) !interogates if x2,y2,z2 are in core !*shell
End if
continue
    end  Subroutine SelectQDT  
!**************************************************************************************************  
!  
Subroutine  QDWLB6(QD,ss,x1,y1,z1,x2,y2,z2, log1,log2,i1IN,i1OUT,i1,i2,i3,WL,RCS)!generates FCC coordinates for IN and OUT QDWL.
!Writes file 140,160,180 which contain coords of Ga, Sb, As, respectively.
!Writes file 38,39 which contain coords of Ga IN, OUT QDWL, respectively.
!
    use DataType
    integer, parameter :: dp = kind(1.0d0)
    Logical:: log1, log2
    Integer:: ss,i1IN,i1OUT,i1,i2,i3,QD
    Real(dp):: x1,y1,z1,x2,y2,z2 ,WL, RCS
!------------------------------------------------------------
!=============================================
!@IF(QD==1 .or. QD==2 .or. QD==3) then 
IF(QD==1 .or. QD==2 .or. QD==21 .or. QD==3) then 
    If(log2==.True.) then ! coords (x2=a2*At(ss,n,3),...) for Ga & Sb, IN QD 
        !"if ((dsqrt(x2**2+y2**2).le.(nx+latshift2)*a1)) then
        if (ss==1)then !for Ga IN (for lat_ct=1 one has At(1,n,5)=z=0,0.5,1,..)
            i1=i1+1 !here, i1 counts #Ga IN QD (in Main, i1 also counts Ga atoms in WL); i1 counts total #Ga atoms
            i1IN=i1IN+1 !here, i1IN counts #Ga IN QD (in Main, i1IN counts Ga atoms in WL)
            write(38,*) x2,y2,z2 !in Main, in file 38 are also written Ga coords in WL
            write(140,*) x2,y2,z2 !in Main, in file 140 are also written Ga coords in WL
        end if
!            
        if (ss==2)then !for Sb IN QD (for lat_ct=1 one has here At(2,n,5)=z=0.25,0.75,..)
            i2=i2+1 !for Sb IN QDWL (in Main, i2 also counts Sb atoms in WL)
            write(160,*) x2,y2,z2 !in Main, in file 160 are also written Sb coords in WL
        end if
        !"End if
    End If        
!----------------  
!        
IF(log2==.False.) then ! coords (x2=a2*At(ss,n,3),...) are OUT of QD
    iF (QD==3.and.(dabs(x2).le.(nx-latshift1)*a1).and. (dabs(y2).le.(ny-latshift1)*a1)) then !for pyramid QD 
        if (ss==1 .and. (z2.gt.0d0))then !Ga capping atom (for lat_ct=1 one has here At(1,n,5)=z=0,0.5,1,..)
            i1=i1+1 !i1 counts #Ga OUT QD            
            write(39,*)  x2,y2,z2 !Ga capping atom
            i1OUT=i1OUT+1 !counts #Ga OUT QDWL             
            write(140,*) x2,y2,z2 !counts total #Ga atoms
        end if
!
        if (ss==2 .and. (z2.gt.0d0))then !for As capping atom (for lat_ct=1 one has At(2,n,5)=z=0.25,0.75,..)
            i3=i3+1
            write(180,*) x2,y2,z2 !As capping atom
        end if
    end iF
!----
!@    iF ((QD==1.or.QD==2).and.(dsqrt(x2**2+y2**2).le.(nx-latshift2)*a1)) then ! for cone & hemi-torus QD
iF ((QD==1.or.QD==2.or.QD==21).and.(dsqrt(x2**2+y2**2).le.(nx-latshift2)*a1)) then ! for cone & hemi-torus QD
        if (ss==1 .and. (z2.gt.0d0))then !Ga capping atom (for lat_ct=1 one has here At(1,n,5)=z=0,0.5,1,..)
            i1=i1+1 !i1 counts #Ga OUT QD            
            write(39,*)  x2,y2,z2 !Ga capping atom
            i1OUT=i1OUT+1 !counts #Ga OUT QDWL             
            write(140,*) x2,y2,z2 !counts total #Ga atoms
        end if
!
        if (ss==2 .and. (z2.gt.0d0))then !for As capping atom (for lat_ct=1 one has At(2,n,5)=z=0.25,0.75,..)
            i3=i3+1
            write(180,*) x2,y2,z2 !As capping atom
        end if
    end iF
!---- 
End IF 
!------------------------------------------------------------
!
IF(log1==.False.) then ! coords (x1=a1*At(ss,n,3),...) are OUT of QD
    iF (QD==3.and.(dabs(x1).le.(nx-latshift1)*a1).and. (dabs(y1).le.(ny-latshift1)*a1)) then !for pyramid QD 
     !----              
        If((z1 .lt. -WL)) then  ! coords (x1=a1*At(ss,n,3),...) are IN SUBSTRATE
            if (ss==1) then !for Ga in SUBSTRATE (for lat_ct=1 here one selects At(1,n,5)=z=-1,-1.5,..)
            i1=i1+1 !here, for Ga in SUBSTRATE (in Main i1 also counts  Ga in WL); i1 counts total # Ga atoms
            i1OUT=i1OUT+1 !counts #Ga OUT QDWL         
            write(39,*) x1,y1,z1 !Ga in SUBSTRATE; in Main, in file 140 are also written Ga coords in WL 
            write(140,*) x1,y1,z1 !Ga in SUBSTRATE; in Main, in file 140 are also written Ga coords in WL
            end if
!           
            if (ss==2)then !for As IN SUBSTRATE (for lat_ct=1 one has At(2,n,5)=z=-0.75, -1.25..)
            i3=i3+1
            write(180,*) x1,y1,z1 !for As in SUBSTRATE.
            end if        
        end If
    !----  
    end iF
!
!@    iF ((QD==1.or.QD==2).and.(dsqrt(x1**2+y1**2).le.(nx-latshift1)*a1)) then ! for cone & hemi-torus QD
    iF ((QD==1.or.QD==2.or.QD==21).and.(dsqrt(x1**2+y1**2).le.(nx-latshift1)*a1)) then ! for cone & hemi-torus QD
    !----              
        If((z1 .lt. -WL)) then  ! coords (x1=a1*At(ss,n,3),...) are IN SUBSTRATE
            if (ss==1) then !for Ga in SUBSTRATE (for lat_ct=1 here one selects At(1,n,5)=z=-1,-1.5,..)
            i1=i1+1 !here, for Ga in SUBSTRATE (in Main i1 also counts  Ga in WL); i1 counts total # Ga atoms
            i1OUT=i1OUT+1 !counts #Ga OUT QDWL         
            write(39,*) x1,y1,z1 !Ga in SUBSTRATE; in Main, in file 140 are also written Ga coords in WL 
            write(140,*) x1,y1,z1 !Ga in SUBSTRATE; in Main, in file 140 are also written Ga coords in WL
            end if
!           
            if (ss==2)then !for As IN SUBSTRATE (for lat_ct=1 one has At(2,n,5)=z=-0.75, -1.25..)
            i3=i3+1
            write(180,*) x1,y1,z1 !for As in SUBSTRATE.
            end if        
        end If
    !----  
    end iF
End IF
END IF !for QD==1 .or. QD==2 .or.QD==21 .or. QD==3
!=============================================
!
IF(QD==4)then !for core-shell
!----------------    
    If(log2==.True.) then ! coords (x2=a2*At(ss,n,3),...) for Ga & Sb, IN QD 
        if (ss==1)then !for AT1(Ga) IN (for lat_ct=1 one has At(1,n,5)=z=..-0.5,0,0.5,1,..).
            i1=i1+1 !i1 counts total (#AT1+#At3)(Ga); i1_max=N13.
            i1IN=i1IN+1 !counts #AT1(Ga) IN core; i1IN_max=N1IN.
            write(38,*) x2,y2,z2 !AT1(Ga) IN core.
            write(140,*) x2,y2,z2 !AT1(Ga) IN core.
        end if
!            
        if (ss==2)then !AT2(Sb) IN core (for lat_ct=1 one has At(2,n,5)=z=..,-0.25,0.25,0.75,..).
            i2=i2+1 !#AT2(Sb); i2_max=Nat2.
            write(160,*) x2,y2,z2 !AT2(Sb) IN core .
        end if
    End If        
!----------------     
IF(log1==.False.) then ! coords (x1=a1*At(ss,n,3),...) are OUT of core.
    If((x1**2+y1**2+z1**2) .le. RCS**2) then !ensures spherical shape of core-shell.
            if (ss==1) then !for AT3(Ga) OUT (for lat_ct=1 here one selects At(1,n,5)=z=..,-1,-1.5,..)
            i1=i1+1 !i1 counts total (#AT1+#At3)(Ga); i1_max=N13.
            i1OUT=i1OUT+1 !counts #AT3(Ga) OUT of core; i1OUT_max=N3OUT.      
            write(39,*) x1,y1,z1 !AT3(Ga) OUT of core. 
            write(140,*) x1,y1,z1 !i1 counts total (#AT1+#At3)(Ga); i1_max=N13.
            end if
!           
            if (ss==2)then !for AT4(As) OUT QD (for lat_ct=1 one has At(2,n,5)=z= -1.25, -0.75..)
            i3=i3+1
            write(180,*) x1,y1,z1 !for AT4(As) OUT QD
            end if
   End If
End IF
!----------------
END IF !for QD==4,for core-shell
!   
End Subroutine QDWLB6
!******************************************************
!   
Subroutine IniConfig6(x) !T1,T2,T3 rearrange files 140,160,180 & generates T23 (140+160); writes 201,202 containing initial x.
!x from file 201,202 has the form: x(1),x(2),x(3),..,x(3*N13), x(3*N13+1),..x(3*(N13+Nat2)), x(3*(N13+Nat2)+1),..x(3*(N13+Nat2+Nat3)), which in Cartesian would be (x1,y1,z1), (x2,y2,z2),...(xN,yN,zN).
!
use DataType
integer, parameter :: dp = kind(1.0d0)
Real(dp)::T1(N13,3),T10(N13,3),T2(N2IN,3),T20(N2IN,3),T160(3*N2IN)
Real(dp)::T3(N4OUT,3),T30(N4OUT,3),T180(3*N4OUT) !,T23((Nat2+Nat3),3)
Real(dp)::T1IN(N1IN,3), T3OUT(N3OUT,3), T2IN(N2IN,3), T4OUT(N4OUT,3)
Real(dp)::T1IN0(N1IN,3), T3OUT0(N3OUT,3), T2IN0(N2IN,3), T4OUT0(N4OUT,3)
Real(dp)::x(3*(N1IN+N3OUT+N2IN+N4OUT)) 
Integer:: N132, N1324 !N13,
!##########################################
!
N13=N1IN+N3OUT;
N132=N1IN+N3OUT+N2IN;
N1324=N1IN+N3OUT+N2IN+N4OUT;
!
rewind(38) !collects At1 IN(Ga)
read(38,*) T1IN0
T1IN=RESHAPE(T1IN0, (/N1IN,3/), ORDER = (/2,1/)) !Rearrange in C style
do i=1,3*N1IN,3
    x(i)=T1IN(i/3+1,1);x(i+1)=T1IN(i/3+1,2);x(i+2)=T1IN(i/3+1,3)
    write(201,*) x(i),x(i+1),x(i+2)
    write(202,*) x(i),x(i+1),x(i+2)
end do
continue
!
rewind(39) !collects At3 OUT(Ga)
read(39,*) T3OUT0
T3OUT=RESHAPE(T3OUT0, (/N3OUT,3/), ORDER = (/2,1/)) !Rearrange in C style
do i=3*N1IN+1, 3*N13,3
    x(i)=T3OUT(i/3-N1IN+1,1); x(i+1)=T3OUT(i/3-N1IN+1,2); x(i+2)=T3OUT(i/3-N1IN+1,3)
    write(201,*) x(i),x(i+1),x(i+2)
    write(202,*) x(i),x(i+1),x(i+2)
end do
continue
!
rewind(160) !collects At2 IN(Sb)
read(160,*) T2IN0
T2IN=RESHAPE(T2IN0, (/N2IN,3/), ORDER = (/2,1/)) !Rearrange in C style
do i=3*N13+1, 3*N132, 3
    x(i)=T2IN(i/3-N13+1,1); x(i+1)=T2IN(i/3-N13+1,2); x(i+2)=T2IN(i/3-N13+1,3)
    write(201,*) x(i),x(i+1),x(i+2)
    write(202,*) x(i),x(i+1),x(i+2)
end do
continue

!collects At4 OUT(As)
rewind(180)
read(180,*) T4OUT0
T4OUT=RESHAPE(T4OUT0, (/N4OUT,3/), ORDER = (/2,1/)) !Rearrange in C style
do i=3*N132+1, 3*N1324, 3
    x(i)=T4OUT(i/3-N132+1,1); x(i+1)=T4OUT(i/3-N132+1,2); x(i+2)=T4OUT(i/3-N132+1,3)
    write(201,*) x(i),x(i+1),x(i+2)
    write(202,*) x(i),x(i+1),x(i+2)
end do
continue
!##########################################
!collects At1(Ga)
!rewind(140)
!read(140,*) T10
!T1=RESHAPE(T10, (/N13,3/), ORDER = (/2,1/)) !Rearranges T10 in C style
!do i=1,3*N13,3
    !x(i)=T1(i/3+1,1);x(i+1)=T1(i/3+1,2);x(i+2)=T1(i/3+1,3)
    !write(201,*) x(i),x(i+1),x(i+2)
    !write(202,*) x(i),x(i+1),x(i+2)
!end do
!continue
!
    end Subroutine IniConfig6
!*******************************************************************
!    
subroutine elastic_ct
use DataType
integer, parameter :: dp = kind(1.0d0) ! SCHIMBA AICI
!
!bond-stretching force constants used in U, DU calculus
!
    ac13=(C11out+3*C12out)*a1/4 !alpha1 GaAs OUT; C OUT
    bc13=(C11out-C12out)*a1/4 !beta1 GaAs;C OUT; C OUT
    ac12=(C11in+3*C12in)*a2/4 !alpha2 GaSb/InAs;Si IN
    bc12=(C11in-C12in)*a2/4   !beta2 GaSb/InAs;Si IN
!
    dc12=a2*dsqrt(3d0)/4 !GaSb/InAs;Si IN
    dc13=a1*dsqrt(3d0)/4 !GaAs;C OUT
!
!a1-GaAs;C (OUT) < a2-GaSb/InAs;Si(IN)
!ac12=ac21,ac13=ac31,bc12=bc21,bc13=bc31,dc12=dc21,dc13=dc31
end Subroutine elastic_ct
!******************************************************************
!    
Subroutine Sorting4X(i1,T1,T2,T12c,Natt1,Natt2,Nb)
!Sorting4X finds the first 4 neighbors of At_i1 belonging to the sequence T1.
!T2-sequence (indexes) of all neighbors of the At_i1. At_i1 can be of type At_1(ga), At_2(Sb), or At_3(As). T2 is formed by the NA indexes corresponding to At_i1 as PA.
!T12c stores ascedingly the first 4 ordered distances (d12) between At_i1 (PA) and its NAs from T2.
!Natt1,Natt2 are dimensions of T1 and T2.
!ind12s-initial index for neighbors which locate inside the sphere of radius = radius.
!ind12p-new index for neighbors which are located inside a sphere of radius = radius in the sequence T2 to be ordered.
!Nb-first 4 neighbor indexes ordered as function of ascendingly  associated  distance d12 & the same 4 atom neighbors with their index in the initial T2 sequence.
!Sorting4X dynamically adapts the magnitude of radius for an efficient ordering; an initial value should be chosen (try between 2 and 5 for radius in the inp file, by default radius=4) for minimum time of calculus.
!    
    use DataType
    implicit none
integer, parameter :: dp = kind(1.0d0)
integer::i1, i2, k, j
real(dp) :: d12, eps, radiusV
integer, intent(in) :: Natt1, Natt2
real(dp) :: T1(Natt1,3), T2(Natt2,3)
! real(dp) :: T2IN(Natt2,3)
real(dp), allocatable :: T12p(:)
integer,allocatable :: ind12p(:),ind12s(:)
real(dp), Dimension (8) :: T12c
integer,Dimension (8) :: Nb
eps=1.0_dp-10 !parameter in SR hpsort_eps_epw(..)
j=0; radiusV=radius
!Default radius=4d0 is the initial radius of a sphere centered on PA At_i1 for searching its NAs. For example, the At_1 (Ga) neighbors are At_2(Sb) and/or At_3(As).
open(unit=224,file='tempx.dat',status='replace')
! https://stackoverflow.com/questions/38176611/overwrite-a-file-using-fortran
11  continue
    do i2=1, Natt2
!*         d12=dsqrt((T1(i1,1)-T2(i2,1))**2+(T1(i1,2)-T2(i2,2))**2+(T1(i1,3)-T2(i2,3))**2) !distance At_i1-At_i2 (if i2<=Natt2) or At_i1-At_i3 (if i2>Natt2)      
       d12=dsqrt((T1(i1,1)-T2(i2,1))**2+(T1(i1,2)-T2(i2,2))**2+(T1(i1,3)-T2(i2,3))**2) !distance At_i1-At_i2 (if i2<=Natt2) or At_i1-At_i3 (if i2>Natt2)      
        if (d12.lt.radiusV) then
          j=j+1 !index counting atoms inside the sphere (of radiusV)
          write(224,*)i2,j,d12 
!224 is a temporary file for storing the data: i2-initial index of NA in T2, j-new ordering index (1,2,3...) given to the NA indexed by i2 in T2; d12 is the distance between At_i1 and its neighbor At_i2 (from T2).
        end if 
    end do
continue
if (j.lt.4) then ! dynamically adaptating the radiusV magnitude.
    radiusV=radiusV+1 !increases the radius of the sphere to find at least 4 NAs of a given PA.
    j=0
    close(224) !delete tempx.dat if contains data for less than 4 NAs.
    goto 11 !loop back to find at least 4 NAs.
end if
!continue
!
allocate (T12p(j),ind12p(j),ind12s(j)) 
!T12p for d12, ind12s for the initial index of NAs, ind12p is the new ordering index (1,2,3...j) of NAs to be sorted. 
rewind(224)
do k=1,j
read(224,*)ind12s(k),ind12p(k),T12p(k)
end do
continue
    call hpsort_eps_epw (j, T12p, ind12p, eps) !ascendengly sorting distances At_1-At_2 or At_1-At_3; j.le.4.
    !ind12p on the above row is the reordered ind12p(k) from the above command 'read(224,*)ind12s(k),ind12p(k),T12p(k)'
!continue 
    T12c(1)=T12p(1); T12c(2)=T12p(2);T12c(3)=T12p(3);T12c(4)=T12p(4); !the first 4 ordered distances d12, PA(i1)-NA1,_4.
    T12c(5)=ind12p(1);T12c(6)=ind12p(2);T12c(7)=ind12p(3);T12c(8)=ind12p(4); 
    !On the above and below rows: PA(i1) has NAs with the indexes 1,2,,,j (j.ge.4). After ascendengly sorting distances PA(i1)-NA1,_j it results a new ordering of NA indexes; one keeps the first ordered 4 NA indexes, ind12p(1),...ind12p(4).
   Nb(1)=ind12p(1);Nb(2)=ind12p(2);Nb(3)=ind12p(3);Nb(4)=ind12p(4)          ! --- II ---
   Nb(5)=ind12s(Nb(1));Nb(6)=ind12s(Nb(2));Nb(7)=ind12s(Nb(3));Nb(8)=ind12s(Nb(4)) !turns the initial index of neighbor atoms from T2
        !Above: For i=5,...8, IF Nb(i)>Natt2 then NA is At_2(Sb), ELSE NA is At_3 (As).
continue
close(224) !tempx.dat is deleted
continue 
    end Subroutine Sorting4X
!******************************************************************************************************
!
    Subroutine hpsort_eps_epw (n, ra, ind, eps)
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .                                                                                         
  ! Adapted from flib/hpsort_eps
  !---------------------------------------------------------------------
  !I insignificantly changed the subroutine 
  !---------------------------------------------------------------------
  ! sort an array ra(1:n) into ascending order using heapsort algorithm,
  ! and considering two elements being equal if their values differ
  ! for less than "eps".
  ! n is input, ra is replaced on output by its sorted rearrangement.
  ! create an index table (ind) by making an exchange in the index array
  ! whenever an exchange is made on the sorted data array (ra).
  ! in case of equal values in the data array (ra) the values in the
  ! index array (ind) are used to order the entries.
  ! if on input ind(1)  = 0 then indices are initialized in the routine,
  ! if on input ind(1) != 0 then indices are assumed to have been
  !                initialized before entering the routine and these
  !                indices are carried around during the sorting process
  !
  ! no work space needed !
  ! free us from machine-dependent sorting-routines !
  !
  ! adapted from Numerical Recipes pg. 329 (new edition)
  !
  implicit none  
  !-input/output variables
integer, parameter    :: dp = kind(1.0d0)
integer, intent(in)   :: n  
real(dp), intent(in)  :: eps
integer               :: ind (n)  
real(dp)              :: ra(n)
!
  !-local variables
integer  :: i, ir, j, l, iind  
real(dp) :: rra
!
  ! initialize index array
  IF (ind (1) .eq.0) then  
     DO i = 1, n  
        ind (i) = i  
     ENDDO
  ENDIF
  ! nothing to order
  IF (n.lt.2) return  
  ! initialize indices for hiring and retirement-promotion phase
  l = n / 2 + 1  

  ir = n  

  sorting: do 
  
    ! still in hiring phase
    IF ( l .gt. 1 ) then  
       l    = l - 1  
       rra  = ra (l)  
       iind = ind (l)  
       ! in retirement-promotion phase.
    ELSE  
       ! clear a space at the end of the array
       rra  = ra (ir)  
       !
       iind = ind (ir)  
       ! retire the top of the heap into it
       ra (ir) = ra (1)  
       !
       ind (ir) = ind (1)  
       ! decrease the size of the corporation
       ir = ir - 1  
       ! done with the last promotion
       IF ( ir .eq. 1 ) then  
          ! the least competent worker at all !
          ra (1)  = rra  
          !
          ind (1) = iind  
          exit sorting  
       ENDIF
    ENDIF
    ! wheter in hiring or promotion phase, we
    i = l  
    ! set up to place rra in its proper level
    j = l + l  
    !
    DO while ( j .le. ir )  
       IF ( j .lt. ir ) then  
          ! compare to better underling
          IF ( hslt( ra (j),  ra (j + 1) ) ) then  
             j = j + 1  
          !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
             ! this means ra(j) == ra(j+1) within tolerance
           !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
          ENDIF
       ENDIF
       ! demote rra
       IF ( hslt( rra, ra (j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
       !else if ( .not. hslt ( ra(j) , rra ) ) then
          !this means rra == ra(j) within tolerance
          ! demote rra
         ! if (iind.lt.ind (j) ) then
         !    ra (i) = ra (j)
         !    ind (i) = ind (j)
         !    i = j
         !    j = j + j
         ! else
             ! set j to terminate do-while loop
         !    j = ir + 1
         ! endif
          ! this is the right place for rra
       ELSE
          ! set j to terminate do-while loop
          j = ir + 1  
       ENDIF
    ENDDO
    ra (i) = rra  
    ind (i) = iind  
continue
  END DO sorting 
!write(41,*) ra, ind
contains 

  !  internal function 
  !  compare two real number and return the result

  logical function hslt( a, b )
  implicit none
    REAL(DP) :: a, b
    IF( abs(a-b) <  eps ) then
      hslt = .false.
    ELSE
      hslt = ( a < b )
    end if
  end function hslt  !
    end subroutine hpsort_eps_epw
!**************************************************************!************************************
!
!Two slightly different constraints imposed to the WL-matrix interface atoms may be chosen by setting the bopt parameter to optimize the evaluations in the line search; bopt=1 or 2 to select between the 2 settings.
Subroutine opt_bound(bopt, x, nbd, l, u, i)
use DataType
integer, parameter :: dp = kind(1.0d0)
integer::bopt,nbd(3*(N13+N2IN+N4OUT)),i
Real(dp)::x(3*(N13+N2IN+N4OUT))
Real(dp)::l(3*(N13+N2IN+N4OUT))
Real(dp)::u(3*(N13+N2IN+N4OUT))

!==========================    
    If  (bopt==1) then ! 1st setting for the bound atoms.
!------------
    if (x(i+2).le.-clt1*a2/4+0.1d0) then ! for the substrate simulation.
    nbd(i)=2 ;nbd(i+1)=2; nbd(i+2)=2 !limited displacement of atom x,y,z coords. (bound coordinates).
    l(i)=x(i)-0.05d0;  l(i+1)=x(i+1)-0.05d0  ;l(i+2)=x(i+2)+0.0d0  !the highest allowed value of atom coords.
    u(i)=x(i)+0.05d0;  u(i+1)=x(i+1)+0.05d0  ; u(i+2)=x(i+2)+0.05d0 !the lowest allowed  value of atom coords.
!-----------
    else if (x(i+2).gt.-a2/3-0.1d0.and.x(i+2).le.-0.1d0) then ! for the WL effect simulation.
    nbd(i)=0 ;nbd(i+1)=0; nbd(i+2)=2 !unbound atoms in xy plane; limited displacement in z planes.
    l(i+2)=x(i+2)-0.2d0 !the highest allowed z value of atom coords.
    u(i+2)=x(i+2)-0.1d0 !the lowest allowed z value of atom coords.
!------------   
    else if (x(i+2).gt.-0.1d0.and.x(i+2).le.a2/4+0.1d0) then ! for the WL effect simulation.
    nbd(i)=0 ;nbd(i+1)=0; nbd(i+2)=2 !unbound atoms in xy plane; limited displacement in z planes. 
    l(i+2)=x(i+2)-0.2d0  !the highest allowed z value of atom coords.
    u(i+2)=x(i+2)-0.1d0 !the lowest allowed z value of atom coords.
!------------       
    else if (x(i+2).gt.a2/4+0.1d0.and.x(i+2).le.clt2*a2+0.1d0) then ! for the WL effect simulation.
    nbd(i)=0 ;nbd(i+1)=0; nbd(i+2)=2 !unbound atoms in xy plane; limited displacement in z planes. 
    l(i+2)=x(i+2)-0.05d0  !the highest allowed z value of atom coords. 
    u(i+2)=x(i+2)+0d0 !the lowest allowed z value of atom coords.
!------------       
    else !for unbound atoms in the QD+cap domain
    nbd(i)=0; nbd(i+1)=0; nbd(i+2)=0 
    end if 
!-------------
    End If
!=========================
!
    If  (bopt==2) then ! 2nd setting for the bound atoms.
!------------
    if (x(i+2).le.-clt1*a2/4+0.1d0) then ! for the substrate simulation.
    nbd(i)=2 ;nbd(i+1)=2; nbd(i+2)=2 !limited displacement of x,y,z atom coords. (bound coordinates).
    l(i)=x(i)-0.05d0; l(i+1)=x(i+1)-0.05d0; l(i+2)=x(i+2)+0.0d0  !the lowest allowed value of x=x(i),y=x(i+1),z=x(i+2) Cartesian coordinates of the Mod(i-1,3)+1 atom label (index).
    u(i)=x(i)+0.05d0; u(i+1)=x(i+1)+0.05d0; u(i+2)=x(i+2)+0.05d0 !the highest allowed...
!-----------
    else if (x(i+2).gt.-0.1d0.and.x(i+2).le.clt2*a2+0.1d0) then ! for the WL effect simulation.
    nbd(i)=0 ;nbd(i+1)=0; nbd(i+2)=2 !unbound atoms in xy plane; limited displacement in z planes. 
    l(i+2)=x(i+2)-0.05d0  !the lowest allowed z value of atom coords. 
    u(i+2)=x(i+2)+0d0 !the highest allowed z value of atom coords.
!------------       
    else !for unbound atoms in the QD+cap domain
    nbd(i)=0; nbd(i+1)=0; nbd(i+2)=0 
    end if 
!-------------
    End If
!==========================
End Subroutine opt_bound
    
