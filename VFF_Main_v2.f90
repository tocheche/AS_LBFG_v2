    
    ! AS_LBFG_v2 calculates the atoms coordinates in half-torus, cone, truncated cone, pyramidal, core-shell quantum dots (QDs) heterostructures (QD=1,2,21,3,4).
    ! In the comments below subroutine is abreviated by SR.
    ! GaSb/GaAs, InAs/GaAs, Si/C are commented as examples of hetero-structures to easy the code reading.
    !   For GaSb/GaAs At1 is Ga, At2 is Sb, At3 is As; For InAs/GaAs At1 is As, At2 is In, At3 is Ga; For Si/C At1&At2 is Si, At3&At4 is C.    
!------------------------------------------------------------------------------------------
!    
! AS_LBFG_v2 - CONTENT
!* SR elastic_ct stores the bond-stretching force constants & bulk bond lengths used in U, DU calculus (GaSb/GaAs, Si/C). 
!* SR FCCbox(..) generates the FCC crystal in a box with 8*na=8*(2*nx+1)*(2*ny+1)*(nzmin+nzmax+1) total # atoms (of type At1 & At2) and latt_ct=1. The z=0 plane contains 2(2*nx+1)*(2*ny+1) At1 of coords (x,y,z) with each atom in the center of square of edge size equal to 1 and one of the atoms has coords (0,0,0). Each xy plane contains (2*nx+1)*(2*ny+1) atoms.
!* In VFF_Main_v2.f90:
    !At11 & At21 generate At1 and At2 coords; At11 & At21 are called in loop  to build two FCC crystal boxes, box1 with latt_ct=a1 (e.g. GaAs for matrix-OUT) and box2 with latt_ct=a2 (e.g. GaSb for QD-IN);
    !QD=1,2,3: latt_ct=a1 (OUT, for the substrate, GaAs), latt_ct=a2 (IN, for QD (GaSb) & capping(GaAs) in THE MODEL USED)
    !QD=4: latt_ct=a1 OUT (shell), latt_ct=a2 IN(core).
    ! At(1,:,:),At(2,:,:) rename At11,At21, respectively.
!* SelectQDT(..) selects the QD type and establishes if the atom is IN or OUT QD+WL(QDWL) for QD=1,2,3, or simply QD for QD=4. It specifies if the atom At(ss,n,..) from the crystal (with latt_ct=1) located either in box1 (coords index 1) or in box2 (coords index 2), is IN or OUT by log1, log2, respectively.
!* In VFF_Main_v2.f90: in CountAtoms loop the atoms are distributed in the QDWL, substrate, and capping domain for QD=1,2,3, or IN and OUT for QD=4.
!* SR QDWLB6(..) puts FCC atoms in the box (matrix) which encloses the QDWL and generates files 38,39,140,160,180. It uses SRs QDHTF(..),  QDConeF(..), QDPyramidF(..), CSQDF(..) which is the geometrical function shaping the QDWL or QD. File 140 collects IN & OUT QD  At1(Ga) coords; File 160  collects IN QD the At2(Sb) coords; File 180 collects OUT QD At3(As) coords.
!* After CountAtoms loop is performed we obtain the atoms coords and # of each kind of atoms (N13-Ga, Nat2-Sb, N4OUT-As). 
!* SR IniConfig6(x) writes the files T1,T2,T3, rearrange files 140,160,180 generated by QDWLB6(..), and writes files 201,202 in which the initial configuration (Cartesian atoms coords), x, is stored.
!* SR Sortx123(x) finds distances and indexes (ordering label(OL)) of each atom and its first 4 neighbor atoms (NAs); it uses the atom coords x generated by SR IniConfig6(x).
    !It contains SR Sorting4X which finds the distances between each pricipal atom (PA) and its first 4 NAs, and the OL of each atom (as a PA) and the OLs of its first 4 NAs.
    !Thus for each PA one finds its NAs indexed as the sequence 1,2,..j. Then one computes all the j distances PA-NA, one ascendingly sorts these j distances and the the sequence 1,2,..j is reordered as function of increasing PA-NA distance (one keeps only the first 4 indexes obtained after the sorting). Example: one finds that a PA has 7 NAs indexed 1,2,..7; after the sorting as function of PA-NA distance the new NAs order for the first 4 distances is, for example, 1,4,2,6 (new OLs). 
    !For At1(Ga): File temp240 contains: the first 4 distances PA-NA (given PA as Ga) & the new OLs denoted as Nb1.
    !For At2(Sb): File temp260 contains: the first 4 distances PA-NA (given PA) & the new OLs denoted as Nb2.
    !For At3(As): File temp280 contains: the first 4 distances PA-NA (given PA) & the new OLs denoted as Nb3.
    !File 111 provides the indexes of each atom and its first 4 NAs for the initial config x. 
!* SR abdiddpn(..) rearranges the data in files temp240.dat,temp260.dat,temp280.dat and collects alfa(q),beta(q),did(q),dpn(q) in file 551.
!* SR EnConfig(..) uses files 111, 551, atom coords (written if wished in file 501) to  calculate the elastic energy Uel and the gradient DUel for a given configuration x.
!* Routine setulb (..) minimizes the elastic energy Uel for a given initial config x contained by file 202. It requires Uel and its gradient DUel in each step of the minimization (Uel and its gradient DUel are computed by SR EnConfig(..)). 
!-----------------------------------------------------------------
!
!COMMENTS
    !1. The neighbors founds in the initial atomic configuration are labeled and these labels are kept during the the LBFG minimization procedure used in this code.
    !2. oneAPI Fortran https://software.intel.com/content/www/us/en/develop/articles/oneapi-standalone-components.html  
    
    !?3. VFF-Main should minimize the elastic energy (Uel) by increasing the # atoms until the strain field convergence is obtained.
    !?4. The code can run (for strtmax.gt.1) for various radii of the cutting-off sphere enclosing the neighbors of an atom.  
!*************************************************

    PROGRAM MAIN
!
!I. Declaration of variables & Definitions
use DataType !for a1,a2,Rt,Rq,ratio,N13,..., N24 type of data
integer, parameter       :: dp = kind(1.0d0)
    real(dp)             :: x1, y1, z1, x2, y2, z2 ! Cartesian atom coordinates
    real(dp)             :: WL, adf(5), RCS 
    integer              :: na, nt, n
    logical              :: log1, log2
    real(dp), allocatable:: At1(:,:), At2(:,:), At11(:,:), At21(:,:), At(:,:,:)
    real(dp), allocatable:: T1(:,:), T2(:,:), T3(:,:), T23(:,:)
    real(dp), allocatable:: T240(:,:) !, T2(:,:), T3(:,:), T23(:,:)
    integer              :: i1IN, i1OUT, i1, i2, i3, ss
    real(dp), allocatable:: x(:), xiter(:)
    real(dp)             :: Uel, epsmch, epscap
    real(dp), allocatable:: DUel(:,:), El_en(:), tollr(:)
    real                 :: time1,Testim2,Testim3
    integer              :: QD, bound, strt, ConvSteps, iter, iterr, iterrcount, bopt,celshape, CSQD
    real(dp)             :: Geom_par1, Geom_par2
!------------------------------------------------------- 
!I.2 Declarations for blas.f, lbfgsb_v2.f, linpack.f, timer.f:       
    integer, parameter   :: m = 5, iprint = 1
    real(dp), parameter  :: pgtol  = 1.0d-5
    real(dp)             :: factr
    character(len=60)    :: task, csave
    logical              :: lsave(4)
    integer              :: isave(44)
    real(dp)             :: f, ff
    real(dp)             :: dsave(29)
    integer, allocatable :: nbd(:), iwa(:)
    real(dp), allocatable:: l(:), u(:), g(:), wa(:)
    integer              :: i, Ntot,lbfg,sortcounter
    integer              :: xfile, fileno
!--------------------------------------------------------
!
    Open (15,file='VFF_v2.inp',status='old') !reading input parameters
    Read (15,*) QD
    Read (15,*) bound ! active command only for QD=1,2,3 (see the input file)
    Read (15,*) bopt
    Read (15,*) clt1, clt2
    Read (15,*) Rt, Rq
    Read (15,*) Rc, h
    Read (15,*) Rct, ht, hct
    Read (15,*) RCO 
    Read (15,*) a1, a2
    Read (15,*) Nw, nx, ny, nzmin, nzmax
    Read (15,*) C11out, C12out, C11in, C12in 
    Read (15,*) radius, radiuscut1 !, radiuscut2 !, radiuscut3, radiuscut4, radiuscut5
    Read (15,*) strtmax, factr
    Read (15,*) latshift1, latshift2
    Read (15,*) htj ! active command only for core-shell (see the input file)
    Read (15,*) celshape ! active command only for core-shell (see the input file)
    Read (15,*) CSQD ! active command only for core-shell (see the input file)
    Close (15)  
    continue
!------------------------------------------------------------
!
    call elastic_ct ! setting the elasticity parameters
!
na=(2*nx+1)*(2*ny+1)*(nzmin+nzmax+1) !a quarter of the total # atoms of one kind (At1 or At2) in an FCC box; there is a common kind of atom of QD and matrix (Ga).    
nt=4*na !there are nt At1 & nt At2 in an FCC box (8*na total # atoms (At1 & At2)).
WL=Nw*a1/2  !WL thickness; Nw # monolayers (Nw=1 means one layer of Sb/In at z=-0.25*a2 & one layer of Ga/As at z=-0.5*a2, for GaSb/InAs QD).
adf=nzmax*a2*(/0.7d0,0.75d0,0.8d0,0.85d0,0.9d0/) ! used to accomodate the diagonal strain tensor components to the bulk value (0 value). 
!-----------------
!
if (QD==4) then ! to ensure the existence of atoms in the shell of core-shell QD
RCS=nx*a1
    if (RCO .ge. RCS) then
    write(6,*) "Input RCO smaller or nx=ny=nzmin=nzmax larger to have atoms in the shell of core-shell QD"
    STOP
    end if
end if
!******************************************************************
    
!II.1 INITIAL CONFIGURATION
!   FCC symmetry atoms are placed in a rectangular box (matrix) which encloses QDWL.
!
if (QD==1) Print*,'Relaxation for Half-Torus QD is starting:'
if (QD==2) Print*,'Relaxation for Cone QD is starting:'
if (QD==21) Print*,'Relaxation for Truncated Cone QD is starting:'
if (QD==3) Print*,'Relaxation for Pyramid QD is starting:'
if (QD==4) Print*,'Relaxation for Core-Shell QD is starting:'
!Print*,'           * * *'
Print*,''
Print*,'A. START the calculus of initial atomic configuration x (file 202):'

!
Call FCCbox(nx,ny,nzmin,nzmax) !generates At1 & At2 coords in a box with latt_ct=1 (files 11,12).
    Allocate (At1(na,20),At2(na,20),At11(nt,5),At21(nt,5))
    Allocate (At(2,nt,5))
!
    Rewind(11)  !11 collects index(ordering label), bulk atom coords (for latt_ct=1) for At1; z layers at 0,-0.5,0.5,...
    Read (11,*) At1
    !If stack overflow in Visual Studio: Project>NameProperties...>Fortran>Optimization>Heap Arrays>0.
    At11=RESHAPE(At1, (/nt,5/), ORDER = (/2,1/)) !Write array At1 in C style: At11(n,1)=n(index), At11(n,2)=1(Atom type), At11(n,3)=x1(n), At11(n,4)=y1(n), At11(n,5)=z1(n).
    At(1,:,:)=At11 !used in the below CountAtoms loop to place At1.  
!
    Rewind(21) !21 collects: type of At2 (As or Sb), index(ordering), FCC bulk atom coords, latt_ct=1; z layers at -0.25,0.25,-0.75,0.75,...
    Read (21,*) At2
    At21=RESHAPE(At2, (/nt,5/), ORDER = (/2,1/)) !Rearrange array At2 in C style: At21(n,1)=n(index), At21(n,2)=2(Atom type), At21(n,3)=x2(n), At21(n,4)=y2(n), At21(n,5)=z2(n) for latt_ct=1.
    At(2,:,:)=At21 !used in the below CountAtoms loop to place At2 (Sb orAs).
!===========================================================================
!
!Next, At11 & At21 are used to build FCC bulk atom coordinates in 2 boxes; in box1 one sets latt_ct=a1, while in box2 one sets latt_ct=a2. Box2 contains in the origin (0,0,0) of Cartesian system At1 (Ga), WL, QD, and capping. Box1 contains underneath WL, that is the substrate (GaAs).
!
! Below loop CountAtoms sets the initial FCC atom coords in the substrate (latt_ct=a1, GaAs), QDWL(latt_ct=a2, GaSb), and capping(latt_ct=a2, GaAs). In CountAtoms one generates the 2 boxes with 8*na atoms in each of them. SRs SelectQDT(..) & QDWLB6(..) inside the loop are used to select/count the atoms from the boxes as belonging to the substrate (for atoms from box1), or QD, WL, and capping (for atoms from box2).
!
i1IN=0; i1OUT=0; i1=0; i2=0; i3=0;
!===========================================================================
!
CountAtoms: DO ss=1,2 ! ss=1 type of At1 (Ga); ss=2 type of At2 (As or Sb);! outer loop of counting atoms
InnerLoop:     Do n=1, nt  !n counts atoms; there are nt At1 and nt At2; inner loop of counting atoms
!"""""""""""""""""""""""""""""""""""""
!
QD123: If (QD.ne.4) then ! atomic positions for half-torus, cone, pyramid
!---------------------------------***
!
!Conditions for atom positions in WL:
!-------------------        
WL1: If (ss==1 .and. At(ss,n,5)==0d0) then !for At1(Ga) WL atoms.
!    x1=a2*At(ss,n,3); y1=a2*At(ss,n,4); z1=a1*At(ss,n,5) !x2,y2,z2 coords in z=-a2/2 for Ga in WL. !P,P2-v1
        x1=a2*At(ss,n,3); y1=a2*At(ss,n,4); z1=a2*At(ss,n,5) !x2,y2,z2 coords in z=-a2/2 for Ga in WL. !P2-v2
!----    
        iF (QD==3.and.(dabs(x1).le.(nx-latshift1)*a1).and. (dabs(y1).le.(ny-latshift1)*a1)) then !for pyramid QD 
                i1IN=i1IN+1 !counts #Ga located only in WL.
                i1=i1+1 !counts #Ga located in WL.        
                write(38,*)  x1,y1,z1 !for in WL At1 (Ga)
                write(140,*) x1,y1,z1 !for At1 (Ga)
        end iF
!----           
        iF ((QD==1.or.QD==2.or.QD==21).and.(dsqrt(x1**2+y1**2).le.nx*a1)) then !for cone & hemi-torus QD 
                i1IN=i1IN+1 !counts #Ga located only in WL.
                i1=i1+1 !counts #Ga located in WL.        
                write(38,*)  x1,y1,z1 !for in WL At1 (Ga)
                write(140,*) x1,y1,z1 !for At1 (Ga)
        end iF
!---- 
        go to 2347 !after atom positions set up cycle out InnerLoop
End If WL1
!-------------------
!
WL2: If (ss==2 .and. At(ss,n,5)==-0.25d0 ) then !for At2 (Sb) WL atoms.
 !       x2=a2*At(ss,n,3); y2=a2*At(ss,n,4); z2=a1*At(ss,n,5) !x2,y2,z2 coords in z=-a2/4 layers for Sb in WL.in G, H,M        
  !              x2=a2*At(ss,n,3); y2=a2*At(ss,n,4); z2=a2*At(ss,n,5) !x2,y2,z2 coords in z=-a2/4 layers for Sb in WL. !P2-v2    
                 x2=a2*At(ss,n,3); y2=a2*At(ss,n,4); z2=(a1+a2)/2*At(ss,n,5) !x2,y2,z2 coords in z=-a2/4 layers for Sb in WL. !P2-v3   
  !  x2=(a1+a2)/2*At(ss,n,3); y2=(a1+a2)/2*At(ss,n,4); z2=(a1+a2)/2*At(ss,n,5) !x2,y2,z2 coords in z=-a2/4 layers for Sb in WL. !P2-v4   
!----    
        iF (QD==3.and.(dabs(x2).le.(nx-latshift1)*a1).and.(dabs(y2).le.(ny-latshift1)*a1))  then !for pyramid QD
            i2=i2+1 !counts #Sb located in WL.
            write(160,*) x2,y2,z2 ! At2 (Sb) coords.
        end iF
! 
        iF ((QD==1.or.QD==2.or.QD==21).and.dsqrt(x2**2+y2**2).le.nx*a1) then ! for cone & hemi-torus QD       
            i2=i2+1 !counts #Sb located in WL.
            write(160,*) x2,y2,z2 ! At2 (Sb) coords.
        end iF
!---- 
        go to 2347  !after the atom positions is set up cycle out InnerLoop     
End If WL2
!-------------------
!
WL3: If (ss==1 .and. At(ss,n,5)==-0.5d0) then !for At1(Ga) WL atoms.
        x1=(a1+a2)/2*At(ss,n,3); y1=(a1+a2)/2*At(ss,n,4); z1=a1*At(ss,n,5) !x2,y2,z2 coords in z=-a2/2 for Ga in WL.
!----    
        iF (QD==3.and.(dabs(x1).le.(nx-latshift1)*a1).and. (dabs(y1).le.(ny-latshift1)*a1)) then !for pyramid QD 
            i1IN=i1IN+1 !counts #Ga located only in WL.
            i1=i1+1 !counts #Ga located in WL.        
            write(38,*)  x1,y1,z1 !for in WL At1 (Ga)
            write(140,*) x1,y1,z1 !for At1 (Ga)
        end iF
!      
        iF ((QD==1.or.QD==2.or.QD==21).and.dsqrt(x1**2+y1**2).le.nx*a1) then        
            i1IN=i1IN+1 !counts #Ga located only in WL.
            i1=i1+1 !counts #Ga located in WL.        
            write(38,*)  x1,y1,z1 !for in WL At1 (Ga)
            write(140,*) x1,y1,z1 !for At1 (Ga)
        end iF
!----
        go to 2347 !after the atom position is set up cycle out InnerLoop
End If WL3
!-------------------
! End of Conditions for atom positions in WL
!--------------------------------------------------------------------------
!
! Conditions for atom positions in two FCC boxes (excluding WL):
!-------------------
!For the 2 boxes (of latt_ct=a1, a2) without the WL domain:
IF ((At(ss,n,5).ne.-0.25d0) .and. (At(ss,n,5).ne.-0.5d0) ) then !for atoms of box excepting WL atoms
!Next command sets x1,y1,z1 as coords in z1-layers at 0,a1/2,-2a1/2,2a1/2,... for Ga (if ss=1) and in z1-layers at a1/4,-3a1/4,3a1/4,... for Sb or As (if ss=2):
    x1=a1*At(ss,n,3); y1=a1*At(ss,n,4); z1=a1*At(ss,n,5)
!Next command sets x2,y2,z2 coords in z2-layers at 0,a2/2,-2a2/2,2a2/2,... for Ga (if ss=1) and in z2-layers at a2/4,-3a2/4,3a2/4,... for Sb or As (if ss=2):
    x2=a2*At(ss,n,3); y2=a2*At(ss,n,4); z2=a2*At(ss,n,5)
End If
! End of Conditions for atom positions in two FCC boxes (excluding WL)
!-------------------
!
End If QD123 ! atom positions for pyramid, cone, half-torus
!---------------------------------***
!
!---------------------------------***
!
QD4: If (QD==4) then !atom positions for core-shell
    x2=a2*At(ss,n,3); y2=a2*At(ss,n,4); z2=a2*At(ss,n,5)
        if (CSQD==1) then ! one sets a1 equal to a2 (lattice contant of QD in both core and shell)
        x1=x2; y1=y2; z1=z2;
        end if
        if (CSQD.ne.1) then ! one sets different lattice constants for core and shell
            x1=a1*At(ss,n,3); y1=a1*At(ss,n,4); z1=a1*At(ss,n,5)
        end if    
!       
    iF(celshape==0) then ! for cubic shape of shell of core-shell
    go to 2346
    enD iF
!    
    iF(celshape==1) then !for spherical shell of core-shell
        if(((x1**2+y1**2+z1**2).le.RCS**2) .and. ((x2**2+y2**2+z2**2).le.RCS**2)) then ! makes the shell spherical
            go to 2346
        else
            go to 2347    
        end if
    enD iF
End If QD4
!---------------------------------***
!"""""""""""""""""""""""""""""""""""""
!
!SelectQDT(..) selects QD type & specifies if the atoms of coords (x1,y1,z1), (x2,y2,z2) are located IN or OUT QDWL by log1,log2,respectively:
2346    Call SelectQDT(x1,y1,z1,x2,y2,z2,log1,log2,QD)
!     
!QDWLB6(..) counts and writes atom coords for IN and OUT QDWL; i1,i2,i3 are loop counters for atoms Ga,Sb,As; files 140,160,180 store unrelaxed atom coords:        
    Call QDWLB6(QD,ss,x1,y1,z1,x2,y2,z2, log1,log2,i1IN,i1OUT,i1,i2,i3,WL,RCS) 
!""""""""""""""""""""""""""""""""""""""
!
2347 End Do InnerLoop ! inner loop of counting atoms
End DO CountAtoms ! outer loop of counting atoms
!=================================================
!QD(At1,At2)/matrix(At3,At4)
    N1IN=i1IN !# At1(Ga) IN QDWL; ss=1.
    N3OUT=i1OUT !# At3(Ga) OUT QDWL; ss=1.
    N13=i1 !#At1-IN + #At3-OUT(Ga); i1=i1IN+i1OUT
    N2IN=i2 !#At2(Sb) IN QDWL; ss=2.
    N4OUT=i3 !# At4(As). N13+N2IN+N4OUT = 8*na.
    N24=N2IN+N4OUT
!
!******************************************************************
!
If(QD==4) then
Print*, 'Core-Shell QD'
Print*, '#atom_kind1-IN QD = ',N1IN
Print*, '#atom_kind2-IN QD = ',N2IN
Print*, '#atom_kind3-OUT QD = ',N3OUT
!Print*, 'Check: #atom1-IN + #atom3-OUT = ',N13
Print*, '#atom_kind4-OUT QD = ',N4OUT
End If
!
If(QD.ne.4) then
Print*, 'WL has thickness (Angstr�m)= ', WL
Print*, '#atom_kind1-IN(Ga/As) QDWL = ',N1IN
Print*, '#atom_kind2-IN(Sb/In) QDWL = ',N2IN
Print*, '#atom_kind3-OUT(Ga/As) QDWL = ',N3OUT
Print*, '#atom_kind4-OUT(As/Ga) QDWL = ',N4OUT
Print*, '#atom_kind1-IN(Ga/As) + #atom_kind3-OUT(Ga/As) = ',N13
End If
!stop
!-----------------------------------
epscap=(a2-a1)/a2
If (QD==1)  then            
    Geom_par1=Rt; Geom_par2=Rq
end if
If (QD==2 .or. QD==3) then
    Geom_par1=Rc; Geom_par2=h
end if
If (QD==21) then
    Geom_par1=Rct; Geom_par2=ht; Geom_par3=hct ! for truncated cone
end if
If (QD==4) then           
    Geom_par1=RCO; Geom_par2=RCS
end if
!* write(124,*) QD, WL, N1IN, N3OUT, N2IN, N4OUT, a1, a2, Geom_par1, Geom_par2, adf,epscap
continue
!STOP
!======================================================
!
Ntot=3*(N13+N2IN+N4OUT)
    allocate (nbd(Ntot), x(Ntot),xiter(Ntot),l(Ntot), u(Ntot), g(Ntot) )  !for lbfgsb.f
    allocate (iwa(3*Ntot)) !for lbfgsb.f
    allocate (wa(2*m*Ntot + 5*Ntot + 11*m*m + 8*m) ) !for lbfgsb.f
    allocate (radiuscut(1)) ! (radiuscut(5)) can be used if in VFF_v2.inp one sets 5 values of vector radiuscut.
    allocate (T1(N13,3),T2(N2IN,3),T3(N4OUT,3),T23(N2IN+N4OUT,3))
    allocate (T240(N13,8))
!
!if stackoverflow in visual studio, see https://stackoverflow.com/questions/44568453/meaning-of-open-shared-in-fortran
Call IniConfig6(x)
!
!radiuscut=(/radiuscut1,radiuscut2,radiuscut3,radiuscut4,radiuscut5/) !Radii for searching the NA of a PA inside a sphere of radius=radiuscut(..) 
radiuscut=(/radiuscut1/) !Radius for searching the NA of a PA inside a sphere of radius=radiuscut1.
    !radius=4 is the initial searching radius centered on PA to find he 1st order NAs. radiuscut1,2.. is distanta PA to the 1st order NAs, a little larger than the ideal structure.
!
! file 124 would be used by code AS_Calc_v2 for the strain field calculus
write(124,*) QD, WL, N1IN, N3OUT, N2IN, N4OUT, a1, a2, Geom_par1, Geom_par2, Geom_par3, adf,epscap,radiuscut,C11out,C12out,C11in,C12in,htj
!STOP
allocate(DUel(N13+N2IN+N4OUT,3)) !allocation for the Uel gradient.
!**************************************************
!
!III. ELASTIC ENERGY MINIMIZATION
!
sortcounter=0 !if sortcounter remains zero then the PA and neighbor indexes in the sorting from initial configurationare are not changed during minimization.
Print*,'           * * *'
Print*,''
!=====================================================
!
!III.1. START MINIMIZATION
Print*,'B. START SOLVING THE MINIMIZATION PROBLEM:'
Minimization: DO strt=1,strtmax ! strtmax is max # of cycles of minimization (=1 for the results in the manuscript)
!strt is the counter for the #strt loop in the minimization process. For strt=q>1 the output file 60(q-1) in the x-th minimization loop becomes with the file 202 which contains the initial coords in the x-th minimization loop (e.g., if strt=2 then 601 becomes the file 202 and 2nd minimization process starts with NAs inside a sphere of radius=radiuscut2). For strt=1 the initial file 202 is identical to 201; at the end of the 'Minimization Process Loop (MPL)', file 202 contains the new x coords after the #strt MPL and it becomes the input for the next MPL step (if any).
!
    Rewind(202)   
        do i=1,3*(N13+N2IN+N4OUT),3
    read(202,*) x(i),x(i+1),x(i+2)   
    end do    
    Print*,'--------------------------------------'
    Print*,''  
    Print*, 'MINIMIZATION PROCESS #', strt
    Print*,''
!
!********************************************   
i=0
DO i=1, Ntot, 3
!========================== 
If(bound==1) then !Displacement of some PAs is limited to simulate the substrate and WL. There are 2 options for location of bound atoms which are set with below 'bopt'.
    Call opt_bound(bopt, x, nbd, l, u, i) !bound atom low and up limits are defined 
End If
!==========================
If(bound==0) then !unbound atoms (no restriction for the PA on the displacement limit)
    nbd(i)=0; nbd(i+1)=0; nbd(i+2)=0   
End If 
!========================== 
END DO 
!********************************************
!
allocate(tollr(0:1000))
allocate (El_en(0:1000))
!
xfile=600+strt;
Print*,'The atomic coordinates x for minimum Uel are enclosed by file',xfile,' .'
!
    task = 'START' !START the iteration by initializing the task.
    lbfg = 0 !Counter for # of the loops in the minimization process (for given strt).
    iterr = 0
    El_en = 0
    iterrcount=0
!   
main_do: DO while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or.task.eq.'START')         
    call setulb (Ntot,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint,&
                 csave,lsave,isave,dsave )
!SR setulb(..) initiates the L-BFGS-B code.
!In the 1st step: in setulb(..) START->FG_START, x=x_ini. EnConfig(..) finds f & g for x_ini.
!In the 2nd step: setulb(..) computes x=X_NEW; EnConfig(..) computes new f & g (f=Uel, g=DUel in EnConfig)//, etc...
    main_if: IF (task(1:2) .eq. 'FG') then !In 1st step: task->FG_START. In next steps: task->FG_LNSRCH or task->X_NEW
    lbfg=lbfg+1 !lbfg counts loop # in the strt-th minimization process.
!    
    if (lbfg==1) then
        rewind(76)
        read(76,*) epsmch !epsmch is the machine precision, which is automatically generated by the code (in lbfgsb_v2.f)
    end if
!
!--------   
    if (lbfg==1)then
        call CPU_TIME (time1)
    Testim2=time1
    else if (lbfg==2)then
    call CPU_TIME (time1)
    Testim3=time1
    endif
!    
if (lbfg==2)then
    Print*,'===================================================='
    Print*,'INFO:'
    ! ConvSteps is the approximate # of iterations 
    if (factr.ge.1.0d+10.and.factr.lt.1.0d+11)then !factr is accuracy parameter, see VFF_v2.inp.
        ConvSteps=1000*strtmax
    elseif (factr.ge.1.0d+11.and.factr.lt.1.0d+12)then
        ConvSteps=300*strtmax
    elseif (factr.ge.1.0d+12.and.factr.lt.1.0d+13)then 
        ConvSteps=200*strtmax
    elseif (factr.ge.1.0d+13.and.factr.lt.1.0d+15)then 
        ConvSteps=50*strtmax 
    end if
    Print*,' The estimated total time calculus is ',(Testim3-Testim2)*ConvSteps,' s'
    Print*,'===================================================='    
end if
continue
!-------------------------------
!
    sortcounter=sortcounter+1
    if (lbfg==1 .and. sortcounter==1)then
        Call Sortx123(x,lbfg) !finds the 1st 4NAs of each PA.
        Call abdiddpn(radiuscut,strt) !generates elastic constants and distances PA-NA
        continue
    end if
!        
if(lbfg.gt.1)then
    rewind(78)
    read(78,*) iterr, ff !iterr counts how many times ff=Uel has been counted.
    El_en(iterr)=ff
    close (78)
end if
!
if (iterr.ne.0 .and. iterr==10*NINT(iterr/10d0).and.iterrcount.ne.13) then !collecting the atomic coords after each 10 iterations
        tollr(iterr)=dabs(El_en(iterr)-El_en(iterr-1))/max(El_en(iterr),El_en(iterr-1),1d0) ! (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
        write(6000+iterr,3008) tollr(iterr)/epsmch
        write(6000+iterr,*) x        
        fileno=6000+iterr
        write(6,3009) fileno, El_en(iterr), tollr(iterr)/epsmch
3008    format('tollr(lbfg)/epsmch=',e12.6)
3009    format ('*** Atomic config from file ',i0.2, ' has Uel=',3p, d16.10, ' and Uel factor accuracy factr=',e12.6,/) 
    iterrcount=13 !avoids successively appending file 6000+iterr (e.g. appending to 1st created config 6010 file a 2nd config file labeled also 6010).
end if
if (iterr.ne.10*NINT(iterr/10d0)) iterrcount=0 !iterrcount is changed back to 0 to avoid successively appending file 6000+iterr.
!        
   Call EnConfig(QD,x,Uel,DUel,lbfg,radiuscut,strt)
!         
call CPU_TIME (time1)
Print*,'Total time for computing neighbors, Uel, grad(Uel) is',time1,'s'
!
    f=Uel !elastic energy
!
continue
    do i=1,Ntot,3
        ii=i/3+1
        g(i)=DUel(ii,1);g(i+1)=DUel(ii,2);g(i+2)=DUel(ii,3) !elastic energy gradient
    end do 
    
open(unit=707,file='temp707.dat',status='replace') 
    do i=1,Ntot-1,3
        write(707,*) x(i),x(i+1),x(i+2)
    end do
close(707)
continue
!
    END IF main_if   
END DO main_do
!====================================
continue
!
!Files 601,602, 603 collect x coords after 1st,2nd,and 3rd 'Minimization' (strt=1,2,3,respectively for strtmax=3, e.g.). File 201 contains initial x coords.
    close(202)
    open(unit=707,file='temp707.dat',status='old') 
        do i=1,3*(N13+N2IN+N4OUT),3 
        read(707,*) x(i),x(i+1),x(i+2)
        write(202,*) x(i),x(i+1),x(i+2) !prepares file to to be read in 'Minimization' cycle.
        write(600+strt,*) x(i),x(i+1),x(i+2)    
    end do
    close(707)
    write(6,3011) 600+strt, factr
3011 format('  File', 1x, i3,' file encloses relaxed atoms coordinates; limit factr accuracy of Uel is ', d10.3,'.')
continue
! 
END DO Minimization
!=====================================================
!    
if(QD==1)then  !Hemi-Torus
        write(6,3004) Rt,Rq,clt1,clt2
        3004 format ('Calculus done for Hemi-Torus QD with',3x,'Rt =',(1x,d10.3),',   Rq =',(1x,d10.3),',   clt1 =',(1x,d10.3),',   clt2 =',(1x,d10.3)) 
elseif (QD==2)then  !Cone
        write(6,3005) Rc,h,clt1,clt2
3005    format (/'Calculus done for Cone QD',3x,'Rc =',(1x,d10.3),',   h =',(1x,d10.3),',   clt1 =',(1x,d10.3),',   clt2 =',(1x,d10.3)) 
elseif (QD==21)then  !Truncated Cone
        write(6,30051) Rct,ht,clt1,clt2
        30051 format (/'Calculus done for Truncated Cone QD',3x,'Rct=',(1x,d10.3),', ht=',(1x,d10.3),', hct=',(1x,d10.3),', clt1 =',(1x,d10.3),',   clt2 =',(1x,d10.3))
elseif (QD==3)then  !Pyramid
        write(6,3006) 2*Rc,h,clt1,clt2 !base=2Rc
        3006    format ('  Calculus done for Pyramid QD with',3x,'base =',(1x,d10.3),',   h =',(1x,d10.3),',   clt1 =',(1x,d10.3),',   clt2 =',(1x,d10.3))
elseif (QD==4)then  !Core-Shell
        write(6,3107) RCO, RCS 
        3107 format ('  Calculus done for Core-Shell QD with',3x,'RCO =',(1x,d10.3),',   RCS =',(1x,d10.3)) 
end if
!
if(bound==0)then
    Print*,' and the atoms are unbound.'
elseif (bound==1)then
    Print*,' and some atoms are bound.'
end if
!
Write(6,3007)
     3007 format (/,'  AS_LBFG_v2 finished.')
!
    End Program Main    