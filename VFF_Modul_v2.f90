!***************************************************************    
    Module DataType
    implicit none 
integer,parameter     :: dpp = kind(1.0d0)
real(dpp),save        :: a1,a2,latshift1,latshift2
real(dpp),save        :: ratio,radius,radiuscut1 !,radiuscut2 ,radiuscut3,radiuscut4, radiuscut5
integer,save          :: strtmax, htj
real(dpp),allocatable :: radiuscut(:)
real(dpp),save        :: Rt,Rq,Rc,h,RCO,clt1,clt2, Rct,ht,hct
real(dpp),save        :: ac12,ac13,bc12,bc13,dc12,dc13
real(dpp),save        :: C11out,C12out,C11in,C12in
integer,save          :: N1IN, N2IN, N3OUT, N4OUT, N13, N24, Nw, nx, ny, nzmin, nzmax
    End Module DataType      
!*****************************************************************  