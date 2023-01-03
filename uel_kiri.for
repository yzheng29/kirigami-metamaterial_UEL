!	  User defined element subroutine (UEL) for planar kirigami
!
!     Material Properties to be passed to the UEL
!     c0   = PROPS(1) 
!     c1   = PROPS(2) 
!     c2   = PROPS(3)   
!     alpha = PROPS(4)
!     beta = PROPS(5)
!
!	  Number of UEL element, 
!	  Offset number of dummy element,  
!	  Number of integration points
!	  are used as global variables to transfer actuation angle from UEL
!	  to UVARM for visualization
!
!************************************************************************
      module global

      integer numElem,ElemOffset,numPts, err	                         ! Number of UEL element, Offset number of dummy element, Number of integration points
      parameter(numElem=1600, ElemOffset=10000, numPts=9)				     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Change to match the input file

      double precision :: globalSdv(numElem,numPts,4)

      end module global

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)
	 
      USE global
      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
	  
	  uvar(1) = globalSdv(NOEL-ElemOffset,NPT,1)
	  uvar(2) = globalSdv(NOEL-ElemOffset,NPT,2)
	  uvar(3) = globalSdv(NOEL-ElemOffset,NPT,3)
	  uvar(4) = globalSdv(NOEL-ElemOffset,NPT,4)
	  
      RETURN
      END SUBROUTINE UVARM
	  
!************************************************************************
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
      USE global
      INCLUDE 'ABA_PARAM.INC'
    
    
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)



    ! Local Variables
      integer      :: i,j,k, n_points,kint

      double precision  ::  xi(2,9)                                      ! Area integration points
      double precision  ::  wt(9)                             			 ! Area integration weights
      double precision  ::  N(9)                             			 ! shape functions
      double precision  ::  dNdxi(9,2)                       			 ! Shape function derivatives
      double precision  ::  dNdX(9,2)                        			 ! X undeformed coords
      double precision  ::  dXdxi(2,2)                       
      double precision  ::  dNdy(9,2) 									 ! y deformed coords	
      double precision  ::  dxidX(2,2)           
      double precision  ::  Jmap
	  
      double precision  ::  F(2,2)                            			 ! Deformation gradient
      double precision  ::  Finv(2,2)                         			 ! Inverse of deformation gradient
      double precision  ::  JF                                			 ! det(F)	  
	  
	  double precision  ::  aa                                       	 ! Actuation angle
	  double precision  ::  daa(1,2)
	  
      double precision  ::  Ru(18,1)                         	 		 ! Displacement residual
      double precision  ::  Ra(9,1)	                             		 ! Actuation angle residual
      double precision  ::  Kuu(18,18)               					 ! Tangent matrix
      double precision  ::	Kaa(9,9)
      double precision  ::	Kua(18,9)	
      double precision  ::	Kau(9,18)
	  
	  double precision  ::  W                                       	 ! Elastic energy
	  double precision  ::  dWda
      double precision  ::  stress(3)                         
      double precision  ::  D(4,4)
	  double precision  ::  ddWdaa
	  double precision  ::  Dua(3)	  
	  double precision  ::  dWdg(2,1)
	  double precision  ::  ddWdgg(2,2)
	  
	  
      double precision  ::  Bvgt(3,18)	  
      double precision  ::  Bmat(4,18) 
      double precision  ::  Nvec(9,1)	

      integer      :: n1, n2, m1, m2	  

	! Number of nodes and number of integration pionts
      if (NNODE == 3) n_points = 1              						 ! Linear triangle
      if (NNODE == 4) n_points = 4               						 ! Linear rectangle
      if (NNODE == 6) n_points = 4              						 ! Quadratic triangle
      if (NNODE == 8) n_points = 9               						 ! Serendipity rectangle
      if (NNODE == 9) n_points = 9             							 ! Quadratic rect


      call abq_UEL_integrationpoints2D(n_points, NNODE, xi, wt)	
	  
	
    ! Initialization  
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      ENERGY(1:8) = 0.d0
	  
	  Ru = 0.d0
	  Ra = 0.d0
	  Kuu = 0.d0
	  Kaa = 0.d0
	  Kua = 0.d0
	  Kau = 0.d0

    !------------------------------
    ! Loop over integration points
      do kint = 1, n_points
	  
        call abq_UEL_shapefunctions2D(xi(1:2,kint),NNODE,N,dNdxi)
        dXdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
		
        call abq_UEL_invert2d(dXdxi,dxidX,Jmap)
        dNdX(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidX)

    ! Caculate the deformation gradient
        do i = 1,2
            k = 3*(NNODE-1)+i
            F(i,1:2) = matmul(U(i:k:3),dNdX(1:NNODE,1:2))
            F(i,i) = F(i,i) + 1.d0
        end do

        
        call abq_UEL_invert2D(F,Finv,JF)
        dNdy(1:NNODE,1:2) = matmul(dNdX(1:NNODE,1:2),Finv)
		
	!   Calculate the actuation angle and its derivatives
		aa = 0.d0
		do i=1,NNODE
			aa = aa + U(3*i)*N(i)
		end do
		daa(1,1:2) = matmul(U(3:3*NNODE:3),dNdX(1:NNODE,1:2))
		
	!   Save the actuation angle for visualization
		globalSdv(JELEM,kint,1) = aa
		
	!	Use continuum field theory of planar kirigami	
        call mat_kirigami(PROPS(1:NPROPS),NPROPS,F,aa,daa,
     1       W, stress, dWda, D, ddWdaa, Dua, dWdg, ddWdgg )
		
		Energy(2) = W

	!   Save the true stress for visualization		
		globalSdv(JELEM,kint,2) = stress(1)/JF
		globalSdv(JELEM,kint,3) = stress(2)/JF
		globalSdv(JELEM,kint,4) = stress(3)/JF
		
        Bvgt = 0.d0
        Bvgt(1,1:2*NNODE-1:2) = dNdy(1:NNODE,1)
        Bvgt(2,2:2*NNODE:2)   = dNdy(1:NNODE,2)
        Bvgt(3,1:2*NNODE-1:2) = dNdy(1:NNODE,2)
        Bvgt(3,2:2*NNODE:2)   = dNdy(1:NNODE,1)		
              
        
        Bmat = 0.d0
        Bmat(1,1:2*NNODE-1:2) = dNdy(1:NNODE,1)
        Bmat(2,2:2*NNODE:2)   = dNdy(1:NNODE,1)
        Bmat(3,1:2*NNODE-1:2) = dNdy(1:NNODE,2)
        Bmat(4,2:2*NNODE:2)   = dNdy(1:NNODE,2)
		
        do i=1,NNODE
			Nvec(i,1) = N(i)
        end do

	! 	Calculate the residual vector
        Ru(1:2*NNODE,1) = Ru(1:2*NNODE,1)
     1   - matmul(transpose(Bvgt(1:3,1:2*NNODE)),stress(1:3))*
     2                                          wt(kint)*Jmap
	 
	    Ra(1:NNODE,1) = Ra(1:NNODE,1) 
     1      +( - Nvec(1:NNODE,1)*dWda  
     2	      - matmul(dNdX(1:NNODE,1:2),dWdg(1:2,1))   )
     3			 *wt(kint)*Jmap

	 
	!	Calculate the tangent matrix  
        Kuu(1:2*NNODE,1:2*NNODE) = Kuu(1:2*NNODE,1:2*NNODE)
     1  + matmul( transpose(Bmat(1:4,1:2*NNODE)),
     2     matmul(D,Bmat(1:4,1:2*NNODE)) )*wt(kint)*Jmap
	    
		Kaa(1:NNODE,1:NNODE) = Kaa(1:NNODE,1:NNODE)
     1  + ( spread(N(1:NNODE),2,NNODE)*spread(N(1:NNODE),1,NNODE)*ddWdaa
     2	  + matmul(   dNdX(1:NNODE,1:2),
     3	    matmul(ddWdgg, transpose(dNdX(1:NNODE,1:2)))   )    
     4     )*wt(kint)*Jmap

	 
	    Kua(1:2*NNODE,1:NNODE) = Kua(1:2*NNODE,1:NNODE)
     1  + matmul( transpose(Bvgt(1:3,1:2*NNODE)),
     2	         ( spread(Dua(1:3),2,NNODE)*spread(N(1:NNODE),1,3) ) )
     3		   *wt(kint)*Jmap
	 
		Kau(1:NNODE,1:2*NNODE) = Kau(1:NNODE,1:2*NNODE)
     1	+ matmul( ( spread(N(1:NNODE),2,3)*spread(Dua(1:3),1,NNODE) ),
     2	 			Bvgt(1:3,1:2*NNODE) )*wt(kint)*Jmap
	 
 
      end do
	! End Loop 
    !------------------------------

	! Assemble the residual vector
	  do i=1,NNODE
	      n1 = 3*(i-1)+1
		  n2 = 2*(i-1)+1
		  
		  RHS(n1:n1+1,1) = Ru(n2:n2+1,1)
		  RHS(n1+2,1) = Ra(i,1)
	  end do
	  
	! Assemble the tangent matrix
	  do i=1,NNODE
	      do j=1,NNODE
		      n1 = 3*(i-1)+1
			  n2 = 2*(i-1)+1
			  m1 = 3*(j-1)+1
			  m2 = 2*(j-1)+1
			  
			  AMATRX(n1:n1+1,m1:m1+1) = Kuu(n2:n2+1,m2:m2+1)
			  AMATRX(n1:n1+1,m1+2) = Kua(n2:n2+1,j)		  
			  AMATRX(n1+2,m1:m1+1) = Kau(i,m2:m2+1)			   
			  AMATRX(n1+2,m1+2) = Kaa(i,j)
		  end do
	  end do

	  

      END SUBROUTINE UEL
!************************************************************************

      subroutine mat_kirigami(props,nprops,F,aa,daa,
     1                    W, stress, dWda, D, ddWdaa, Dua, dWdg, ddWdgg )

       implicit none

       integer, intent(in)           :: nprops
       double precision, intent(in)  :: props(nprops)
       double precision, intent(in)  :: F(2,2)
       double precision, intent(in)  :: aa
       double precision, intent(in)  :: daa(1,2)
	  
       double precision, intent(out) :: stress(3)
       double precision, intent(out) :: D(4,4)
       double precision, intent(out) :: W
       double precision, intent(out) :: dWda
       double precision, intent(out) :: ddWdaa
       double precision, intent(out) :: Dua(3)
       double precision, intent(out) :: dWdg(2,1)
       double precision, intent(out) :: ddWdgg(2,2)
	   
       double precision :: c0, c1, c2, alpha, beta
	   double precision :: A(2,2), dA(2,2), ddA(2,2)
	   double precision :: Ainv(2,2), dAinv(2,2), ddAinv(2,2)
	   double precision :: JA, dJA, ddJA
       double precision :: Finv(2,2), JF
       double precision :: B(2,2), trB
       double precision :: R(2,2), JR, FdAinv(2,2)
	   double precision :: B0(2,2), B1(2,2), B2(2,2), B3(2,2)
	   double precision :: trB0, trB1, trB2, trB3
	   
	   double precision :: dWdu(2,2)
       double precision :: ddWduu(2,2,2,2)
       double precision :: ddWdua(2,2)
	    
       double precision :: Iden(2,2)
       integer :: i,n,k,m


    !  Obtain material parameters
       c0 = props(1)
       c1 = props(2)
       c2 = props(3)
       alpha = props(4)
       beta  = props(5)
	   
	   Iden=0.d0
	   Iden(1,1)=1.d0
	   Iden(2,2)=1.d0
	   
	!  Calculate matrix A
	   A = 0.d0
	   A(1,1) = cos(aa) - alpha*sin(aa)
	   A(2,2) = cos(aa) + beta*sin(aa)
	   
	   dA = 0.d0
	   dA(1,1) = -sin(aa) - alpha*cos(aa)
	   dA(2,2) = -sin(aa) + beta*cos(aa)
	   
	   ddA = 0.d0
	   ddA(1,1) = -cos(aa) + alpha*sin(aa)
	   ddA(2,2) = -cos(aa) - beta*sin(aa)
	   
	   call abq_UEL_invert2D(A,Ainv,JA)
	   
	   dAinv = 0.d0
	   dAinv(1,1) = (sin(aa) + alpha*cos(aa))
     1              /(cos(aa) - alpha*sin(aa))**2
       dAinv(2,2) = (sin(aa) - beta*cos(aa))
     1              /(cos(aa) + beta*sin(aa))**2
	 
       ddAinv = 0.d0
	   ddAinv(1,1) = 1.d0/(cos(aa) - alpha*sin(aa))
     1               + (2*(sin(aa) + alpha*cos(aa))**2)
     2                   /(cos(aa) - alpha*sin(aa))**3
       ddAinv(2,2) = 1.d0/(cos(aa) + beta*sin(aa))
     1               + (2*(sin(aa) - beta*cos(aa))**2)
     2                   /(cos(aa) + beta*sin(aa))**3
	   
	   dJA = - (sin(aa) + alpha*cos(aa))*(cos(aa) + beta*sin(aa))
     1       - (cos(aa) - alpha*sin(aa))*(sin(aa) - beta*cos(aa))
	 
	   ddJA = 2*(sin(aa) + alpha*cos(aa))*(sin(aa) - beta*cos(aa))
     1      - 2*(cos(aa) - alpha*sin(aa))*(cos(aa) + beta*sin(aa))
	   
	  
    !  Calculate Green Tensor
       call abq_UEL_invert2D(F,Finv,JF)      
       B = matmul(F,transpose(F))
       trB = B(1,1) + B(2,2)
	   
	   R = matmul(F,Ainv)
	   B0 = matmul(R,transpose(R))
	   trB0 = B0(1,1) + B0(2,2)
	   JR = JF/JA
	   
	   FdAinv = matmul( F,dAinv )
	   B1 = matmul( FdAinv,transpose(R)  )
	   B2 = matmul( FdAinv,transpose(FdAinv)  )
	   B3 = matmul( matmul(F,ddAinv),transpose(R) )
	   trB1 = B1(1,1) + B1(2,2)
       trB2 = B2(1,1) + B2(2,2)
       trB3 = B3(1,1) + B3(2,2)
	   
	!  Calculate elastic energy
       W = c0*(trB0/JR-2) + c0*(JR-1.d0)**2 
     1   + c1*aa**2 + c2*norm2(daa)**2
	   
	! Calculate first derivatives
	   dWdu = 2.d0*c0/JR*B0 -c0/JR*trB0*Iden
     1      + 2.d0*c0*(JR-1)*JR*Iden

       
	   stress(1) = dWdu(1,1)
	   stress(2) = dWdu(2,2)
	   stress(3) = dWdu(1,2)	 
	   
	   dWda = 2.d0*c0*trB1/JR +c0*trB0/JF*dJA
     1      - 2.d0*c0*(JF**2/JA**3 - JF/JA**2)*dJA
     2      + 2.d0*c1*aa
	

	! Calculate second derivatives
	   ddWduu = 0.d0
       do i=1,2
          do n = 1,2
             do k = 1,2
                do m = 1,2
                   ddWduu(i,n,k,m) = 2.d0*c0/JR*Iden(i,k)*B0(m,n)
     1                - 2.d0*c0/JR*( Iden(i,n)*B0(k,m)
     2                             + Iden(k,m)*B0(i,n) )
     3                + c0*trB0/JR*( Iden(i,n)*Iden(k,m)
     3                             + Iden(i,m)*Iden(k,n) ) 
     5                + 2.d0*c0*(2.d0*JR**2-JR)*Iden(i,n)*Iden(k,m)
     6                - 2.d0*c0*( JR**2-JR )*Iden(i,m)*Iden(k,n)

                enddo
            enddo
          enddo
       enddo	    
       
	   D = 0.d0
       D(1,1) = ddWduu(1,1,1,1)
       D(1,2) = ddWduu(1,1,2,1)
       D(1,3) = ddWduu(1,1,1,2)
       D(1,4) = ddWduu(1,1,2,2)
       D(2,1) = ddWduu(2,1,1,1)
       D(2,2) = ddWduu(2,1,2,1)
       D(2,3) = ddWduu(2,1,1,2)
       D(2,4) = ddWduu(2,1,2,2)
       D(3,1) = ddWduu(1,2,1,1)
       D(3,2) = ddWduu(1,2,2,1)
       D(3,3) = ddWduu(1,2,1,2)
       D(3,4) = ddWduu(1,2,2,2)
       D(4,1) = ddWduu(2,2,1,1)
       D(4,2) = ddWduu(2,2,2,1)
       D(4,3) = ddWduu(2,2,1,2)
       D(4,4) = ddWduu(2,2,2,2)
	   
	   ddWdaa = 2.d0*c0*trB2/JR +2.d0*c0*trB3/JR
     1          +4.d0*c0*trB1/JF*dJA
     2          +c0*trB0/JF*ddJA
     3          +2*c0*(3.d0*JF**2/JA**4 - 2.d0*JF/JA**3)*dJA**2
     4          -2*c0*(JF**2/JA**3 - JF/JA**2)*ddJA
     5          +2.d0*c1
	 
	   ddWdua = 2.d0*c0/JR*( B1 +transpose(B1) )
     1   	    + 2.d0*c0*B0/JF*dJA
     2          - 2.d0*c0*trB1/JR*Iden
     3          - c0*trB0/JF*dJA*Iden
     4          + 2.d0*c0*(-2.d0*JF**2/JA**3 + JF/JA**2)*dJA*Iden
	 
       Dua(1) = ddWdua(1,1)
	   Dua(2) = ddWdua(2,2)
	   Dua(3) = ddWdua(1,2)
	   
	   dWdg(1,1) = 2.d0*c2*daa(1,1)
	   dWdg(2,1) = 2.d0*c2*daa(1,2)
	   
	   ddWdgg = 2.d0*c2*Iden
	   
	  
      return
      end subroutine mat_kirigami
!************************************************************************
! Element Utilities
! modified from codes for course EN2340, Brown University 
! (c) A.F. Bower (2017)
!************************************************************************	
      subroutine abq_UEL_integrationpoints2D(n_points, n_nodes, xi, w)

      implicit none
      integer, intent(in) :: n_points
      integer, intent(in) :: n_nodes

      double precision, intent(out) :: xi(2,*)
      double precision, intent(out) :: w(*)

      integer :: i,j,k,n

      double precision :: cn,w1,w2,w11,w12,w22

    !         Defines integration points and weights for 2D continuum elements

      if ( n_points==1 ) then
        if ( n_nodes==4 .or. n_nodes==9 ) then    !     ---   4 or 9 noded quad
            xi(1, 1) = 0.D0
            xi(2, 1) = 0.D0
            w(1) = 4.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then !     ---   3 or 6 noded triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = 1.D0/3.D0
            w(1) = 1.D0/2.D0
        end if
      else if ( n_points==3 ) then
        xi(1, 1) = 0.5D0
        xi(2, 1) = 0.5D0
        w(1) = 1.D0/6.D0
        xi(1, 2) = 0.D0
        xi(2, 2) = 0.5D0
        w(2) = w(1)
        xi(1, 3) = 0.5D0
        xi(2, 3) = 0.D0
        w(3) = w(1)
      else if ( n_points==4 ) then
        if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
            !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
            !     43
            !     12
            cn = 0.5773502691896260D0
            xi(1, 1) = -cn
            xi(1, 2) = cn
            xi(1, 3) = cn
            xi(1, 4) = -cn
            xi(2, 1) = -cn
            xi(2, 2) = -cn
            xi(2, 3) = cn
            xi(2, 4) = cn
            w(1) = 1.D0
            w(2) = 1.D0
            w(3) = 1.D0
            w(4) = 1.D0
        else if ( n_nodes==3 .or. n_nodes==6 ) then
            !     xi integration points for triangle
            xi(1, 1) = 1.D0/3.D0
            xi(2, 1) = xi(1, 1)
            w(1) = -27.D0/96.D0
            xi(1, 2) = 0.6D0
            xi(2, 2) = 0.2D0
            w(2) = 25.D0/96.D0
            xi(1, 3) = 0.2D0
            xi(2, 3) = 0.6D0
            w(3) = w(2)
            xi(1, 4) = 0.2D0
            xi(2, 4) = 0.2D0
            w(4) = w(2)
        end if

      else if ( n_points==7 ) then
        ! Quintic integration for triangle
        xi(1,1) = 1.d0/3.d0
        xi(2,1) = xi(1,1)
        w(1) = 0.1125d0
        xi(1,2) = 0.0597158717d0
        xi(2,2) = 0.4701420641d0
        w(2) = 0.0661970763d0
        xi(1,3) = xi(2,2)
        xi(2,3) = xi(1,2)
        w(3) = w(2)
        xi(1,4) = xi(2,2)
        xi(2,4) = xi(2,2)
        w(4) = w(2)
        xi(1,5) = 0.7974269853d0
        xi(2,5) = 0.1012865073d0
        w(5) = 0.0629695902d0
        xi(1,6) = xi(2,5)
        xi(2,6) = xi(1,5)
        w(6) = w(5)
        xi(1,7) = xi(2,5)
        xi(2,7) = xi(2,5)
        w(7) = w(5)
      else if ( n_points==9 ) then
        !     3X3 GAUSS INTEGRATION POINTS
        !     789
        !     456
        !     123
        cn = 0.7745966692414830D0
        xi(1, 1) = -cn
        xi(1, 2) = 0.D0
        xi(1, 3) = cn
        xi(1, 4) = -cn
        xi(1, 5) = 0.D0
        xi(1, 6) = cn
        xi(1, 7) = -cn
        xi(1, 8) = 0.D0
        xi(1, 9) = cn
        xi(2, 1) = -cn
        xi(2, 2) = -cn
        xi(2, 3) = -cn
        xi(2, 4) = 0.D0
        xi(2, 5) = 0.D0
        xi(2, 6) = 0.D0
        xi(2, 7) = cn
        xi(2, 8) = cn
        xi(2, 9) = cn
        w1 = 0.5555555555555560D0
        w2 = 0.8888888888888890D0
        w11 = w1*w1
        w12 = w1*w2
        w22 = w2*w2
        w(1) = w11
        w(2) = w12
        w(3) = w11
        w(4) = w12
        w(5) = w22
        w(6) = w12
        w(7) = w11
        w(8) = w12
        w(9) = w11
      end if

      return

      end subroutine abq_UEL_integrationpoints2D


      subroutine abq_UEL_shapefunctions2D(xi,n_nodes,f,df)

      implicit none
      integer, intent(in) :: n_nodes

      double precision, intent(in) :: xi(2)
      double precision, intent(out) :: f(*)
      double precision, intent(out) :: df(9,2)
      double precision g1, g2, g3, dg1, dg2, dg3
      double precision h1, h2, h3, dh1, dh2, dh3
      double precision z,dzdp, dzdq

            if ( n_nodes==3 ) then        !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2

            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                dg1 = xi(1) - 0.5d0
                dg2 = -2.d0*xi(1)
                dg3 = xi(1) + 0.5d0
                dh1 = xi(2)-0.5d0
                dh2 = -2.d0*xi(2)
                dh3 = xi(2) + 0.5d0
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
                df(1,1) = dg1*h1
                df(1,2) = g1*dh1
                df(2,1) = dg2*h1
                df(2,2) = g2*dh1
                df(3,1) = dg3*h1
                df(3,2) = g3*dh1
                df(4,1) = dg1*h2
                df(4,2) = g1*dh2
                df(5,1) = dg2*h2
                df(5,2) = g2*dh2
                df(6,1) = dg3*h2
                df(6,2) = g3*dh2
                df(7,1) = dg1*h3
                df(7,2) = g1*dh3
                df(8,1) = dg2*h3
                df(8,2) = g2*dh3
                df(9,1) = dg3*h3
                df(9,2) = g3*dh3
            end if

      end subroutine abq_UEL_shapefunctions2D



      subroutine abq_UEL_invert2D(A,Ainv,detA)

      implicit none
      double precision, intent(in) :: A(2,2)
      double precision, intent(out) :: Ainv(2,2), detA

      double precision  :: invdetA


      detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
      invdetA = 1.d0/detA

      Ainv(1,1) =  invdetA*A(2,2)
      Ainv(1,2) = -invdetA*A(1,2)
      Ainv(2,1) = -invdetA*A(2,1)
      Ainv(2,2) =  invdetA*A(1,1)


      return
      end subroutine abq_UEL_invert2D  
