!------------------------------------
! Subsurface module
!------------------------------------

# if defined (SUBSURFACE)

MODULE SUBSURFACE_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,  &
                    Mglob,Nglob, &
                    Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,  &
                    DX,DY, DX2,DY2,DXDY, &
                    H,H0,H_u,H_v,ETA,U,V,P,Q,PQ,  &
                    MASK,MASKu,MASKv,MaxDepth, &
                    MinDepth,MASK,DT,tmp4preview, &
                    dt_over_dx,dt_over_dy
                 
  USE INPUT_Util
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY,n_west,n_east,n_suth,n_nrth
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

  LOGICAL :: SATUATION =.TRUE.
  LOGICAL :: OUT_ETA_GROUND = .FALSE.
  LOGICAL :: OUT_P_GROUND = .FALSE.
  LOGICAL :: OUT_Q_GROUND = .FALSE.
  LOGICAL :: OUT_MASK_GROUND = .FALSE.
  LOGICAL :: OUT_MASKu_GROUND = .FALSE.
  LOGICAL :: OUT_MASKv_GROUND = .FALSE.

  REAL(SP):: PLOT_INTV_ground,PLOT_COUNT_ground
  REAL(SP),DIMENSION(:,:),ALLOCATABLE::PORO_ground,ETA_ground,P_ground,Q_ground, &
                         Du_ground,Su_ground,Dv_ground,Sv_ground,  &
                         KAu_ground,KAv_ground,DepthX,DepthY, &
                         ETAu_ground,ETAv_ground,POROgrn
  INTEGER, DIMENSION(:,:),ALLOCATABLE::MASK_ground,MASKu_ground,MASKv_ground
  REAL(SP) :: S_constant,KA_constant,POR_constant

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif 

CONTAINS

! ------------------------------  
SUBROUTINE SUBSURFACE_INITIAL
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5
                    
  USE Input_Util
  IMPLICIT NONE

  CHARACTER(LEN=80) FILE_NAME
  INTEGER::LINE
  INTEGER :: ierr

! read parameters from input.txt
      FILE_NAME='input.txt'

      CALL GET_LOGICAL_VAL(SATUATION,FILE_NAME,'SATUATION',line,ierr)

      IF(ierr==1)THEN
        SATUATION = .TRUE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'No SATUATION defined, USE defalt: TRUE'
         WRITE(3,'(A80)')'No SATUATION defined, USE defalt: TRUE'
      endif
# else
         WRITE(*,'(A80)')'No SATUATION defined, USE defalt: TRUE'
         WRITE(3,'(A80)')'No SATUATION defined, USE defalt: TRUE'
# endif
       ENDIF

      CALL GET_LOGICAL_VAL(OUT_P_GROUND,FILE_NAME,'P_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_Q_GROUND,FILE_NAME,'Q_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_ETA_GROUND,FILE_NAME,'ETA_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASK_GROUND,FILE_NAME,'MASK_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASKu_GROUND,FILE_NAME,'MASKu_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASKv_GROUND,FILE_NAME,'MASKv_GROUND',line,ierr)

      CALL GET_Float_VAL(KA_constant,FILE_NAME,'Kappa',line,ierr)
      IF(ierr==1)THEN
        KA_constant = 0.01
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Kappa: Use default,', '0.01'
         WRITE(3,'(A80)')'Kappa: Use default,', '0.01'
      endif
# else
         WRITE(*,'(A80)')'Kappa: Use default,', '0.01'
         WRITE(3,'(A80)')'Kappa: Use default,', '0.01'
# endif
      ENDIF

! constant porosity
      CALL GET_Float_VAL(POR_constant,FILE_NAME,'POROSITY_GROUND',line,ierr)
      IF(ierr==1)THEN
        POR_constant = 0.3
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'POROgrn: Use default,', '0.3'
         WRITE(3,'(A80)')'POROgrn: Use default,', '0.3'
      endif
# else
         WRITE(*,'(A80)')'POROgrn: Use default,', '0.3'
         WRITE(3,'(A80)')'POROgrn: Use default,', '0.3'
# endif
      ENDIF

      CALL GET_Float_VAL(S_constant,FILE_NAME,'HardBottom',line,ierr)
      IF(ierr==1)THEN
        S_constant = MaxDepth
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'S_constant: Use default,', 'MaxDepth'
         WRITE(3,'(A80)')'S_constant: Use default,', 'MaxDepth'
      endif
# else
         WRITE(*,'(A80)')'S_constant: Use default,', 'MaxDepth'
         WRITE(3,'(A80)')'S_constant: Use default,', 'MaxDepth'
# endif
      ENDIF


! ALLOCATION

        ALLOCATE(ETA_ground(Mloc,Nloc),P_ground(Mloc1,Nloc),Q_ground(Mloc,Nloc1), &
                 Du_ground(Mloc1,Nloc),Dv_ground(Mloc,Nloc1), &
                 Su_ground(Mloc1,Nloc),Sv_ground(Mloc,Nloc1), &
                 KAu_ground(Mloc1,Nloc),KAv_ground(Mloc,Nloc1), &
                 MASK_ground(Mloc,Nloc),MASKu_ground(Mloc1,Nloc), &
                 POROgrn(Mloc,Nloc), &
                 MASKv_ground(Mloc,Nloc1), &
                 DepthX(Mloc1,Nloc),DepthY(Mloc,Nloc1), &
                 ETAu_ground(Mloc1,Nloc),ETAv_ground(Mloc,Nloc1))

! assign variables

    MASKu_ground = 0
    MASKv_ground = 0
    MASK_ground = 0
    ETA_ground = ZERO
    ETAu_ground = ZERO
    ETAv_ground = ZERO
    P_ground = ZERO
    Q_ground = ZERO
    Du_ground = ZERO
    Dv_ground = ZERO

    Su_ground = S_constant
    Sv_ground = S_constant
    KAu_ground = KA_constant 
    KAv_ground = KA_constant
    POROgrn = POR_constant

! depthx and depthy

! when first time calculate H_u and H_v, they represent depthx and depthy
     DepthX(I,J) = H_u(I,J)
     DepthY(I,J) = H_v(I,J)

    CALL UPDATE_GROUND_MASK_ETA



END SUBROUTINE SUBSURFACE_INITIAL

!------------------------

SUBROUTINE UPDATE_GROUND_MASK_ETA

! eta at u v points

   DO J=1,Nloc
   DO I=Ibeg,Iend1
       ETAu_ground(I,J)=0.5_SP*(ETA_ground(I-1,J)+ETA_ground(I,J))
   ENDDO
   ENDDO

   DO J=Jbeg,Jend1
   DO I=1,Mloc
       ETAv_ground(I,J)=0.5_SP*(ETA_ground(I,J-1)+ETA_ground(I,J))
   ENDDO
   ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH_2(ETAu_ground)
     CALL PHI_EXCH_3(ETAv_ground)
# endif

! thickness

    DO J=1,Nloc
    DO I=1,Mloc1
      IF(MASKu(I,J) == 0)THEN
        Du_ground(I,J) = Su_ground(I,J) +ETAu_ground(I,J)
      ELSE
        Du_ground(I,J) = Su_ground(I,J)-DepthX(I,J)
      ENDIF
    ENDDO
    ENDDO

    DO J=1,Nloc1
    DO I=1,Mloc
      IF(MASKv(I,J) == 0)THEN
        Dv_ground(I,J) = Sv_ground(I,J) +ETAv_ground(I,J)
      ELSE
        Dv_ground(I,J) = Sv_ground(I,J)-DepthY(I,J)
      ENDIF
    ENDDO
    ENDDO


END SUBROUTINE UPDATE_GROUND_MASK_ETA

SUBROUTINE UPDATE_SUBSURFACE
  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL

  IMPLICIT NONE
  INTEGER::ISTEP,ISTAGE,IVAR

     CALL UPDATE_GROUND_MASK_ETA

     CALL DISPERSION_ground

     CALL FLUXES_ground

     
     DO J=Jbeg,Jend
     DO I=Ibeg,Iend1
       P_ground(I,J) = - KAu_ground(I,J)*Du_ground(I,J) &
                    *(ETA_ground(I,J)-ETA_ground(I-1,J))/DX
     ENDDO
     ENDDO

     DO J=Jbeg,Jend1
     DO I=Ibeg,Iend
       Q_ground(I,J) = - KAv_ground(I,J)*Dv_ground(I,J) &
                    *(ETA_ground(I,J)-ETA_ground(I,J-1))/DY
     ENDDO
     ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH_2(P_ground)
     CALL PHI_EXCH_3(Q_ground)
# endif     


     DO J=Jbeg,Jend
     DO I=Ibeg,Iend
       IF(MASK(I,J)==0)THEN
         ETA_ground(I,J) = ETA_ground(I,J)  &
              -1.0_SP/POROgrn(I,J)*((P_ground(I+1,J)-P_ground(I,J))*dt_over_dx  &
              +(Q_ground(I,J+1)-Q_ground(I,J))*dt_over_dy)
       ELSE
         ETA_ground(I,J) = ETA(I,J)
       ENDIF
     ENDDO
     ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH_1(ETA_ground)
# endif

    CALL FILLIN_GRN_GHOST

END SUBROUTINE UPDATE_SUBSURFACE


SUBROUTINE FILLIN_GRN_GHOST

!	west
#	if defined(PARALLEL)
	if(	n_west .eq.	MPI_PROC_NULL	)then
#	endif

		DO J=1,Nloc
		DO I=1,Nghost
                        ETA_ground(I,J) = ETA_ground(Ibeg,J)
			P_ground(I,J)=P_ground(Ibeg,J)
			Du_ground(I,J)=Du_ground(Ibeg,J)
		ENDDO
		ENDDO

		DO J=Jbeg,Jend1
		DO I=1,Nghost
			Q_ground(I,J)=Q_ground(Ibeg,J)
			Dv_ground(I,J)=Dv_ground(Ibeg,J)
		ENDDO
		ENDDO

#	if defined(PARALLEL)
	endif
#	endif

!	east
#	if defined(PARALLEL)
	if(	n_east .eq.	MPI_PROC_NULL	)then
#	endif

		DO J=1,Nloc
		DO I=Iend1+1,Mloc1
                        ETA_ground(I-1,J) = ETA_ground(Iend,J)
			P_ground(I,J)=P_ground(Iend1,J)
			Du_ground(I,J)=Du_ground(Iend1,J)
		ENDDO
		ENDDO

		DO J=Jbeg,Jend1
		DO I=Iend1,Mloc
			Q_ground(I,J)=Q_ground(Iend,J)
			Dv_ground(I,J)=Dv_ground(Iend,J)
		ENDDO
		ENDDO

#	if defined(PARALLEL)
	endif
#	endif

!	south
#	if defined(PARALLEL)
	if(	n_suth .eq.	MPI_PROC_NULL	)then
#	endif

		DO J=1,Nghost
		DO I=1,Mloc
                        ETA_ground(I,J) = ETA_ground(I,Jbeg)
			Q_ground(I,J)=Q_ground(I,Jbeg)
			Dv_ground(I,J)=Dv_ground(I,Jbeg)
		ENDDO
		ENDDO

		DO J=1,Nghost
		DO I=Ibeg,Iend1
			P_ground(I,J)=P_ground(I,Jbeg)
			Du_ground(I,J)=Du_ground(I,Jbeg)
		ENDDO
		ENDDO

#	if defined(PARALLEL)
	endif
#	endif

!	north
#	if defined(PARALLEL)
	if(	n_nrth .eq.	MPI_PROC_NULL	)then
#	endif

		DO J=Jend1+1,Nloc1
		DO I=1,Mloc
                        ETA_ground(I,J-1) = ETA_ground(I,Jend)
			Q_ground(I,J)=Q_ground(I,Jend1)
			Dv_ground(I,J)=Dv_ground(I,Jend1)
		ENDDO
		ENDDO

		DO J=Jend1,Nloc
		DO I=Ibeg,Iend1
			P_ground(I,J)=P_ground(I,Jend)
			Du_ground(I,J)=Du_ground(I,Jend)
		ENDDO
		ENDDO

#	if defined(PARALLEL)
	endif
#	endif


END SUBROUTINE FILLIN_GRN_GHOST

SUBROUTINE DISPERSION_ground
  USE GLOBAL
  IMPLICIT NONE

     REAL(SP),Dimension(Mloc,Nloc) :: DU,DV,DUt,DVt
     REAL(SP) :: UxxVxy,UxyVyy,HUxxHVxy,HUxyHVyy, &
                 UxxVxy_x,UxxVxy_y,UxyVyy_x,UxyVyy_y, &
                 HUxxHVxy_x,HUxxHVxy_y,HUxyHVyy_x,HUxyHVyy_y, &
                 rh,rhx,rhy,reta,ken1,ken2,ken3,ken4,ken5
# if defined (CARTESIAN)
     REAL(SP) :: omega_0,omega_1
     REAL(SP),Dimension(Mloc,Nloc) :: omega            
# endif
    
! uxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,U,Uxx)
! uxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,U,Uxy)
! vxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,V,Vxy)
! vyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,V,Vyy)

# if defined (CARTESIAN)
! gamma2.ne.0
    IF(Gamma2>ZERO)THEN
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,U,Ux)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,V,Vx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,U,Uy)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,V,Vy)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Eta,ETAx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Eta,ETAy)
    ELSEIF(SHOW_BREAKING)THEN
# else
    IF(SHOW_BREAKING)THEN
# endif
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Eta,ETAx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Eta,ETAy)
    ENDIF

! DU DV
     DO J=1,Nloc-1
     DO I=1,Mloc-1
       DU(I,J)=Max(Depth(I,J),MinDepthFrc)*U(I,J)
       DV(I,J)=Max(Depth(I,J),MinDepthFrc)*V(I,J)


# if defined (CARTESIAN)
       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX-(Q(I,J+1)-Q(I,J))/DY
# else
       ETAT(I,J)=-(P(I+1,J)-P(I,J))/DX(I,J)-(Q(I,J+1)-Q(I,J))/DY(I,J)
# endif
     ENDDO
     ENDDO 

! ETAT
# if defined (CARTESIAN)
    IF(Gamma2>ZERO)THEN

       DO J=1,Nloc
       DO I=1,Mloc

	  Ut(I,J) = (U(I,J)-U0(I,J)) / DT   !ykchoi
	  Vt(I,J) = (V(I,J)-V0(I,J)) / DT   !ykchoi

        DUt(I,J)=Max(Depth(I,J),MinDepthFrc)*Ut(I,J)
        DVt(I,J)=Max(Depth(I,J),MinDepthFrc)*Vt(I,J)
       ENDDO
       ENDDO

    ELSEIF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN
# else
    IF(SHOW_BREAKING .OR. WAVEMAKER_VIS)THEN
# endif

    ENDIF

! DUxx
    CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DU,DUxx)
! DUxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DU,DUxy)
! DVxy
    CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DV,DVxy)
! DVyy
    CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DV,DVyy)
      
# if defined (CARTESIAN)
    IF(Gamma2>ZERO)THEN
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DU,DUx)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DV,DVx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DU,DUy)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DV,DVy)
     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Ut,Utx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Vt,Vty)

     CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,Ut,Utxx)
     CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,Vt,Vtyy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,Ut,Utxy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,Vt,Vtxy)

     CALL DERIVATIVE_X(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DUt,DUtx)
     CALL DERIVATIVE_Y(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DVt,DVty)

     CALL DERIVATIVE_XX(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DUt,DUtxx)
     CALL DERIVATIVE_YY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DY,DVt,DVtyy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DUt,DUtxy)
     CALL DERIVATIVE_XY(Mloc,Nloc,Ibeg,Iend,Jbeg,Jend,MASK9,DX,DY,DVt,DVtxy)
    ENDIF
# endif
      
!  boundary conditions
# if defined (PARALLEL)
    if(n_west.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_WEST)THEN
# endif
    IF (WaveMaker(1:11)=='LEFT_BC_IRR')THEN
     ! do nothing
    ELSE
     DO J=1,Nloc
       Uxy(Ibeg,J)=ZERO
       DUxy(Ibeg,J)=ZERO
       Vxy(Ibeg,J)=ZERO
       DVxy(Ibeg,J)=ZERO
       Utxy(Ibeg,J)=ZERO
       DUtxy(Ibeg,J)=ZERO
       Vtxy(Ibeg,J)=ZERO
       DVtxy(Ibeg,J)=ZERO
     ENDDO
    ENDIF ! left_bc wavemaker
# if defined(COUPLING)  
   ENDIF 
# endif
  
# if defined (PARALLEL)
    endif  
# endif

# if defined (PARALLEL)
    if(n_east.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_EAST)THEN
# endif
     DO J=1,Nloc
       Uxy(Iend,J)=ZERO
       DUxy(Iend,J)=ZERO
       Vxy(Iend,J)=ZERO
       DVxy(Iend,J)=ZERO
       Utxy(Iend,J)=ZERO
       DUtxy(Iend,J)=ZERO
       Vtxy(Iend,J)=ZERO
       DVtxy(Iend,J)=ZERO
     ENDDO 
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif
  
# if defined (PARALLEL)
    if(n_suth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_SOUTH)THEN
# endif
     DO I=1,Mloc
       Uxy(I,Jbeg)=ZERO
       DUxy(I,Jbeg)=ZERO
       Vxy(I,Jbeg)=ZERO
       DVxy(I,Jbeg)=ZERO
       Utxy(I,Jbeg)=ZERO
       DUtxy(I,Jbeg)=ZERO
       Vtxy(I,Jbeg)=ZERO
       DVtxy(I,Jbeg)=ZERO
     ENDDO   
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif

# if defined (PARALLEL)
    if(n_nrth.eq.MPI_PROC_NULL) then
# endif
# if defined(COUPLING)
    IF(.NOT.IN_DOMAIN_NORTH)THEN
# endif
     DO I=1,Mloc
       Uxy(I,Jend)=ZERO
       DUxy(I,Jend)=ZERO
       Vxy(I,Jend)=ZERO
       DVxy(I,Jend)=ZERO
       Utxy(I,Jend)=ZERO
       DUtxy(I,Jend)=ZERO
       Vtxy(I,Jend)=ZERO
       DVtxy(I,Jend)=ZERO
     ENDDO 
# if defined(COUPLING)
    ENDIF
# endif
# if defined (PARALLEL)
    endif  
# endif
     
    CALL EXCHANGE_DISPERSION
     
! calculate V1p  without nonlinear dispersion
     DO J=1,Nloc
     DO I=1,Mloc
# if defined (CARTESIAN)
       U4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       V4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*DEPTH(I,J)*DEPTH(I,J)*(Uxy(I,J)+Vyy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))      

       IF(gamma2>ZERO)THEN
	   UxxVxy = Uxx(I,J) + Vxy(I,J)
	   UxyVyy = Uxy(I,J) + Vyy(I,J)

	   HUxxHVxy = DUxx(I,J) + DVxy(I,J)
	   HUxyHVyy = DUxy(I,J) + DVyy(I,J)

	   rh = Depth(I,J)
	   reta = Eta(I,J)

  	   ken1 = ( 1.0_SP/6.0_SP - Beta_1 + Beta_1*Beta_1 )*rh*reta*Beta_2   &
	         + ( 1.0_SP/2.0_SP*Beta_1*Beta_1 - 1.0_SP/6.0_SP )*reta*reta*Beta_2*Beta_2
	   ken2 = ( Beta_1 - 1.0_SP/2.0_SP )*reta*Beta_2

	   U4(I,J) = U4(I,J) + gamma2*MASK9(I,J)*( ken1*UxxVxy + ken2*HUxxHVxy )
	   V4(I,J) = V4(I,J) + gamma2*MASK9(I,J)*( ken1*UxyVyy + ken2*HUxyHVyy )
	 ENDIF

# endif

# if defined (ZALPHA)
       U4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       V4(I,J)=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)  &
               *DEPTH(I,J)*DEPTH(I,J)*(Uxy(I,J)+Vyy(I,J)) &
                +(Beta_1-1.0_SP/2.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# endif

               
# if defined (CARTESIAN)
       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))
       
	 !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# else
       !U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
	 U1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)*(Uxx(I,J)+Vxy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxx(I,J)+DVxy(I,J))

       !V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)  & !ykchoi
       V1p(I,J)=0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)  &
                *(Uxy(I,J)+Vyy(I,J)) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(DUxy(I,J)+DVyy(I,J))
# if defined(SPH_EXTRA_DISP)
 
       U1p(I,J)=U1p(I,J)+0.5_SP*(1.0_SP-Beta_1)  &
                *DEPTH(I,J)*DEPTH(I,J)*(-1.0_SP/R_earth*  &
                TAN(Lat_theta(I,J))*0.5_SP*(VxL(I,J)+VxR(I,J))) &
               +(Beta_1-1.0_SP)*DEPTH(I,J)*(-1.0_SP/R_earth* &
                TAN(Lat_theta(I,J))*0.5_SP*(HVxL(I,J)+HVxR(I,J)))


       V1p(I,J)=V1p(I,J)+0.5_SP*(1.0_SP-Beta_1) &
                 *DEPTH(I,J)*DEPTH(I,J) &
                /R_earth*( TAN(Lat_theta(I,J))*0.5_SP*(UxL(I,J)+UxR(I,J)) &
                          -TAN(Lat_theta(I,J))*0.5_SP*(VyL(I,J)+VyR(I,J)) &
                          -V(I,J)/R_earth/COS(Lat_theta(I,J))/COS(Lat_theta(I,J)) ) &
               +(Beta_1-1.0_SP)*DEPTH(I,J) &
                /R_earth*( 0.5_SP*(HUxL(I,J)+HUxR(I,J))*TAN(Lat_theta(I,J)) &
                          -0.5_SP*(HVyL(I,J)+HVyR(I,J))*TAN(Lat_theta(I,J)) &
                          -HV(I,J)/R_earth/COS(Lat_theta(I,J))/COS(Lat_theta(I,J)) &        
                 )
# endif

# endif

     ENDDO
     ENDDO

# if defined (CARTESIAN)
     IF(gamma2>ZERO)THEN
       DO J=Jbeg,Jend
       DO I=Ibeg,Iend
! 
        UxxVxy=Uxx(I,J)+Vxy(I,J)
        UxyVyy=Uxy(I,J)+Vyy(I,J)
        UxxVxy_x=(Uxx(I+1,J)+Vxy(I+1,J)-Uxx(I-1,J)-Vxy(I-1,J))/2.0_SP/DX
        UxxVxy_y=(Uxx(I,J+1)+Vxy(I,J+1)-Uxx(I,J-1)-Vxy(I,J-1))/2.0_SP/DY
        UxyVyy_x=(Uxy(I+1,J)+Vyy(I+1,J)-Uxy(I-1,J)-Vyy(I-1,J))/2.0_SP/DX
        UxyVyy_y=(Uxy(I,J+1)+Vyy(I,J+1)-Uxy(I,J-1)-Vyy(I,J-1))/2.0_SP/DY

        HUxxHVxy=DUxx(I,J)+DVxy(I,J)
        HUxyHVyy=DUxy(I,J)+DVyy(I,J)
        HUxxHVxy_x=(DUxx(I+1,J)+DVxy(I+1,J)-DUxx(I-1,J)-DVxy(I-1,J))/2.0_SP/DX
        HUxxHVxy_y=(DUxx(I,J+1)+DVxy(I,J+1)-DUxx(I,J-1)-DVxy(I,J-1))/2.0_SP/DY
        HUxyHVyy_x=(DUxy(I+1,J)+DVyy(I+1,J)-DUxy(I-1,J)-DVyy(I-1,J))/2.0_SP/DX
        HUxyHVyy_y=(DUxy(I,J+1)+DVyy(I,J+1)-DUxy(I,J-1)-DVyy(I,J-1))/2.0_SP/DY

        rh=Depth(I,J)
        rhx=(Depth(I+1,J)-Depth(I-1,J))/2.0_SP/DX
        rhy=(Depth(I,J+1)-Depth(I,J-1))/2.0_SP/DY
        reta=Eta(I,J)

        U1pp(I,J)=-reta*Beta_2*ETAx(I,J)*Beta_2*(Utx(I,J)+Vty(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxx(I,J)+Vtxy(I,J))&
                  -ETAx(I,J)*Beta_2*(DUtx(I,J)+DVty(I,J)) -reta*Beta_2*(DUtxx(I,J)+DVtxy(I,J))

        V1pp(I,J)=-reta*Beta_2*ETAy(I,J)*Beta_2*(Utx(I,J)+Vty(I,J)) - 0.5_SP*reta*reta*Beta_2*Beta_2*(Utxy(I,J)+Vtyy(I,J))&
                  -ETAy(I,J)*Beta_2*(DUtx(I,J)+DVty(I,J)) -reta*Beta_2*(DUtxy(I,J)+DVtyy(I,J))
        
  
	  ken1 = Beta_1*( 1.0_SP - Beta_1 )*rh*ETAT(I,J)*Beta_2 - Beta_1*Beta_1*reta*Beta_2*ETAT(I,J)*Beta_2
	  ken2 = Beta_1*( 1.0_SP - Beta_1 )*rh*reta*Beta_2 - 0.5*Beta_1*Beta_1*reta*reta*Beta_2*Beta_2
	  ken3 = Beta_1*ETAT(I,J)*Beta_2
	  ken4 = Beta_1*reta*Beta_2

	  U1pp(I,J) = U1pp(I,J) - ken1*UxxVxy - ken2*( Utxx(I,J)+Vtxy(I,J) ) + ken3*HUxxHVxy + ken4*(DUtxx(I,J)+DVtxy(I,J))
	  V1pp(I,J) = V1pp(I,J) - ken1*UxyVyy - ken2*( Utxy(I,J)+Vtyy(I,J) ) + ken3*HUxyHVyy + ken4*(DUtxy(I,J)+DVtyy(I,J))

         ken1=(Beta_1-1.0_SP)*(rhx+ETAx(I,J))*Beta_2
         ken2=(Beta_1-1.0_SP)*(rh+reta)*Beta_2
         ken3=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhx*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhx*reta*Beta_2+rh*ETAx(I,J)*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*ETAx(I,J)*Beta_2*Beta_2 )
         ken4=( 0.5_SP*(1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rh*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*rh*reta*Beta_2 &
                      +0.5_SP*(Beta_1*Beta_1-1.0_SP)*reta*reta*Beta_2*Beta_2 )
         ken5=( (1.0_SP-Beta_1)*(1.0_SP-Beta_1)*rh*rhy*Beta_2*Beta_2-Beta_1*(1.0_SP-Beta_1)*(rhy*reta*Beta_2+rh*ETAy(I,J)*Beta_2) &
                    +(Beta_1*Beta_1-1.0_SP)*reta*ETAy(I,J)*Beta_2*Beta_2 )

        U2(I,J)=ken1*(U(I,J)*HUxxHVxy+V(I,J)*HUxyHVyy) &
                +ken2*(Ux(I,J)*HUxxHVxy+U(I,J)*HUxxHVxy_x &
                    +Vx(I,J)*HUxyHVyy+V(I,J)*HUxyHVyy_x) &
                +ken3 & 
                   *(U(I,J)*UxxVxy+V(I,J)*UxyVyy) &
                +ken4  &
                   *(Ux(I,J)*UxxVxy+U(I,J)*UxxVxy_x+Vx(I,J)*UxyVyy+V(I,J)*UxyVyy_x) &
                +Beta_2*Beta_2*(DUx(I,J)+DVy(I,J)+reta*Beta_2*(Ux(I,J)+Vy(I,J)))  &
                   *(HUxxHVxy+ETAx(I,J)*Beta_2*(Ux(I,J)+Vy(I,J))+reta*Beta_2*UxxVxy)

        ken1=(Beta_1-1.0_SP)*(rhy+ETAy(I,J))*Beta_2

        V2(I,J)=ken1*(U(I,J)*HUxxHVxy+V(I,J)*HUxyHVyy) &
                +ken2*(Uy(I,J)*HUxxHVxy+U(I,J)*HUxxHVxy_y &
                    +Vy(I,J)*HUxyHVyy+V(I,J)*HUxyHVyy_y) &
                +ken5 & 
                   *(U(I,J)*UxxVxy+V(I,J)*UxyVyy) &
                +ken4  &
                   *(Uy(I,J)*UxxVxy+U(I,J)*UxxVxy_y+Vy(I,J)*UxyVyy+V(I,J)*UxyVyy_y) &
                +Beta_2*Beta_2*(DUx(I,J)+DVy(I,J)+reta*Beta_2*(Ux(I,J)+Vy(I,J)))  &
                   *(HUxyHVyy+ETAy(I,J)*Beta_2*(Ux(I,J)+Vy(I,J))+reta*Beta_2*UxyVyy)

        omega_0=Vx(I,J)-Uy(I,J)

        omega_1=( b2*rhx + Beta_1*ETAx(I,J) )*Beta_2*( HUxyHVyy + (b2*rh+Beta_1*reta)*Beta_2*UxyVyy )  &
              - ( b2*rhy + Beta_1*ETAy(I,J) )*Beta_2*( HUxxHVxy + (b2*rh+Beta_1*reta)*Beta_2*UxxVxy )

	omega(I,J)=omega_0+omega_1

        IF(OUT_VORmax) THEN
        IF(abs(omega(I,J)).GT.VorticityMax(I,J)) THEN
        VorticityMax(I,J)=omega(I,J)
        ENDIF
        ENDIF

       ken1=((Beta_1-1.0_SP/2.0_SP)*(reta+rh)*Beta_2)
       ken2=(1.0_SP/3.0_SP-Beta_1+0.5_SP*Beta_1*Beta_1)*rh*rh*Beta_2*Beta_2  &
               + (1.0_SP/6.0_SP-Beta_1+Beta_1*Beta_1)*rh*reta*Beta_2 &
               +(1.0_SP/2.0_SP*Beta_1*Beta_1-1.0_SP/6.0_SP)*reta*reta*Beta_2*Beta_2

       U3(I,J)=-V(I,J)*omega_1 - omega_0 &
                 *(ken1*HUxyHVyy &
                   +ken2*UxyVyy)

	 V3(I,J) = U(I,J)*omega_1 + omega_0 &
	            *(ken1*HUxxHVxy &
                    +ken2*UxxVxy)
       ENDDO
       ENDDO            

     ENDIF  
# endif


END SUBROUTINE DISPERSION_ground
     USE GLOBAL
     IMPLICIT NONE

      IF(HIGH_ORDER(1:3)=='FOU') THEN    

       CALL CONSTRUCTION_HO   

       CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,UxL,UxR,VyL,VyR,HxL,HxR,HyL,HyR, &
            SxL,SxR,SyL,SyR)

      ELSE                                     
       CALL DelxyFun
       CALL CONSTRUCTION 
       CALL WAVE_SPEED(Mloc,Nloc,Mloc1,Nloc1,UxL,UxR,VyL,VyR,HxL,HxR,HyL,HyR, &
            SxL,SxR,SyL,SyR)

     ENDIF

     IF(CONSTR(1:3)=='HLL')THEN
       CALL FLUX_AT_INTERFACE_HLL
     ELSE
       CALL FLUX_AT_INTERFACE
     ENDIF

!     CALL BOUNDARY_CONDITION ! will add later


SUBROUTINE FLUXES_ground

END SUBROUTINE FLUXES_ground

END MODULE SUBSURFACE_MODULE
# endif
! end the module

