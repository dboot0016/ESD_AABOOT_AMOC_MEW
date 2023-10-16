!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Cimatoribus coupled to CC model v1.0
!   Start: 26-10-2022
!   Author: Amber Boot
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION CO2, Ct, Cts, Cn, Cs, Cd, At, Ats, An, As, Ad, Pt, Pts, Pn, Ps, Pd, TC, TA, TP, dCtot
      DOUBLE PRECISION Ht, Hts, Hn, Hs, Hd, CO3t, CO3ts, CO3n, CO3s, CO3d, pCO2t, pCO2ts, pCO2n, pCO2s
      DOUBLE PRECISION K0t, K1t,K2t, K0ts, K1ts, K2ts, K0n, K1n,K2n, K0s, K1s, K2s, K0d, K1d, K2d, BTt, BTts, BTn, BTs, BTd
      DOUBLE PRECISION Kt0, Kts0, Kn0, Ks0, Kd0, K11t, K11ts, K11n, K11s, K11d, Kt1, Kts1, Kn1, Ks1, Kd1
      DOUBLE PRECISION dVCt, dVCts, dVCn, dVCs, dVCd, dKCt, dKCts, dKCn, dKCs, dKCd, Kspt, Kspts, Kspn
      DOUBLE PRECISION Ksps, Kspd, Kspds, Ksprest, Ksprests, Kspresn, Kspress, Kspresd
      DOUBLE PRECISION T0, Tt, Tts, Tn, Ts, Td, S0, St, Sts, Sn, Ss, Sd, rho0, rhot, rhots, rhon, rhos, rhod, alpha, beta
      DOUBLE PRECISION V0, Vt, Vts, Vn, Vs, Vd, VAt, SAt, SAts, SAn, SAs, dt, dts, dn, ds, D
      DOUBLE PRECISION d0, LxA, LxS , Ly, A, Ar, t_dim, S_dim, D_dim
      DOUBLE PRECISION Cphys_t, Cphys_ts, Cphys_n, Cphys_s, Aphys_t, Aphys_ts, Aphys_n, Aphys_s, Pphys_t, Pphys_ts, Pphys_n, Pphys_s
      DOUBLE PRECISION Ccarb_t, Ccarb_ts, Ccarb_n, Ccarb_s, Ccarb_d, Acarb_t, Acarb_ts, Acarb_n, Acarb_s
      DOUBLE PRECISION Ccarb_t1, Ccarb_ts1, Ccarb_n1, Ccarb_s1, Ccarb_d1, Ccarb_t2, Ccarb_ts2, Ccarb_n2, Ccarb_s2, Ccarb_d2
      DOUBLE PRECISION Cbio_t, Cbio_ts, Cbio_n, Cbio_s, Pbio_t, Pbio_ts, Pbio_n, Pbio_s
      DOUBLE PRECISION Cair_t, Cair_ts, Cair_n, Cair_s, Criver_t, Ariver_t, Priver_t
      DOUBLE PRECISION qe, qs, qu, qn, qEk, QN1, QS1, QS2, rn, rs, eta, kappa, Agm, tau, Ea, Es, fs
      DOUBLE PRECISION Zt, Zts, Zn, Zs, b, rcp, FCAt, FCAts, FCAn, FCAs, Cat, Cats, Can, Cas, Cad
      DOUBLE PRECISION OmegaC_t, OmegaC_ts, OmegaC_n, OmegaC_s, OmegaC_d, OmegaC_ds
      DOUBLE PRECISION PVt, PVts, PVn, PVs, A1_C, A2_C, A3_C, B1_C, B2_C, B3_C, n, kCa, PerC, DC, Priver, Csi, Vsi, Vca
      DOUBLE PRECISION RAt, RAts, RAn, RAs, RAd, DMt, DMts, DMn, DMs, DMd, Kd
      DOUBLE PRECISION TC1, TA1, TP1, TS1, b1, b2, b3, b4, b5, b6, b7, Kt2, Kts2, Kn2, Ks2, Kd2
      DOUBLE PRECISION Srs,Srd,gamma1,psi1, Ep,drs,dtot,Vrs,Vrd,SArs,SArd,Ep_to_s,rp

        ! Switches

        ! Start with defining the state variables
        St = U(1)         ! Salinity thermocline box
        Sts = U(2)         ! Salinity ts box
        Sn = U(3)         ! Salinity northern box
        Ss = U(4)         ! Salinity southern box
        D = U(5)         ! Thermocline depth
        Srs = U(6)
        Srd = U(7)
        
        t_dim = 100*365*86400
        D_dim = 1000
        S_dim = 35

        ! Parameters physical model (not from Cimatoribus: dn, ds, Tt, Ts)
        V0 = 3d17       ! Total volume of the basin [m3]
        Vn = 3d15       ! Volume northern box [m3]
        Vs = 9d15       ! Volume southern box [m3]
        A = 1d14        ! Horizontal area of Atlantic pycnocline [m2]
        LxA = 1d7       ! Zonal exten of the Atlantic at its southern end [m]
        LxS = 3d7       ! Zonal extent of the Southern Ocean [m]
        Ly = 1d6        ! Meridional extent of frontal region Southern Ocean [m]
        Ar = 4d14       ! Area of thermocline + ts [m2]
        Vts = LxA*Ly*D*D_dim/2          ! Volume ts box [m3]
        Vt = A*D*D_dim           ! Volume thermocline box

        dt = D*D_dim  ! Depth thermocline box
        dts = D*D_dim/2   ! Depth ts box
        dn = 300   ! Depth northern box [m] (first guess for now)
        ds = 300   ! Depth southern box [m] (first guess for now)
        SAn = Vn/dn      ! Surface area northern box [m2]
        SAs = Vs/ds      ! Surfacer area southern box [m2]
        SAt = Vt/dt      ! Surdace area thermocline box [m2]
        SAts = LxA*Ly   ! Surface area ts box [m2]

        Vd = V0-Vn-Vs-Vts-Vt           ! Volume deep box [m3]
   
        tau = 0.1       ! Average zonal wind stress amplitude [Nm-2]
        Agm = 1700      ! Eddy diffusivity [m2s-1]
        kappa = 1d-5    ! Vertical diffusivity [m2s-1]
        eta = 3d4       ! Hydraulic constant [ms-1]
        fs = -1d-4      ! Coriolis paramter [s-1]
        Es = PAR(2)*1d6     ! Symmetric freshwater flux [m3s-1]
        rn = 5d6        ! Transport by northern subtropical gyre [m3s-1]
        rs = 1d7        ! Transport by southern subtropical gyre
        Ea = par(1)     ! Asymmetric freshater flux [m3s-1]

        rho0 = 1027.5   ! Reference density [kgm-3]
        S0 = 35         ! Reference salinity [psu]
        T0 = 5          ! Reference temperature [degree C]
        Tn = 5          ! Temperature northern box [degree C]
        Tts = 10         ! Temperature ts box [degree C]
        alpha = 2d-4    ! Thermal expansion coefficient [K-1]
        beta = 8d-4     ! Haline contraction coefficient [psu-1]

        Tt = 23.44       ! Temperature thermocline box [degree C] (SCP-M, box 1)
        Ts = 0.93   ! Temperature southern box [degree C] (SCP-M, box 5)
        Td =  1.8  ! Temperature deep box [degree C] (SCP-M, box 6)
        
        dtot = 4000
        drs = 300
        Vrs = drs/dtot*4*V0
        Vrd = 4*V0-Vrs
        SArs = Vrs/drs
        SArd = SArs
        
        Ep = 0
        Ep_to_s = 0.855d6
        psi1 = 18d6
        gamma1 = 30d6
        rp = 90d6

        ! Determine ocean circulation
        qEk = tau*Lxs/(rho0*ABS(fs))
        qe = Agm*D*D_dim*LxA/Ly
        qs = qEk - qe
        qu = kappa*A/(D*D_dim)
        qn = eta*(D*D_dim*D*D_dim)*(alpha*(Tts-Tn)+beta*(Sn*S_dim-Sts*S_dim))

        ! Conservation laws for S, C, A, and P
        Sd= ((S0*V0*5.0-(St*Vt+Sts*Vts+Sn*Vn+Ss*Vs+Srs*Vrs+Srd*Vrd)*S_dim)/Vd)/S_dim !-(b1 + b2*Sn + b5*D*St + b6*D*Sts + b7*Ss)/(b3+b4*D)

        ! Heaviside functions
        IF (qn.lt.0) THEN
            QN1 = 0
        ELSE
            QN1 = qn
        ENDIF

        IF (qs.gt.0) THEN
            QS2 = 0
            QS1 = qs
        ELSE
            QS2 = qs
            QS1 = 0
        ENDIF

        ! State the ODEs
        F(1) = ((QS1*Sts+QS2*St)+qu*Sd-QN1*St+rs*(Sts-St)+rn*(Sn-St)+2*Es+Ep)*t_dim/Vt-(St*(qu+qs-QN1)*t_dim/(Ar*D_dim))/D                                     ! Salinity thermocline box
        F(2) = (qEk*Ss-qe*Sts-(QS1*Sts+QS2*St)+rs*(St-Sts))*t_dim/Vts-(Sts*(qu+qs-QN1)*t_dim/(Ar*D_dim))/D                                            ! Salinity ts box
        F(3) = ((QN1+rn)*(St-Sn)-(Es+Ea*1e6))*t_dim/Vn            ! Salinity northern box
        F(4) = ((QS1*Sd+QS2*Ss)+qe*Sts-qEk*Ss-(Es-Ea*1e6+Ep_to_s)+(rp+psi1)*(Srs-Ss))*t_dim/Vs   ! Salinity southern box
        F(5) = (qu+qs-QN1)*t_dim/(Ar*D_dim)                  ! Thermocline depth
		F(6) = ((gamma1+psi1)*(Srd-Srs)/Vrs+(-Ep+Ep_to_s+rp*(Ss-Srs))/Vrs)*t_dim
        F(7) = (gamma1*(Srs-Srd)/Vrd+psi1*1.0*(Sd-Srd)/Vrd)*t_dim
        
      END SUBROUTINE FUNC
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      INTEGER I 
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      

! Parameters
        PAR(1) = 0.
        PAR(2) = 0.6
 
! Special solution from transient code 
        U(1) = 1.0058376418610575    ! Salinity thermocline
        U(2) = 0.994840665484994    ! Salinity ts
        U(3) = 0.9984467757381668     ! Salinity northern
        U(4) = 0.9910742010761921     ! Salinity southern
        U(5) = 0.7382511281434748     ! Thermocline depth
        U(6) = 1
        U(7) = 1
        
        !U(1) = 9.96811E-01    ! Salinity thermocline
        !U(2) = 1.00423E+00    ! Salinity ts
        !U(3) = 9.12811E-01     ! Salinity northern
        !U(4) = 1.00692E+00     ! Salinity southern
        !U(5) = 1.75107E+00     ! Thermocline depth

      END SUBROUTINE STPNT
!----------------------------------------------------------------------

      SUBROUTINE BCND
      END SUBROUTINE BCND

    SUBROUTINE ICND(NDIM,PAR,ICP,NINT,U,UOLD,UDOT,UPOLD,FI,IJAC,DINT)
      !---------- ----

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NDIM, ICP(*), NINT, IJAC
    DOUBLE PRECISION, INTENT(IN) :: PAR(*)
    DOUBLE PRECISION, INTENT(IN) :: U(NDIM), UOLD(NDIM), UDOT(NDIM), UPOLD(NDIM)
    DOUBLE PRECISION, INTENT(OUT) :: FI(NINT)
    DOUBLE PRECISION, INTENT(INOUT) :: DINT(NINT,*)
    

    END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
