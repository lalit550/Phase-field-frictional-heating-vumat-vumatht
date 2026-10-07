
c ======== SHARED MODULE (put in a separate file if required: pf_shared.f) ========
             module pf_shared
      implicit none

      integer, parameter :: MAX_ELEM = 5000000

      real*8 H_store(MAX_ELEM)
      real*8 phi_store(MAX_ELEM)
      real*8 Tgrad(MAX_ELEM,3)
	  
      logical initialized(MAX_ELEM)

      data initialized /MAX_ELEM*.false./

      end module pf_shared
c ======================================================================
c =============================== VUMAT =================================
c ======================================================================
      subroutine vumat(
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, jInfoArray,
     2  stepTime, totalTime, dtArray, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
     7  stressNew, stateNew, enerInternNew, enerInelasNew )

      use pf_shared
      include 'vaba_param.inc'

      parameter (i_info_ElemNumStartLoc = 6)
     
      pointer (ptrjElemNum, jElemNum)
      dimension jElemNum(nblock)

      integer nblock, ndir, nshr, nstatev, nfieldv, nprops
      character*80 partName(nblock)
      integer elementNumber(nblock)
      integer intnum, jtyp, locnum, jrcd
      
      integer jInfoArray(nblock) ! Correctly define this as an array



c      real*8 stepTime, totalTime, dtArray(*)
      character*80 cmname
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),dtArray(*),
     4  stretchOld(*), defgradOld(nblock,ndir+nshr+1), fieldOld(nblock,nfieldv),
     5  stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(nblock,ndir+nshr+1), fieldNew(nblock,nfieldv),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     &  enerInternNew(nblock), enerInelasNew(nblock)
      dimension jels(nblock)
c --- local variables ---
c --- Local scalars ---
      
      real*8   lambda, mu
      real*8   H_old, H_new
      real*8   dt, rho

c --- Strain ---
      real*8    eps12_eng
      real*8   eps12_t                          ! tensor shear = eng/2

c --- Eigenvalue decomposition (2D + out-of-plane) ---
      real*8   tr_eps                           ! trace = eps11+eps22+eps33
      real*8   half_sum, half_diff, rad
      real*8   e1, e2, e3                       ! principal strains
      real*8   e1p, e2p, e3p                    ! positive parts
      real*8   e1n, e2n, e3n                    ! negative parts
      
     

c --- Stress in principal frame ---
      real*8   sig1p, sig2p, sig3p             ! tensile contrib (degraded)
      real*8   sig1n, sig2n, sig3n             ! compressive (NOT degraded)
      real*8   sig1, sig2, sig3                ! total principal

c --- Stress in global frame ---
      real*8   s11, s22, s33, s12

c --- Energy ---
      real*8   psi_plus                        ! tensile strain energy density
      

c --- Numerical tolerances ---
      real*8   ZERO_TOL, EPS_TOL
      parameter (ZERO_TOL = 1.d-20)
      parameter (EPS_TOL  = 1.d-14)
      integer i,j,constitutive_type,kGlobal
      integer elemStart
      real*8 E,xnu,xl,Gc
      real*8 twomu,eg,sixmu,alamda,lambda_lame,mu_lame
      real*8 stran(nblock,ndir+nshr)
      real*8 eps11,eps22,eps33,eps12
      real*8 trace2,diff,eig1,eig2,eig3
      real*8 pos1,pos2,pos3,neg1,neg2,neg3,tr_pos,tr_neg
      real*8 sig1_tens,sig2_tens,sig3_tens
      real*8 sig1_comp,sig2_comp,sig3_comp
      real*8 sig1_tot,sig2_tot,sig3_tot
      real*8 theta_p,cth,sth
      real*8 s11_tot,s22_tot,s12_tot,s33_tot
      real*8 eigVal(3),trp1,trp2,psip,psit,H
      real*8 stressPower,dtloc
      real*8 phi,g
      real*8 x,y
      real*8 s1,s2,s3,s4
      real*8 pi,base
      real*8 del_status, J_est, eps1_lim, J_lim, gamma_lim
      character*80 cpname
 
      parameter (eps1_lim = 0.35d0)   ! max principal strain limit
      parameter (J_lim    = 0.05d0)   ! Jacobian collapse limit
      parameter (gamma_lim= 0.50d0)   ! shear strain limit
      parameter (pi=3.141592653589d0)
      real*8 L_CAP
      parameter (L_CAP = 1.0d8)
      real*8 L11, L12, L21, L22
	  
	  
c ====== FRICTION ADDITIONS: local variables ======
      real*8 mu_fric, Cv_mech
      real*8 grad1,grad2,grad3,gradnorm2,gradnorm
      real*8 n1,n2,n3,t1,t2,vs,vt,v1n,v2n
      real*8 F_old(2,2),F_new(2,2),Fdot(2,2),Finv(2,2),L(2,2),detF
      real*8 tn,tn_comp,gamma_pf,qdot_fric,deltaT_fric
	  

      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum = loc(jInfoArray(iElemNumStartLoc))

c ======================================================================
c  Read material properties ONCE outside the loop
c ======================================================================
      E   = props(1)
      xnu = props(2)
      xl  = props(3)    ! length scale l0
      Gc  = props(4)

c  Lame constants
      mu     = E / (2.d0*(1.d0 + xnu))
      lambda = E*xnu / ((1.d0 + xnu)*(1.d0 - 2.d0*xnu))

c  Guard against incompressibility (xnu -> 0.5)
      if (abs(1.d0 - 2.d0*xnu) .lt. EPS_TOL) then
        lambda = 1.d20
      endif


c ======================================================================
c  Main loop over material points
c ======================================================================
      do i = 1, nblock
          
         dt  = dtArray(1)
         rho = density(i)
		 dtloc = dt

c --- 1. Read phase field phi from global array(pf_shared) ---
              intnum = jElemNum(i) 
              g = 1.d0
			  phi=0
              call vgetpartinfo(intnum, 1, cpname, locnum, jrcd)
              phi =    phi_store(locnum)
			  g      = (1.d0-phi)**2.d0 + 1.d-6
		
             

c --- 3. Read history variable for Irreversibility ---
          H_old = stateOld(i,7)

c ======================================================================
c  4. Update total strains  (accumulate from SDVs)
c     NOTE: strainInc from Abaqus uses engineering shear for nshr comp.
c ======================================================================
         eps11     = stateOld(i,1) + strainInc(i,1)
         eps22     = stateOld(i,2) + strainInc(i,2)
         eps33     = stateOld(i,3) + strainInc(i,3)  ! out-of-plane
         eps12_eng = stateOld(i,4) + strainInc(i,4)  ! tensor shear

c        Store updated strains in SDVs 1-4
         stateNew(i,1) = eps11
         stateNew(i,2) = eps22
         stateNew(i,3) = eps33
         stateNew(i,4) = eps12_eng

c        Tensor shear (VUMAT by defult give tensor shear) 
         eps12_t =  eps12_eng

c ======================================================================
c  5. Spectral decomposition of strain tensor (plane strain, 2D in-plane)
c     In-plane eigenvalues from 2x2 symmetric matrix:
c       [eps11  eps12_t]
c       [eps12_t  eps22 ]
c     Out-of-plane e3 = eps33 (plane strain => eps33=0 kinematically,
c     but kept general here for plane stress or user override)
c ======================================================================
         half_sum  = 0.5d0*(eps11 + eps22)
         half_diff = 0.5d0*(eps11 - eps22)
         rad       = sqrt(max(0.d0, half_diff**2 + eps12_t**2))

         e1 = half_sum + rad    ! larger in-plane principal strain
         e2 = half_sum - rad    ! smaller in-plane principal strain
         e3 = eps33             ! out-of-plane

c        Principal angle (angle of e1 axis w.r.t. global x)
c        Guard: if strain state is isotropic, angle is undefined -> 0
         if (abs(half_diff) .lt. EPS_TOL .and.
     1       abs(eps12_t)   .lt. EPS_TOL) then
            theta_p = 0.d0
         else
            theta_p = 0.5d0 * atan2(eps12_t, half_diff)
         endif

         cth = cos(theta_p)
         sth = sin(theta_p)

c ======================================================================
c  6. Positive/negative projection of principal strains
c     Miehe 2010: split eigenvalues into tensile (+) and compressive (-)
c ======================================================================
         e1p = max(e1, 0.d0)
         e2p = max(e2, 0.d0)
         e3p = max(e3, 0.d0)

         e1n = min(e1, 0.d0)
         e2n = min(e2, 0.d0)
         e3n = min(e3, 0.d0)

c        Trace decomposition — NOTE: clamp the SUM, not individual parts
c        <tr(eps)>+ = max(tr, 0),  <tr(eps)>- = min(tr, 0)
         tr_eps = e1 + e2 + e3
         tr_pos = max(tr_eps, 0.d0)
         tr_neg = min(tr_eps, 0.d0)

c ======================================================================
c  7. Principal stresses — Miehe spectral split
c
c   sigma_+ = lambda*<tr>+*I + 2*mu*eps_+   <- TENSILE  (degraded)
c   sigma_- = lambda*<tr>-*I + 2*mu*eps_-   <- COMPRESSIVE (NOT degraded)
c   sigma   = g * sigma_+ + sigma_-
c ======================================================================
         sig1p = lambda*tr_pos + 2.d0*mu*e1p
         sig2p = lambda*tr_pos + 2.d0*mu*e2p
         sig3p = lambda*tr_pos + 2.d0*mu*e3p

         sig1n = lambda*tr_neg + 2.d0*mu*e1n
         sig2n = lambda*tr_neg + 2.d0*mu*e2n
         sig3n = lambda*tr_neg + 2.d0*mu*e3n

c        Total principal stress
         sig1 = g*sig1p + sig1n
         sig2 = g*sig2p + sig2n
         sig3 = g*sig3p + sig3n

c ======================================================================
c  8. Rotate stress back to global frame
c     [s11 s12] = R * diag(sig1,sig2) * R^T
c     [s12 s22]
c     s11 = sig1*cth^2 + sig2*sth^2
c     s22 = sig1*sth^2 + sig2*cth^2
c     s12 = (sig1 - sig2)*cth*sth
c     s33 = sig3  (out-of-plane, no rotation)
c ======================================================================
         s11 = sig1*cth*cth + sig2*sth*sth
         s22 = sig1*sth*sth + sig2*cth*cth
         s12 = (sig1 - sig2)*cth*sth
         s33 = sig3

c        Store stresses (Abaqus order: 11,22,33,12 for plane strain)
         stressNew(i,1) = s11
         stressNew(i,2) = s22
         stressNew(i,3) = s33
         stressNew(i,4) = s12

c ======================================================================
c  9. Tensile strain energy density psi_plus (Miehe 2010, Eq. 18)
c
c   psi+ = (lambda/2)*<tr(eps)>+^2 + mu*(e1+^2 + e2+^2 + e3+^2)
c
c   IMPORTANT:
c   - Volumetric term uses tr_pos = <tr(eps)>+ = max(tr,0)  [clamp SUM]
c   - Deviatoric term uses e1p,e2p,e3p = individual clamps
c   These two operations are DIFFERENT and must NOT be confused.
c ======================================================================
         psi_plus = 0.5d0*lambda*tr_pos**2
     1            + mu*(e1p**2 + e2p**2 + e3p**2)

c ======================================================================
c  10. Update crack driving history  H = max(H_old, psi_plus)
c      Irreversibility: H never decreases
c ======================================================================
         H_new = max(H_old, min(psi_plus,1.0d10))
         stateNew(i,7) = H_new
         H_store(locnum) = H_new
		 
		 
		 
c ===================================================================
c =============== FRICTIONAL HEATING → (PART1) ================
c ===================================================================
            qdot_fric   = 0.d0
            deltaT_fric = 0.d0
            mu_fric = 0.5d0
      
        if (mu_fric .gt. 0.d0 .and. dtloc .gt. 0.d0
     &      .and. density(i) .gt. 0.d0) then

c ------- crack normal from ∇φ stored in Tgrad (phase-field) -------
             grad1 = 0.d0
             grad2 = 0.d0
             grad3 = 0.d0
             grad1 = Tgrad(locnum,1)
             grad2 = Tgrad(locnum,2)
             grad3 = Tgrad(locnum,3)
        

          gradnorm2 = grad1*grad1 + grad2*grad2 + grad3*grad3

          if (gradnorm2 .gt. 1.d-20) then
             gradnorm = sqrt(gradnorm2)
             n1 = -grad1/gradnorm
             n2 = -grad2/gradnorm
             n3 = -grad3/gradnorm
             t1 = -n2
             t2 =  n1

c ======================================================================
c  VELOCITY GRADIENT via F_new * F_old^{-1}
c  L = Fdot * F_inv  where  Fdot ≈ (F_new - F_old) / dt
c
c  Abaqus defgrad layout for 2D plane strain (ndir=3, nshr=1):
c    col 1=F11, 2=F22, 3=F33, 4=F12, 5=F21
c  Only the in-plane 2x2 block is used for the slip velocity calculations
c ======================================================================


         if (dtloc .gt. 0.d0) then

c --- extract in-plane components only ---
         F_new(1,1) = defgradNew(i,1)
         F_new(2,2) = defgradNew(i,2)
         F_new(1,2) = defgradNew(i,4)
         F_new(2,1) = defgradNew(i,5)

c --- FORCE F33=1 for plane strain regardless of what Abaqus returns ---
c --- CPE4R+nu=0.48 gives corrupted F33; plane strain mandates F33=1 ---
         F_old(1,1) = defgradOld(i,1)
         F_old(2,2) = defgradOld(i,2)
         F_old(1,2) = defgradOld(i,4)
         F_old(2,1) = defgradOld(i,5)

c --- in-plane determinant (F33=1 does not affect 2x2 block) ---
         detF = F_old(1,1)*F_old(2,2) - F_old(1,2)*F_old(2,1)

c --- additional sanity check for locking-corrupted F ---
         if (detF .lt. 0.5d0 .or. detF .gt. 2.0d0) then
c --- fallback to strainInc if F is untrustworthy ---
            stateNew(i,34)=detF
            L11 = strainInc(i,1) / dtloc
            L22 = strainInc(i,2) / dtloc
            L12 = strainInc(i,4)/(dtloc)+relSpinInc(i)/dtloc
            L21 = strainInc(i,4)/(dtloc)-relSpinInc(i)/dtloc
         else
            Finv(1,1) =  F_old(2,2)/detF
            Finv(1,2) = -F_old(1,2)/detF
            Finv(2,1) = -F_old(2,1)/detF
            Finv(2,2) =  F_old(1,1)/detF

            Fdot(1,1) = (F_new(1,1)-F_old(1,1))/dtloc
            Fdot(1,2) = (F_new(1,2)-F_old(1,2))/dtloc
            Fdot(2,1) = (F_new(2,1)-F_old(2,1))/dtloc
            Fdot(2,2) = (F_new(2,2)-F_old(2,2))/dtloc
		
            stateNew(i, 8) = F_new(1,1)
            stateNew(i, 9) = F_new(2,2)
            stateNew(i,10) = defgradNew(i,3)   ! F33 (out-of-plane)
            stateNew(i,11) = F_new(1,2)
            stateNew(i,12) = F_new(2,1)

            L11 = Fdot(1,1)*Finv(1,1) + Fdot(1,2)*Finv(2,1)
            L12 = Fdot(1,1)*Finv(1,2) + Fdot(1,2)*Finv(2,2)
            L21 = Fdot(2,1)*Finv(1,1) + Fdot(2,2)*Finv(2,1)
            L22 = Fdot(2,1)*Finv(1,2) + Fdot(2,2)*Finv(2,2)
         endif
		 endif

c --------- normal velocity gradient and slip speed vs -------------
                v1n = L11*n1 + L12*n2
                v2n = L21*n1 + L22*n2
                vt  = t1*v1n + t2*v2n
                vs  = abs(vt)*xl
			

c --------- normal traction and compression-only part --------------
                tn =  stressNew(i,1)*n1*n1
     &            +  stressNew(i,2)*n2*n2
     &            +  stressNew(i,3)*n3*n3
     &            + 2.d0*stressNew(i,4)*n1*n2
                tn_comp = max(-tn,0.d0)

c --------- phase-field crack density γ(c,∇c) ----------------------
c   γ = 1/(2ℓ) φ^2 + ℓ/2 |∇φ|^2   (area→volume conversion)
                gamma_pf = 0.5d0*( (phi*phi)/xl + xl*gradnorm2 )

c --------- volumetric frictional power qdot_fric ------------------
                qdot_fric =  mu_fric * tn_comp * gamma_pf * vs
				
          endif
         endif
             stateNew(i, 32) = gamma_pf
			 stateNew(i, 31) = tn_comp
			 stateNew(i, 30) =  vs 
			         if (dtloc .gt. 0.d0) then
            enerInelasNew(i) = enerInelasOld(i)
     &                        + qdot_fric*dtloc/density(i)
         else
            enerInelasNew(i) = enerInelasOld(i)
         endif


c ======================================================================
c  12. Internal energy (stress power)
c      For 2D plane strain: contributions from 11, 22, 33, 12
c      no Factor of 2 on shear because strainInc is tensor shear
c      but stress*strain requires tensor strain -> divide by 2 already
c      handled since Abaqus stressNew*strainInc convention is consistent
c ======================================================================
         stressPower =
     1     0.5d0*( stressOld(i,1)+stressNew(i,1) )*strainInc(i,1)
     2   + 0.5d0*( stressOld(i,2)+stressNew(i,2) )*strainInc(i,2)
     3   + 0.5d0*( stressOld(i,3)+stressNew(i,3) )*strainInc(i,3)
     4   +       ( stressOld(i,4)+stressNew(i,4) )*strainInc(i,4)

         if (totalTime .le. 1.d-12) then
            enerInternNew(i) = enerInternOld(i)
         else
            enerInternNew(i) = enerInternOld(i) + stressPower/rho
         endif

c ======================================================================
c  13. Carry forward any unused SDVs (5,6 reserved)
c ======================================================================
         stateNew(i,5) = stateOld(i,5)
         stateNew(i,6) = stateOld(i,6)
		 

      enddo   ! end loop over material points

      return
      end
      
      
      
c ======================================================================
c ========== VUMATHT - PHASE FIELD FRACTURE (φ AS TEMPERATURE) =========
c ======================================================================


      subroutine vumatht(
     *   nblock, nElem, nIntPt, nLayer, nSectPt,
     *   ntgrad, nstatev, nfieldv, nprops,
     *   cmname, stepTime, totalTime, dt,
     *   coordMp, density, props,
     *   tempOld, fieldOld, stateOld, enerThermOld,
     *   tempNew, tempgradNew, fieldNew,
     *   stateNew, fluxNew, enerThermNew, dEnerThDTemp, condEff )
      use pf_shared
      include 'vaba_param.inc'

      integer nblock, ntgrad, nstatev, nfieldv, nprops
      integer km, j
      character*80 cmname

      dimension nElem(nblock)
      dimension coordMp(nblock,*), density(nblock)
      dimension props(nprops)
      dimension tempOld(nblock), tempNew(nblock)
      dimension tempgradNew(nblock,ntgrad)
      dimension fieldOld(nblock,nfieldv), fieldNew(nblock,nfieldv)
      dimension stateOld(nblock,nstatev), stateNew(nblock,nstatev)
      dimension enerThermOld(nblock), enerThermNew(nblock)
      dimension fluxNew(nblock,ntgrad)
      dimension dEnerThDTemp(nblock,2), condEff(nblock)

c --- phase-field "thermal" parameters
      real*8 mobility, Cp
      real*8 phi, phi_old, g, gprime, H, xl, Gc
      real*8 Rsrc, dphi, dU_cap, dU_gen, dEnerTh
      real*8 rho , heat
      real*8 x,y,cpp, MISES_LIMIT
      integer kGlobal
      integer kOld, kNew
      real*8 LdotN1, LdotN2, LdotN_mag
c ====== FRICTION ADDITIONS: local variables ======
      real*8 mu_fric, Cv_mech,MISES
      real*8 grad1,grad2,grad3,gradnorm2,gradnorm
      real*8 n1,n2,n3,t1,t2,vs,vt,v1n,v2n
      real*8 F_old(2,2),F_new(2,2),Fdot(2,2),Finv(2,2),L(2,2),detF
      real*8 tn,tn_comp,gamma_pf,qdot_fric,deltaT_fric
c --- Crack surface density (per unit volume) ---
      real*8 phi_sq_term, grad_sq, crack_surf_density
      integer base
	  real*8 p_full, sdev11, sdev22, sdev33, sdev12
      real*8 v_normal, v_normal_sq, gradv_sq,dtloc  ! projection math

c --- for VGETPARTINFO debug ---
      character*80 cpname
      integer locnum, jrcd, intnum, jtyp



c ---------- PHASE-FIELD PARAMETERS (FROM MATERIAL NAME OR PROPS) ----------
c mobility = D_φ, Cp = capacity per mass C_p

           mobility = (1.0d0)
           Cp       =   props(1)
           Gc       =   props(2)   
           xl       =  4e-6
      
c ---------- LOOP OVER MATERIAL POINTS ----------
          do km = 1, nblock
           if (cmname(1:4).eq.'MAT4' .or.cmname(1:4).eq.'MAT5'.or.cmname(1:4).eq.'MAT6') then
		        phi_old = tempOld(km)
                phi = tempNew(km)
                dphi = phi - phi_old

c ----- get history variables and mechanical SDVs from module -----
c --- use analytic g and g'(φ) consistent with mechanical law ---
c g(φ) = (1-φ)^2 + tiny  => g'(φ) = -2(1-φ)

             g        = (1.d0-phi)**2.d0 + 1.d-6
             gprime   = -2.d0*(1.d0-phi)
             rho = density(km)
         
c ---------- FLUX: f = -mobility ∇φ -----------(mobility * tempgradNew(km,j))

               if (phi .ge. 1.d0) then
                  phi = 1.0d0
               endif
			   
         if (ntgrad .ge. 1) then
            do j = 1, ntgrad
               fluxNew(km,j) = -(mobility *tempgradNew(km,j))
            end do
         endif
         condEff(km) = mobility
         
c          intnum = nElem(km)
c         call vgetpartinfo(intnum, 1, cpname, locnum, jrcd), this also can be used inplace of nElem(km)-1000000
           
           locnum = nElem(km)-1000000
           H= H_store(locnum) 

c ---------- ENERGY INCREMENT (PER MASS) ----------
c U̇ = Cp φ̇ + [ φ/ℓ^2 + g'(φ) Gc H / ℓ ] / ρ
c → ΔU ≈ CpΔφ + [ (φ/ℓ^2 + g'(φ)GcH/ℓ) / ρ ] Δt

			     gprime = -2.d0*(1.d0-phi)
                    
			    if (phi .ge. 1.d0) then
                     Rsrc = phi/(xl*xl)
					 phi = 1.d0
                else
                     Rsrc = phi/(xl*xl) + gprime*H/(xl*Gc)
                endif
             
         dU_cap = Cp * dphi
         dU_gen = ((Rsrc) * dt)/(rho)
         dEnerTh =  dU_gen

         enerThermNew(km) = enerThermOld(km) + dEnerTh

c ---------- CAPACITY (PER MASS) ----------
         dEnerThDTemp(km,1) = Cp 
         dEnerThDTemp(km,2) = 0.d0
		 
c ---------- Update phase field gradient for VUMAT  ----------	
		 	  phi_store(locnum)=phi
	    	  Tgrad(locnum,1) = tempgradNew(km,1)
              Tgrad(locnum,2) = tempgradNew(km,2)
			  
      endif 
      
      end do

      return
      end
