c ======== SHARED MODULE (put in a separate file if required: pf_shared.f) ========
             module pf_shared
      implicit none

      integer, parameter :: MAX_ELEM = 5000000

      real*8 H_store(MAX_ELEM)
      real*8 g_store(MAX_ELEM)
      logical initialized(MAX_ELEM)
	     real*8 phi_store(MAX_ELEM)
      real*8 Tgrad(MAX_ELEM,3)
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
      integer jInfoArray(nblock)

      character*80 cmname
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),dtArray(*),
     4  stretchOld(*), defgradOld(nblock,ndir+nshr+1),
     5  fieldOld(nblock,nfieldv),
     6  stressOld(nblock,ndir+nshr),
     7  stateOld(nblock,nstatev), enerInternOld(nblock),
     8  enerInelasOld(nblock), tempNew(*),
     9  stretchNew(*), defgradNew(nblock,ndir+nshr+1),
     &  fieldNew(nblock,nfieldv),
     &  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     &  enerInternNew(nblock), enerInelasNew(nblock)

      dimension jels(nblock)

c --- integer locals ---
      integer i,j,constitutive_type,kGlobal
      integer elemStart

c --- material constants ---
      real*8 E,xnu,xl,Gc
      real*8 twomu,eg,sixmu,alamda,lambda_lame,mu_lame
      real*8 xK

c --- strain array ---
      real*8 stran(nblock,ndir+nshr)

c --- eigenvalue / spectral work ---
      real*8 eps11,eps22,eps33,eps12
      real*8 trace2,diff,rad,eig1,eig2,eig3
      real*8 pos1,pos2,pos3,neg1,neg2,neg3,tr_pos,tr_neg
      real*8 sig1_tens,sig2_tens,sig3_tens
      real*8 sig1_comp,sig2_comp,sig3_comp
      real*8 sig1_tot,sig2_tot,sig3_tot
      real*8 theta_p,cth,sth
      real*8 s11_tot,s22_tot,s12_tot,s33_tot
      real*8 eigVal(3),trp1,trp2,psip,psit,H

c --- energy / time ---
      real*8 stressPower,dtloc,dt,rho

c --- phase field ---
      real*8 phi,g,p
      real*8 x,y
      real*8 s1,s2,s3,s4
      real*8 p_tens,p_comp
      real*8 pi,base
      real*8 del_status,J_est,eps1_lim,J_lim,gamma_lim
      character*80 cpname

c --- friction ---
      real*8 mu_fric,Cv_mech
      real*8 grad1,grad2,grad3,gradnorm2,gradnorm
      real*8 n1,n2,n3,t1,t2,vs,vt,v1n,v2n
      real*8 F_old(2,2),F_new(2,2),Fdot(2,2),Finv(2,2)
      real*8 L(2,2),detF
      real*8 tn,tn_comp,gamma_pf,qdot_fric,deltaT_fric

c --- stress/strain work (undamaged, deviatoric) ---
      real*8 sig0_11,sig0_22,sig0_33,sig0_12
      real*8 p_full
      real*8 sdev11,sdev22,sdev33,sdev12
      real*8 eps12t
      real*8 eDevSq
      real*8 trStrain,trPos,trNeg
      real*8 epsVol,edev11,edev22,edev33
      real*8 psiPlus,psiMinus

c --- deviatoric tensile/compressive split in global frame ---
      real*8 s1_tens,s2_tens,s4_tens
      real*8 s1_comp,s2_comp,s4_comp
c --- new local variable needed ---
      real*8 xKvol
    

c --- numerical parameters ---
      real*8 L11,L12,L21,L22
      real*8 L_CAP
      parameter (L_CAP    = 1.0d8)
      parameter (eps1_lim = 0.35d0)
      parameter (J_lim    = 0.05d0)
      parameter (gamma_lim= 0.50d0)
      parameter (pi       = 3.141592653589d0)
      parameter (oneThird = 1.d0/3.d0)

c ======================================================================
c  Pointer setup and material constants
c ======================================================================
      iElemNumStartLoc = jInfoArray(i_info_ElemNumStartLoc)
      ptrjElemNum      = loc(jInfoArray(iElemNumStartLoc))

      E    = props(1)
      xnu  = props(2)
      xl   = props(3)
      Gc   = props(4)

      twomu  = E / (1.d0 + xnu)
      eg     = twomu / 2.d0
      alamda = twomu*(E - twomu) / (3.d0*twomu - 2.d0*E)
      xK     = alamda + twomu*oneThird
	  
c ======================================================================
c  Select volumetric coefficient by material:
c   MAT1        -> normal/compressible material -> use true bulk modulus K
c   MAT2, MAT3  -> near-incompressible material  -> lambda ~ K, use lambda
c ======================================================================
      if (cmname(1:4).eq.'MAT2' .or. cmname(1:4).eq.'MAT3') then
         xKvol = alamda
      else
         xKvol = xK
      endif

c ======================================================================
c  Main loop over material points
c ======================================================================
      do i = 1, nblock

         psit  = stateOld(i,7)
         dt    = dtArray(1)
         rho   = density(i)
         dtloc = dt
         intnum = jElemNum(i)

         call vgetpartinfo(intnum, 1, cpname, locnum, jrcd)
c --- g and phi: single guarded block, no redundant reads ---
         if (dtArray(1) .le. 1.d-14) then
            g   = 1.d0
            phi = 0.d0
         else
            phi = phi_store(locnum)
            g   = (1.d0 - phi)**2 + 1.d-6
         endif

C ======================================================================
C  1. Accumulate total strains
C ======================================================================
         stran(i,1) = stateOld(i,8)  + strainInc(i,1)
         stran(i,2) = stateOld(i,9)  + strainInc(i,2)
         stran(i,3) = stateOld(i,10) + strainInc(i,3)
         stran(i,4) = stateOld(i,11) + strainInc(i,4)

         stateNew(i,8)  = stran(i,1)
         stateNew(i,9)  = stran(i,2)
         stateNew(i,10) = stran(i,3)
         stateNew(i,11) = stran(i,4)
         stateNew(i,12) = stateOld(i,12)
         stateNew(i,13) = stateOld(i,13)

C ======================================================================
C  2. Undamaged stress from TOTAL strain — lambda,mu only
C     sigma0 = lambda*tr(eps)*I + 2*mu*eps
C     NO incremental update, NO K used anywhere
C ======================================================================
         trStrain = stran(i,1) + stran(i,2) + stran(i,3)

         sig0_11 = alamda*trStrain + twomu*stran(i,1)
         sig0_22 = alamda*trStrain + twomu*stran(i,2)
         sig0_33 = alamda*trStrain + twomu*stran(i,3)
         sig0_12 = twomu * stran(i,4)

C ======================================================================
C  3. Strain energy split — consistent with deviatoric tensile/comp split
C
C     Standard Amor:
C       psi+ = 0.5*lambda*<tr>+^2  +  mu*(e':e')   <- ALL dev in psi+
C
C     Modified (consistent with stress degradation):
C       psi+ = 0.5*lambda*<tr>+^2  +  mu*(e'_pos:e'_pos)
C       psi- = 0.5*lambda*<tr>-^2  +  mu*(e'_neg:e'_neg)
C
C     where e'_pos and e'_neg are the tensile and compressive
C     principal deviatoric strains rotated back — exactly matching
C     what you preserve/degrade in the stress update.
C
C     This makes H physically consistent:
C       H only grows when tensile volumetric OR tensile deviatoric
C       strain energy increases — compressive shear does NOT drive H.
C ======================================================================

C --- Macaulay on volumetric strain
         trPos = max(trStrain, 0.d0)
         trNeg = min(trStrain, 0.d0)

C --- Deviatoric strain tensor components (tensor shear)
         edev11 = stran(i,1) - oneThird*trStrain
         edev22 = stran(i,2) - oneThird*trStrain
         edev33 = stran(i,3) - oneThird*trStrain
         eps12t = stran(i,4)

C --- Eigendecompose deviatoric STRAIN tensor [edev11 eps12t; eps12t edev22]
C     Same procedure as stress split — applied to strain here
         trace2 = 0.5d0*(edev11 + edev22)
         diff   = sqrt( max(0.d0,
     1            (0.5d0*(edev11-edev22))**2 + eps12t**2) )

         eig1 = trace2 + diff     ! larger  in-plane deviatoric strain
         eig2 = trace2 - diff     ! smaller in-plane deviatoric strain
         eig3 = edev33            ! out-of-plane (already principal)

C --- Principal angle of deviatoric strain tensor
         if (diff .gt. 1.d-20) then
            theta_p = 0.5d0 * atan2(eps12t, 0.5d0*(edev11-edev22))
         else
            theta_p = 0.d0
         endif
         cth = cos(theta_p)
         sth = sin(theta_p)

C --- Split deviatoric strain eigenvalues
         pos1 = max(eig1, 0.d0)
         pos2 = max(eig2, 0.d0)
         pos3 = max(eig3, 0.d0)
         neg1 = min(eig1, 0.d0)
         neg2 = min(eig2, 0.d0)
         neg3 = min(eig3, 0.d0)

C --- Tensile deviatoric strain double contraction e'_pos:e'_pos
C     = pos1^2 + pos2^2 + pos3^2
C     (eigenvalues are orthonormal so no cross terms)

         eDevSq = pos1**2 + pos2**2 + pos3**2

C --- Compressive deviatoric strain double contraction e'_neg:e'_neg
C     (stored separately — used only for psiMinus)

         epsVol = neg1**2 + neg2**2 + neg3**2

C --- Consistent strain energy split
C     psi+ : tensile vol + tensile dev only  -> drives H
C     psi- : compressive vol + compressive dev -> preserved
    
		 psiPlus  = 0.5d0*xKvol*trPos**2 + eg*eDevSq
         psiMinus = 0.5d0*xKvol*trNeg**2 + eg*epsVol

         psip = psiPlus
         H = max(psit, min(psip, 1.0d10))
         H_store(locnum) = H

C ======================================================================
C  4. History field H
C     Now H only grows from tensile volumetric + tensile deviatoric
C     Compressive shear does NOT drive crack growth
C ======================================================================
         H = max(psit, min(psip, 1.0d10))
         H_store(locnum) = H

C ======================================================================
C  5. Stress degradation
C     Reuse pos1,pos2,pos3,neg1,neg2,neg3,cth,sth from strain split
C     above — same principal frame for stress since linear elasticity
C     means stress and strain share eigenvectors
C ======================================================================

         stateNew(i,1) = sig0_11
         stateNew(i,2) = sig0_22
         stateNew(i,3) = sig0_33
         stateNew(i,4) = sig0_12
         stateNew(i,5) = 0.d0
         stateNew(i,6) = 0.d0

         p_full = ( sig0_11 + sig0_22 + sig0_33 ) * oneThird
         p_tens = max(p_full, 0.d0)
         p_comp = min(p_full, 0.d0)

         s1 = sig0_11 - p_full
         s2 = sig0_22 - p_full
         s3 = sig0_33 - p_full
         s4 = sig0_12

C --- Eigendecompose deviatoric STRESS [s1 s4; s4 s2]
C     In linear elasticity: stress = 2*mu*strain (deviatoric part)
C     so eigenvectors are identical — only eigenvalues differ by 2*mu
C     We still decompose stress directly for exactness
         trace2 = 0.5d0*(s1 + s2)
         diff   = sqrt( max(0.d0, (0.5d0*(s1-s2))**2 + s4**2) )
         eig1   = trace2 + diff
         eig2   = trace2 - diff
         eig3   = s3

C --- Principal angle for stress deviatoric
         if (diff .gt. 1.d-20) then
            theta_p = 0.5d0 * atan2(s4, 0.5d0*(s1-s2))
         else
            theta_p = 0.d0
         endif
         cth = cos(theta_p)
         sth = sin(theta_p)

         pos1 = max(eig1, 0.d0)
         pos2 = max(eig2, 0.d0)
         pos3 = max(eig3, 0.d0)
         neg1 = min(eig1, 0.d0)
         neg2 = min(eig2, 0.d0)
         neg3 = min(eig3, 0.d0)

         s1_tens = pos1*cth*cth + pos2*sth*sth
         s2_tens = pos1*sth*sth + pos2*cth*cth
         s4_tens = (pos1 - pos2)*cth*sth

         s1_comp = neg1*cth*cth + neg2*sth*sth
         s2_comp = neg1*sth*sth + neg2*cth*cth
         s4_comp = (neg1 - neg2)*cth*sth

         stressNew(i,1) = p_comp + s1_comp + g*(p_tens + s1_tens)
         stressNew(i,2) = p_comp + s2_comp + g*(p_tens + s2_tens)
         stressNew(i,3) = p_comp + neg3    + g*(p_tens + pos3   )
         stressNew(i,4) =          s4_comp + g* s4_tens

         stateNew(i,7) = H

C ======================================================================
C  6. Internal energy — stress power per unit mass
C     Trapezoidal: 0.5*(sig_old+sig_new)*deps
C     s33 included — nonzero in plane strain
C     First increment: stressOld=0 -> stressPower=0.5*sig_new*deps
C ======================================================================
         stressPower =
     1     0.5d0*(stressOld(i,1)+stressNew(i,1))*strainInc(i,1)
     2   + 0.5d0*(stressOld(i,2)+stressNew(i,2))*strainInc(i,2)
     3   + 0.5d0*(stressOld(i,3)+stressNew(i,3))*strainInc(i,3)
     4   + (stressOld(i,4)+stressNew(i,4))*strainInc(i,4)

          if (dt .le. 1.d-13) then
            enerInternNew(i) = enerInternOld(i)
         else
            enerInternNew(i) = enerInternOld(i) + stressPower/rho
         endif

c ======================================================================
c  FRICTIONAL HEATING — only for Part-1 
c ======================================================================
               qdot_fric   = 0.d0
               deltaT_fric = 0.d0
			   mu_fric = 0.5 
         if (mu_fric .gt. 0.d0 .and. dtloc .gt. 0.d0
     &       .and. density(i) .gt. 0.d0) then

            grad1 = Tgrad(locnum,1)
            grad2 = Tgrad(locnum,2)
            grad3 = Tgrad(locnum,3)

            gradnorm2 = grad1*grad1 + grad2*grad2 + grad3*grad3

            if (gradnorm2 .gt. 1.d-20) then
               gradnorm = sqrt(gradnorm2)
               n1 =  -grad1/gradnorm
               n2 =  -grad2/gradnorm
               n3 =  -grad3/gradnorm
               t1 =  -n2
               t2 =   n1

               if (dtloc .gt. 0.d0) then

                  F_new(1,1) = defgradNew(i,1)
                  F_new(2,2) = defgradNew(i,2)
                  F_new(1,2) = defgradNew(i,4)
                  F_new(2,1) = defgradNew(i,5)

                  F_old(1,1) = defgradOld(i,1)
                  F_old(2,2) = defgradOld(i,2)
                  F_old(1,2) = defgradOld(i,4)
                  F_old(2,1) = defgradOld(i,5)

                  detF = F_old(1,1)*F_old(2,2)
     &                 - F_old(1,2)*F_old(2,1)

                  if (detF .lt. 0.5d0 .or. detF .gt. 2.0d0) then
                     L11 = strainInc(i,1)/dtloc
                     L22 = strainInc(i,2)/dtloc
                     L12 = strainInc(i,4)/(dtloc)
     &                   + relSpinInc(i)/dtloc
                     L21 = strainInc(i,4)/(dtloc)
     &                   - relSpinInc(i)/dtloc
                  else
                     Finv(1,1) =  F_old(2,2)/detF
                     Finv(1,2) = -F_old(1,2)/detF
                     Finv(2,1) = -F_old(2,1)/detF
                     Finv(2,2) =  F_old(1,1)/detF

                     Fdot(1,1) = (F_new(1,1)-F_old(1,1))/dtloc
                     Fdot(1,2) = (F_new(1,2)-F_old(1,2))/dtloc
                     Fdot(2,1) = (F_new(2,1)-F_old(2,1))/dtloc
                     Fdot(2,2) = (F_new(2,2)-F_old(2,2))/dtloc

                     L11 = Fdot(1,1)*Finv(1,1)+Fdot(1,2)*Finv(2,1)
                     L12 = Fdot(1,1)*Finv(1,2)+Fdot(1,2)*Finv(2,2)
                     L21 = Fdot(2,1)*Finv(1,1)+Fdot(2,2)*Finv(2,1)
                     L22 = Fdot(2,1)*Finv(1,2)+Fdot(2,2)*Finv(2,2)
                  endif
               endif

               v1n = L11*n1 + L12*n2
               v2n = L21*n1 + L22*n2
               vt  = t1*v1n + t2*v2n
               vs  = abs(vt)*xl

               tn = stressNew(i,1)*n1*n1
     &            + stressNew(i,2)*n2*n2
     &            + stressNew(i,3)*n3*n3
     &            + 2.d0*stressNew(i,4)*n1*n2
               tn_comp = max(-tn, 0.d0)

               gamma_pf = 0.5d0*((phi*phi)/xl + xl*gradnorm2)

               qdot_fric = mu_fric * tn_comp * gamma_pf * vs

            endif
         endif
         stateNew(i,32)  = gamma_pf
         stateNew(i,31)  = tn_comp
         stateNew(i,30)  = vs
      
		          if (dtloc .gt. 0.d0) then
            enerInelasNew(i) = enerInelasOld(i)
     &                        + qdot_fric*dtloc/density(i)
         else
            enerInelasNew(i) = enerInelasOld(i)
         endif

      enddo

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
            xl =4e-6
		 
		 
c ---------- LOOP OVER MATERIAL POINTS ----------
          do km = 1, nblock
           if (cmname(1:4).eq.'MAT4' .or.cmname(1:4).eq.'MAT5'.or.cmname(1:4).eq.'MAT6') then
             phi_old = tempOld(km)
             phi = tempNew(km)
             dphi = phi - phi_old

c ----- get history variables and mechanical SDVs from module -----
c --- use analytic g and g'(φ) consistent with mechanical law ---
c g(φ) = (1-φ)^2 + tiny  => g'(φ) = -2(1-φ)
          
		  g      = (1.d0-phi)**2.d0 + 1.d-6
          gprime = -2.d0*(1.d0-phi)
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
c         call vgetpartinfo(intnum, 1, cpname, locnum, jrcd)
           
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

c ---------- Upadte these varibale for frictional heat in VUMAT ----------
		 	   phi_store(locnum)=phi
	    	   Tgrad(locnum,1) = tempgradNew(km,1)
               Tgrad(locnum,2) = tempgradNew(km,2)
		    
      endif 
      
      end do

      return
      end
