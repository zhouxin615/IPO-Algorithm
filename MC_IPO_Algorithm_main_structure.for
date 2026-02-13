!======================================================================
! ModelParam - global model and algorithm parameters
!======================================================================
      module ModelParam
      implicit none

      !----------------------------------------------------------------
      ! Material (model) parameters (global)
      !----------------------------------------------------------------
      double precision EYoung, nu, phi_fric, pi, softenCoeff, Beta, cd
      double precision cohesion0, eps_smooth, beta_smooth, psi_dilat,
     &                 sqrt3,
     &                 sin_phi_fric, cos_phi_fric, tan_phi_fric,
     &                 tan_psi_dilat,
     &                 M_slope, K0_slope, gamw, gam, aw
      double precision dstrain_global(6,1)

      !----------------------------------------------------------------
      ! Algorithm parameters (global)
      !----------------------------------------------------------------
      double precision Ftol, hsize, rho, zeta
      integer nten, nvar, Newton, MaxLSM

      end module ModelParam


!======================================================================
! InitMaterial - read material parameters from PROPS and preprocess
!======================================================================
      subroutine InitMaterial(PROPS)
      use ModelParam
      implicit none

      double precision PROPS(*)

      ! 1. mathematical constants
      pi    = 3.1415926535897934d0
      sqrt3 = dsqrt(3.d0)

      ! 2. read raw material parameters from PROPS
      EYoung      = PROPS(1)
      nu          = PROPS(2)
      cohesion0   = PROPS(3)
      phi_fric    = PROPS(4)*pi/180.d0
      psi_dilat   = PROPS(5)*pi/180.d0

      ! numerical regularization parameter
      softenCoeff = 1.d-8

      ! 3. trig functions and derived constants
      sin_phi_fric = dsin(phi_fric)
      cos_phi_fric = dcos(phi_fric)
      tan_phi_fric = dtan(phi_fric)
      tan_psi_dilat= dtan(psi_dilat)

      M_slope  = 6.d0*sin_phi_fric/sqrt3/(3.d0 - sin_phi_fric)
      K0_slope = 6.d0*cos_phi_fric/sqrt3/(3.d0 - sin_phi_fric)

      gamw = 6.d0/pi*datan(sin_phi_fric/sqrt3)
      gam  = 1.d0 - gamw
      aw   = 1.d0/(dcos((gamw+1.d0)*pi/6.d0))

      ! 4. smooth yield surface parameters
      eps_smooth  = 0.1d0
      beta_smooth = 0.9999d0

      return
      end subroutine InitMaterial


!======================================================================
! UMAT - Abaqus user material subroutine
!======================================================================
      SUBROUTINE UMAT (STRESS,statev,DDSDDE,SSE,SPD,SCD,RPL,
     1 DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,
     2 PREDEF,DPRED,CMNAME,NDI,NSHR,ntens,NSTATV,PROPS,NPROPS,
     3 COORDS,DROT,PNEWDT,CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,
     4 KSPT,KSTEP,KINC)

      use ModelParam
      implicit none

      CHARACTER*80 CMNAME

      !----------------------------------------------------------------
      ! Abaqus interface variables (MUST keep the declaration pattern)
      !----------------------------------------------------------------
      integer ntens, NDI, NSHR, NSTATV, NPROPS, NOEL, NPT,
     1        LAYER, KSPT, KSTEP, KINC

      double precision STRESS(ntens), statev(NSTATV),
     1    DDSDDE(ntens,ntens), DDSDDT(ntens), DRPLDE(ntens),
     2    STRAN(ntens), DSTRAN(ntens), TIME(2), PREDEF(1), DPRED(1),
     3    PROPS(NPROPS), COORDS(3), DROT(3,3),
     4    DFGRD0(3,3), DFGRD1(3,3),
     5    SSE, SPD, SCD, RPL, DRPLDT, TEMP, DTEMP,
     6    PNEWDT, CELENT, DTIME

      !----------------------------------------------------------------
      ! Internal variables for constitutive integration
      !----------------------------------------------------------------
      double precision epn, normg, p, q, theta, J2, J3
      double precision dt_sub, t_subsum, dtmin
      double precision CTO(ntens,ntens), g(ntens+2,1),
     &                 De(ntens,ntens)
      double precision r(ntens,1), sigma(ntens,1),
     &                 sigma_trial(ntens,1),
     &                 sigma_dev(ntens,1), dsigma(ntens,1)
      double precision xold(ntens+2,1), xnew(ntens+2,1),
     &                 xnewold(ntens+2,1)
      double precision Jacobian(ntens+2,ntens+2),
     &                 invJacobian(ntens+2,ntens+2)
      double precision dep(ntens,1), dstrain_sum(6,1), Ape_index

      integer Iteration, Isucc, Isubstep, MaxSub, Iter_sum

      !----------------------------------------------------------------
      ! 0. Basic dimensions
      !----------------------------------------------------------------
      nten = ntens
      nvar = ntens + 2

      !----------------------------------------------------------------
      ! 1. Algorithmic parameters
      !----------------------------------------------------------------
      Newton    = 50
      MaxLSM    = 30
      Iteration = 0
      Iter_sum = 0
      
      hsize = 1.d-6
      rho   = 1.d-4
      zeta  = 0.1d0

      Beta   = 1.d-20
      cd     = 1.d0
      Ftol   = 1.d-6
      dtmin  = 1.d-8
      MaxSub = 200

      !----------------------------------------------------------------
      ! 2. Read material parameters and precompute trigonometric terms
      !----------------------------------------------------------------
      call InitMaterial(PROPS)

      !----------------------------------------------------------------
      ! 3. Read stress from previous step and total strain increment
      !----------------------------------------------------------------
      sigma          = 0.d0
      dstrain_global = 0.d0
      call DimIn(STRESS, DSTRAN, sigma)

      ! Initial step: apply a very small compressive stress
      if (TIME(2) .LT. 1.d-10) then
        statev      = 0.d0
        statev(1:3) = 1.d-8/3.d0
        statev(7)   = dsqrt((1.d-8/3.d0)**2.d0 * 3.d0)
        statev(13)  = 0.d0
      end if

      ! Equivalent plastic strain at previous step
      epn = statev(7)

      ! Unknown vector at previous step:
      ! xold = [sigma ; epn ; dphi] with dphi = 0 at previous step
      xold(1:ntens,1) = sigma(1:ntens,1)
      xold(ntens+1,1) = epn
      xold(ntens+2,1) = 0.d0

      ! Elastic stiffness matrix De at current step
      call Defun(xold, De)

      !----------------------------------------------------------------
      ! 4. Elastic trial stress increment
      !----------------------------------------------------------------
      dsigma      = matmul(De, dstrain_global)
      sigma_trial = sigma + dsigma

      xnew             = xold
      xnew(1:ntens,1)  = sigma_trial(1:ntens,1)

      !----------------------------------------------------------------
      ! 5. Initialization for substepping control
      !----------------------------------------------------------------
      Isucc = 0
      Isubstep    = 0
      dt_sub      = 1.d0
      t_subsum    = 0.d0

      dstrain_sum    = dstrain_global
      dstrain_global = 0.d0

      !----------------------------------------------------------------
      ! 6. Substep loop: nonlinear equation solving
      !----------------------------------------------------------------
      do while ((t_subsum .LT. 1.d0) .and. (Isubstep .LE. MaxSub))

        Isubstep = Isubstep + 1
        t_subsum = t_subsum + dt_sub

        ! Current substep strain increment (proportional split)
        dstrain_global = t_subsum * dstrain_sum

        ! Store a "good" initial guess for this substep
        xnewold = xnew
        if (Isubstep .EQ. 1) xnewold = xold

        ! Line-search Newton iteration
        call LSMfun(xnew, xold, Iteration, De, g, normg)

        ! Convergence check and adaptive substepping
        if ((Iteration .GE. Newton) .or. isnan(normg)) then
          ! Substep fails: roll back and reduce dt
          t_subsum = t_subsum - dt_sub
          dt_sub   = dt_sub * 0.25d0
          xnew     = xnewold
          !Isucc = 0
        else
          ! Substep succeeds: gradually increase dt
          Isucc = Isucc + 1
          Iter_sum = Iter_sum + Iteration
          !if (Isucc .GE. 2) dt_sub = dt_sub * 1.5d0
          dt_sub = dt_sub * 1.5d0
          
        end if

        ! Apply lower bound and remaining fraction limits
        dt_sub = max(dt_sub, dtmin)
        dt_sub = min(dt_sub, (1.d0 - t_subsum))

      end do

      !----------------------------------------------------------------
      ! 7. Compute consistent tangent operator, 只有在计算切线刚度矩阵时，使用中心差分
      !----------------------------------------------------------------
      call Jacobfun(xnew, xold, De, Jacobian)
      call CTO_LU(Jacobian, De, CTO, ntens, nvar)

      !----------------------------------------------------------------
      ! 8. Update stress and internal variables
      !----------------------------------------------------------------
      call rfun(xnew, r)

      STRESS(1:ntens) = xnew(1:ntens,1)

      dep(1:ntens,1)  = r(1:ntens,1) * xnew(ntens+2,1)
      statev(1:ntens) = statev(1:ntens) + dep(1:ntens,1)
      statev(15) = xnew(ntens+1,1) - statev(7) ! 等效塑性剪应变增量
      statev(7)  = xnew(ntens+1,1)
      statev(8)  = Iteration
      statev(9)  = normg
      statev(10) = statev(10) + Isubstep
      statev(11) = Isubstep
      statev(12) = Isucc
      statev(13) = Iter_sum
      statev(14) = Ape_index(xnew)
      
      
      !----------------------------------------------------------------
      ! 9. If Newton iteration diverges, suggest Abaqus to reduce step
      !----------------------------------------------------------------
      if ((Iteration .GE. Newton) .or. isnan(normg) .or.
     2    isnan(norm2(CTO))     .or. isnan(norm2(xnew))) then

        PNEWDT = 0.25d0

        !write(*,*) 'UMAT warning: step difficulty detected.'
        !write(*,*) 'time(2) = ', TIME(2)
        !write(*,*) 'tsum    = ', t_subsum
        !write(*,*) 'Iteration, normg = ', Iteration, normg
        !write(*,*) 'dstrain_global = ', dstrain_global
        !write(*,*) 'xnew    = ', xnew
        !write(*,*) 'g       = ', g
        !
        !! Print p-q-theta information for current and previous states
        !sigma(1:ntens,1) = xnew(1:ntens,1)
        !call pqJ2J3(sigma, p, q, sigma_dev, theta, J2, J3)
        !write(*,*) 'new  p, q, theta(deg) = ', p, q, theta*180.d0/pi
        !
        !write(*,*) 'xold   = ', xold
        !sigma(1:ntens,1) = xold(1:ntens,1)
        !call pqJ2J3(sigma, p, q, sigma_dev, theta, J2, J3)
        !write(*,*) 'old  p, q, theta(deg) = ', p, q, theta*180.d0/pi

        CTO  = De
        xnew = xold
        !CTO  = 1.d100
c       call XIT   ! Uncomment to force Abaqus to stop immediately
      end if

      !----------------------------------------------------------------
      ! 10. Dimension conversion and return to Abaqus
      !----------------------------------------------------------------
      call DimOut(DDSDDE, STRESS, CTO, xnew)

      return
      end SUBROUTINE UMAT


!======================================================================
! DimIn - convert stress and strain increment to internal format
!======================================================================
      subroutine DimIn(stress, dstran, sigma)
      use ModelParam
      implicit none

      double precision stress(nten), dstran(nten), sigma(nten,1)

      if ((nten .eq. 6) .or. (nten .eq. 4)) then
        sigma(1:nten,1)          = stress(1:nten)
        dstrain_global(1:nten,1) = dstran(1:nten)
      elseif (nten .eq. 3) then
        write(*,*) 'UMAT error: plane stress (ntens=3) not supported.'
        !call XIT
      end if

      return
      end subroutine DimIn


!======================================================================
! DimOut - write stress and stiffness back to Abaqus in proper size
!======================================================================
      subroutine DimOut(DDSDDE, stress, CTO, xnew)
      use ModelParam
      implicit none

      double precision DDSDDE(nten,nten), stress(nten),
     &                 CTO(nten,nten), xnew(nvar,1)

      if ((nten .eq. 6) .or. (nten .eq. 4)) then
        stress(1:nten) = xnew(1:nten,1)
        DDSDDE         = CTO
      elseif (nten .eq. 3) then
        write(*,*) 'UMAT error: plane stress (ntens=3) not supported.'
        !call XIT
      end if

      return
      end subroutine DimOut
      
!======================================================================
! The rest of the code will be uploaded after the paper is published
!======================================================================
