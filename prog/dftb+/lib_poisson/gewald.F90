!**************************************************************************
!  Copyright (c) 2004 by Univ. Rome 'Tor Vergata'. All rights reserved.   *  
!  Authors: A. Pecchia, L. Latessa, A. Di Carlo                           *
!                                                                         *
!  Permission is hereby granted to use, copy or redistribute this program * 
!  under the LGPL licence.                                                *
!**************************************************************************

#:include 'error.fypp'

module gewald

  use dftbp_accuracy, only : dp
  use dftbp_constants, only : pi
  implicit none

  private
  
  public :: short_pot
  public :: long_pot
  public :: getalpha
  public :: cross,rezvol,phi,phi1


contains

  !======================================================================
  subroutine short_pot(distR,basis,uhatm,deltaQ,tol,sh_pot, iError)

    real(kind=dp) :: distR(3), uhatm, deltaQ, basis(3,3), tol, sh_pot
    integer, intent(out), optional :: iError

    integer i,j,k,nreal,nmax,nmin 
    real(kind=dp) :: rvec(3),R(3),lastshell,tmp,norm

    if (present(iError)) then
      iError = 0
    end if

    sh_pot = 0.0_dp
    nmax = 50
    nmin = 3
    nreal = 0
    lastshell = tol+1e-8
    ! /* sum over R until tolerance is reached */
    DO WHILE ((nreal .le. nmax) .and. ((abs(lastshell) .gt. tol).or. & 
        (nreal .le. nmin)  ))
      lastshell = 0.0
      DO i = -nreal,nreal
        DO j = -nreal,nreal
          DO k = -nreal,nreal
            ! /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j))  .or. &
              &         (nreal.eq. abs(k)) ) THEN

              R(1)=i*basis(1,1)+j*basis(2,1)+k*basis(3,1) 
              R(2)=i*basis(1,2)+j*basis(2,2)+k*basis(3,2) 
              R(3)=i*basis(1,3)+j*basis(2,3)+k*basis(3,3) 

              rvec(1) = distR(1) - R(1)  
              rvec(2) = distR(2) - R(2)   
              rvec(3) = distR(3) - R(3)  

              norm = sqrt(rvec(1)**2+rvec(2)**2+rvec(3)**2)

              IF (norm.gt.1.0e-20) THEN

                tmp = deltaQ*(exp(-3.2_dp*uhatm*norm))*(1.6_dp*uhatm+1/norm)
                sh_pot = sh_pot + tmp
                lastshell = lastshell + tmp

              ELSE 

                lastshell = tol+1e-8

              END IF

            END IF
          END DO
        END DO
      END DO

      nreal = nreal + 1

    END DO

    IF(abs(lastshell) .gt. tol) THEN
      @:ERROR_HANDLING(iError, -1, 'tolerance in subroutine short_pot not reached')
    END IF

  END subroutine short_pot
  !----------------------------------------------------------------     


  subroutine long_pot(r,basis,recbasis,alpha,vol,tol,potential, iError)

    real(kind=dp) ::  r(3), basis(3,3), recbasis(3,3), alpha, vol, tol
    real(kind=dp) ::  potential
    integer, intent(out), optional :: iError

    real(kind=dp) ::  reciprocal,rcspace,cterm
    real(kind=dp) ::  G(3),rh(3),help,norm,lastshell
    integer nrezi, nreal, nmax, nmin
    integer i,j,k
    nmax = 20
    nmin = 2

    if (present(iError)) then
      iError = 0
    end if

    !evaluate reciprocal space term ( sum over G <> 0) ...  
    !/* sum over G until tolerance is reached */
    nrezi = 1
    lastshell = tol+1.d-8  
    reciprocal = 0.0_dp
    DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or. &
      &   (abs(lastshell) .gt.  tol)))
      lastshell = 0.0_dp
      DO i=-nrezi,nrezi
        DO j=-nrezi,nrezi
          DO k=-nrezi,nrezi
            !/*only G belonging to outer shells are new ones */
            IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or. &
              &        (nrezi.eq. abs(k)) ) THEN
              G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
              G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
              G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)

              help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
              help = exp(-help/(4._dp*alpha*alpha))/help
              help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help

              reciprocal = reciprocal + help
              lastshell = lastshell + help/vol       
            END IF
          END DO
        END DO
      END DO
      nrezi = nrezi + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -1, 'tolerance in long_pot not reached in reciprocal space')
    END IF

    reciprocal=(4._dp*Pi*reciprocal)/vol


    !evaluate  real space term (sum over R)   
    !/* sum over R until tolerance is reached */
    rcspace = 0._dp
    nreal = 0
    lastshell = tol+1e-8  
    DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) &
      &            .or. (abs(lastshell) .gt.  tol)))
      lastshell = 0._dp
      DO i=-nreal,nreal
        DO j=-nreal,nreal
          DO k=-nreal,nreal
            !/*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or. &
              &         (nreal.eq. abs(k)) ) THEN
              rh(1)=r(1)-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1))
              rh(2)=r(2)-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2))
              rh(3)=r(3)-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3))
              norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3))
              IF (norm .gt. 1.d-20) THEN 
                !erfc=1-erf   
                help   = terfc(alpha*norm)/norm
                rcspace = rcspace + help
                lastshell = lastshell + help
              ELSE 
                lastshell = tol+1e-8  
              END IF
            END IF
          END DO
        END DO
      END DO
      nreal = nreal + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -2, 'tolerance in long_pot not reached in real space')
    END IF


    !evaluate constant term pi/(Omega*alpha^2) 
    cterm = -Pi/(vol*alpha*alpha)

    !if r = 0 there is another constant to be added
    IF ((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) .lt. 1.0-20_dp) THEN
      cterm = cterm -2._dp*alpha/sqrt(Pi)
    END IF

    !set the value of the potential
    potential = reciprocal + rcspace + cterm

  END subroutine long_pot
  !----------------------------------------------------------------------

  
!  =====================================================================
!  evaluation of the potential phi 
!
!   phi(r) = 4*pi/Omega ( Summe{G neq 0} e^{-G^2/{4 alpha^2}}/{G^2} cos(G r)
!            +Summe{R, R neq r} (1-erf(alpha*|R-r|))/|R-r|
!            -pi/(Omega*alpha^2)
!
!
!   INPUT Parameter:
!   real(kind=dp) ::  r(3)           position of evaluation of potential
!   real(kind=dp) ::  basis(3,3)     basis of cell
!   real(kind=dp) ::  recbasis(3,3)      basis of reciprocal cell
!   real(kind=dp) ::  alpha          convergence parameter
!   real(kind=dp) ::  vol            cell volume
!   real(kind=dp) ::  tol            tolerance for convergence of "last shell"
!                       (can often be used as a criterion for global
!                        convergence) 
!   OUTPUT:
!   real(kind=dp) ::  potential      value of Ewald potential
!  
!  ======================================================================

  subroutine phi(r,basis,recbasis,alpha,vol,tol,potential, iError)

    real(kind=dp) ::  r(3), basis(3,3), recbasis(3,3), alpha, vol, tol
    real(kind=dp) ::  potential
    integer, intent(out), optional :: iError

    real(kind=dp) ::  reciprocal,rcspace,cterm
    real(kind=dp) ::  G(3),rh(3),help,norm,lastshell
    integer nrezi, nreal, nmax, nmin
    integer i,j,k

    nmax = 20
    nmin = 2

    if (present(iError)) then
      iError = 0
    end if

    !  evaluate reciprocal space term ( sum over G <> 0) ...  
    !   /* sum over G until tolerance is reached */
    nrezi = 1
    lastshell = tol+1d-8  
    reciprocal = 0._dp
    DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or. &
      &    (abs(lastshell) .gt.  tol)))
      lastshell = 0._dp
      DO i=-nrezi,nrezi
        DO j=-nrezi,nrezi
          DO k=-nrezi,nrezi
            !  /*only G belonging to outer shells are new ones */
            IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or. &
              &        (nrezi.eq. abs(k)) ) THEN
              G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1)
              G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2)
              G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3)

              help = G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
              help = exp(-help/(4._dp*alpha*alpha))/help
              help = cos(G(1)*r(1)+G(2)*r(2)+G(3)*r(3)) * help

              reciprocal = reciprocal + help
              lastshell = lastshell + help/vol       
            END IF
          END DO
        END DO
      END DO
      nrezi = nrezi + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -1, "tolerance in phi not reached in reciprocal space")
    END IF

    reciprocal=(4._dp*Pi*reciprocal)/vol


    !    evaluate  real space term (sum over R)   
    !   /* sum over R until tolerance is reached */
    rcspace = 0._dp
    nreal = 0
    lastshell = tol+1e-8  
    DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) &
      &            .or. (abs(lastshell) .gt.  tol)))
      lastshell = 0._dp
      DO i=-nreal,nreal
        DO j=-nreal,nreal
          DO k=-nreal,nreal
            ! /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or. &
              &         (nreal.eq. abs(k)) ) THEN
              rh(1)=r(1)-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1))
              rh(2)=r(2)-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2))
              rh(3)=r(3)-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3))
              norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3))
              IF (norm .gt. 1.d-20) THEN 
                !erfc=1-erf   
                help   = terfc(alpha*norm)/norm
                rcspace = rcspace + help
                lastshell = lastshell + help
              ELSE 
                lastshell = tol+1e-8  
              END IF
            END IF
          END DO
        END DO
      END DO
      nreal = nreal + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -2, "tolerance in phi not reached in real space")
    END IF


    !  evaluate constant term pi/(Omega*alpha^2) 
    cterm = -Pi/(vol*alpha*alpha)

    !  if r = 0 there is another constant to be added
    IF ((r(1)*r(1)+r(2)*r(2)+r(3)*r(3)) .lt. 1.d-20) THEN
      cterm = cterm -2._dp*alpha/sqrt(Pi)
    END IF

    !  set the value of the potential
    potential = reciprocal + rcspace + cterm

  END subroutine phi



  !  =====================================================================
  !  evaluation of the derivative of the potential phi 
  !
  !   INPUT Parameter:
  !   real(kind=dp) ::  r(3)           position of evaluation of potential
  !   real(kind=dp) ::  basis(3,3)     basis of cell
  !   real(kind=dp) ::  recbasis(3,3)      basis of reciprocal cell
  !   real(kind=dp) ::  alpha          convergence parameter
  !   real(kind=dp) ::  vol            cell volume
  !   real(kind=dp) ::  tol            tolerance for convergence of "last shell"
  !                       (can often be used as a criterion for global
  !                        convergence) 
  !   OUTPUT:
  !   real(kind=dp) ::  deriv(3)       derivative of ewlad potential
  !  
  !  ======================================================================


  subroutine phi1(r,basis,recbasis, alpha,vol,tol,deriv, iError)

    real(kind=dp) ::  r(3), basis(3,3), recbasis(3,3), alpha, vol, deriv(3)
    integer, intent(out), optional :: iError

    real(kind=dp) ::  reciprocal(3),rcspace(3)
    real(kind=dp) ::  G(3),rh(3),norm,help,tol,lastshell 
    integer i,j,k, nrezi, nreal, nmax, nmin 
    nmax = 20
    nmin = 2

    IF((r(1)*r(1) + r(2)*r(2) + r(3)*r(3)) .lt. 1.d-20) THEN
      deriv(1) = 0._dp
      deriv(2) = 0._dp
      deriv(3) = 0._dp
      return
    END IF

    !       /* evaluate reciprocal space term (sum over G <> 0) ...  */
    nrezi = 1
    lastshell = tol+1.0E-8_dp
    reciprocal(1) = 0._dp
    reciprocal(2) = 0._dp
    reciprocal(3) = 0._dp
    DO WHILE ((nrezi .le. nmax) .and. ((nrezi .le. nmin) .or. &
      &    (abs(lastshell) .gt.  tol)))
      lastshell = 0._dp
      DO i=-nrezi,nrezi
        DO j=-nrezi,nrezi
          DO k=-nrezi,nrezi
            !             /*only G belonging to outer shells are new ones */
            IF((nrezi .eq. abs(i)) .or. (nrezi .eq. abs(j)) .or. &
              &        (nrezi.eq. abs(k)) ) THEN

              G(1)=i*recbasis(1,1)+j*recbasis(2,1)+k*recbasis(3,1) 
              G(2)=i*recbasis(1,2)+j*recbasis(2,2)+k*recbasis(3,2) 
              G(3)=i*recbasis(1,3)+j*recbasis(2,3)+k*recbasis(3,3) 

              help=G(1)*G(1)+G(2)*G(2)+G(3)*G(3) 

              help=exp(-help/(4._dp*alpha*alpha))/help

              help=-sin(G(1)*r(1)+G(2)*r(2)+G(3)*r(3))*help 

              reciprocal(1)=help*G(1) + reciprocal(1)
              reciprocal(2)=help*G(2) + reciprocal(2)
              reciprocal(3)=help*G(3) + reciprocal(3)

              lastshell = lastshell + help/vol
            ENDIF
          END DO
        END DO
      END DO
      nrezi = nrezi + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -1, "tolerance in phi1 not reached in reciprocal space")
    END IF

    reciprocal(1)=(4._dp*Pi*reciprocal(1))/vol
    reciprocal(2)=(4._dp*Pi*reciprocal(2))/vol
    reciprocal(3)=(4._dp*Pi*reciprocal(3))/vol


    !       /* evaluate  real space term (sum over R) */
    !       /* sum over R until tolerance is reached */
    rcspace(1) = 0._dp
    rcspace(2) = 0._dp
    rcspace(3) = 0._dp
    nreal = 0
    lastshell = tol+1e-8
    DO WHILE ((nreal .le. nmax) .and. ((nreal .le. nmin) &
      &            .or. (abs(lastshell) .gt.  tol)))
      lastshell = 0._dp
      DO i=-nreal,nreal
        DO j=-nreal,nreal
          DO k=-nreal,nreal
            !            /*only R belonging to outer shells are new ones */
            IF((nreal .eq. abs(i)) .or. (nreal .eq. abs(j)) .or. &
              &         (nreal.eq. abs(k)) ) THEN
              rh(1)=r(1)-(i*basis(1,1)+j*basis(2,1)+k*basis(3,1)) 
              rh(2)=r(2)-(i*basis(1,2)+j*basis(2,2)+k*basis(3,2)) 
              rh(3)=r(3)-(i*basis(1,3)+j*basis(2,3)+k*basis(3,3)) 

              norm=sqrt(rh(1)*rh(1)+rh(2)*rh(2)+rh(3)*rh(3)) 

              help = (-2/sqrt(Pi)*exp(-alpha*alpha*norm*norm)* &
                  &        alpha*norm - terfc(alpha*norm))/(norm*norm*norm)
              
              rcspace(1) = rh(1)*help + rcspace(1)  
              rcspace(2) = rh(2)*help + rcspace(2)
              rcspace(3) = rh(3)*help + rcspace(3)

              lastshell = lastshell + help
            END IF
          END DO
        END DO
      END DO
      nreal = nreal + 1
    END DO

    ! stop if tolerance not reached
    IF ( abs(lastshell)  .gt. tol ) THEN
      @:ERROR_HANDLING(iError, -2, "tolerance in phi1 not reached in real space")
    END IF


    !       /* add real and reciprocal parts */
    deriv(1) = rcspace(1)  + reciprocal(1) 
    deriv(2) = rcspace(2)  + reciprocal(2) 
    deriv(3) = rcspace(3)  + reciprocal(3) 

  END subroutine phi1



  !==============================================================================
  !
  !       evaluate the cross product of A x B
  !
  !==============================================================================

  subroutine CROSS( A, B, C) 

    real(kind=dp) ::  A(3), B(3), C(3)

    C(1)=A(2)*B(3)-A(3)*B(2)
    C(2)=A(3)*B(1)-A(1)*B(3)
    C(3)=A(1)*B(2)-A(2)*B(1)
  END subroutine CROSS


  !       get reciprocal lattice vectors and volume of unit cell       

  subroutine REZVOL(basis,recbasis,vol)
    

    real(kind=dp) ::   basis(3,3), recbasis(3,3), vol
    real(kind=dp) ::   hv1(3), hv2(3), hv3(3), hv4(3), fac
    integer i
    

    DO i=1,3, 1
      hv1(i)=basis(1,i)
      hv2(i)=basis(2,i)
      hv3(i)=basis(3,i)
    END DO

    call CROSS(hv2,hv3,hv4)
    vol = abs(basis(1,1)*hv4(1)+basis(1,2)*hv4(2)+basis(1,3)*hv4(3))
    fac = 2._dp*Pi/vol

    recbasis(1,1)=hv4(1)*fac
    recbasis(1,2)=hv4(2)*fac
    recbasis(1,3)=hv4(3)*fac

    call CROSS(hv3,hv1,hv4)
    recbasis(2,1)=hv4(1)*fac 
    recbasis(2,2)=hv4(2)*fac 
    recbasis(2,3)=hv4(3)*fac

    call CROSS(hv1,hv2,hv4)
    recbasis(3,1)=hv4(1)*fac
    recbasis(3,2)=hv4(2)*fac
    recbasis(3,3)=hv4(3)*fac

  END subroutine REZVOL



  !==============================================================================
  !
  !     Returns the (tabulated) complementary error function erfc(x) 
  !     with fractional error everywhere less than 1.2 x 10^-7
  !
  !==============================================================================

  FUNCTION terfc(x)
    real(kind=dp) ::  terfc,x
    real(kind=dp) ::  t,z

    z = abs(x)
    t = 1._dp/(1._dp+0.5_dp*z)
    terfc = t*exp(-z*z-1.26551223_dp+t*(1.00002368_dp+t*(0.37409196_dp+&
    &       t*(0.09678418_dp+t*(-0.18628806_dp+t*(0.27886807_dp+&
    &       t*(-1.13520398_dp+t*(1.48851587+t*(-0.82215223+&
    &       t*0.17087277)))))))))
    if (x .lt. 0._dp) terfc=2._dp-terfc

    RETURN
  END FUNCTION terfc


  !==============================================================================
  !      get optimal alpha for Ewald potential phi
  !
  !      INPUT:
  !      real(kind=dp) ::  basis(3,3)     basis of lattice                     
  !
  !      RETURNS:
  !      real(kind=dp) ::                 optimal alpha
  !==============================================================================

  function getalpha(basis)

    real(kind=dp) ::  basis(3,3)
    real(kind=dp) ::  getalpha
    real(kind=dp) ::  alpha, alphal, alphar
    integer nopt
    real(kind=dp) ::  recbasis(3,3), vol, tol
    real(kind=dp) ::  G, R, help1, help2, help3

    tol  = 1e-5

    !       get reciprocal lattice vectors and cell volume

    CALL REZVOL(basis,recbasis,vol)

    !       get sqnorm of smallest vector in reciprocal space 
    help1 = recbasis(1,1)**2+recbasis(1,2)**2+recbasis(1,3)**2
    help2 = recbasis(2,1)**2+recbasis(2,2)**2+recbasis(2,3)**2
    help3 = recbasis(3,1)**2+recbasis(3,2)**2+recbasis(3,3)**2
    G    = sqrt(min(help1,help2,help3))

    !       get norm of smallest vector in real space 
    help1 = basis(1,1)**2 + basis(1,2)**2 + basis(1,3)**2
    help2 = basis(2,1)**2 + basis(2,2)**2 + basis(2,3)**2
    help3 = basis(3,1)**2 + basis(3,2)**2 + basis(3,3)**2
    R    = sqrt(min(help1,help2,help3))


    !       optimise alpha

    !       if (reciprocalspace decline - realspace decline) < 0 convergence
    !       in real space too slow: increase alpha
    !       set starting alphal
    alpha = 1e-5
    DO WHILE( diffrecreal(alpha,G,R,vol) .lt. tol )
      alphal = alpha
      alpha = alpha*2._dp
    END DO

    !       if (reciprocalspace decline - realspace decline) > 0 convergence
    !       in reciprocal space too slow: decrease alpha
    !       set starting alphar
    alpha = 1e+5
    DO WHILE( diffrecreal(alpha,G,R,vol) .gt. tol )
      alphar = alpha
      alpha = alpha / 2._dp
    END DO

    !       now find best value by refining the interval [alphal,alphar]
    alpha = (alphal + alphar)/2._dp
    nopt  = 0
    DO WHILE ( abs(diffrecreal(alpha,G,R,vol)) .gt. tol .and. &
      &             nopt .le. 20)
      IF ( diffrecreal(alpha,G,R,vol) .lt. tol ) THEN
        alphal = alpha
      END IF

      IF ( diffrecreal(alpha,G,R,vol) .gt. tol ) THEN
        alphar = alpha
      END IF
      alpha = (alphal + alphar)/2._dp
      nopt  = nopt + 1
    END DO

    IF (nopt .gt. 20 ) THEN 
      alpha=exp(-0.310104*log(vol)+0.786382)/2._dp
      PRINT*, "WARNING: NO OPTIMISED ALPHA FOUND: "
      PRINT*, "STANDARD ALPHA USED. ALPHA SET TO", alpha 
    END IF

    getalpha = alpha

  END function getalpha


  !==============================================================================
  !      
  !      This subroutine returns the norm of the largest reciprocal space
  !      vector and the norm of the largest real space vector needed for 
  !      a converged Ewald summation
  !
  !      INPUT:
  !      real(kind=dp) ::   alpha         chosen convergence parameter
  !      real(kind=dp) ::   tol           disired tolerance of last contribution
  !      real(kind=dp) ::   vol           cell volume
  !
  !      OUTPUT:
  !      real(kind=dp) ::   Gmax          norm of largest vector in reciprocal space
  !      real(kind=dp) ::   Rmax          norm of largest vector in real space 
  !==============================================================================


  subroutine getGRmax(alpha, tol, vol, Gmax, Rmax, iError)

    real(kind=dp) ::  alpha, tol, vol, Gmax, Rmax
    integer, intent(out), optional :: iError

    integer :: nopt
    real(kind=dp) ::  G, R, Gl, Gr, Rl, Rr


    !       set starting Gl and Rl
    G = 1.0e-5_dp
    R = 1.0e-5_dp
    DO WHILE( Gspace(G,alpha,vol) .gt. tol + 1.0e-10_dp)
      Gl = G
      G = G * 2.0_dp
    END DO

    DO WHILE( Rspace(R,alpha) .gt. tol + 1.0e-10_dp)
      Rl = R
      R = R * 2.0_dp
    END DO

    !       set starting Gr and Rr
    G = 1e+5
    R = 1e+5
    DO WHILE( Gspace(G,alpha,vol) .lt. tol - 1e-10)
      Gr = G
      G    = G/2._dp
    END DO

    DO WHILE( Rspace(R,alpha) .lt. tol - 1e-10)
      Rr = R
      R    = R/2._dp
    END DO



    ! now find best value of G by refining the interval [Gl,Gr]
    G = (Gl + Gr)/2._dp
    nopt  = 0
    DO WHILE ( Gspace(G,alpha,vol) .gt. tol .and. nopt .le. 20)
      IF ( Gspace(G,alpha,vol) .lt. (tol - 1d-10) ) THEN
        Gr = G
      END IF

      IF ( Gspace(alpha,G,vol) .gt. (tol + 1d-10) ) THEN
        Gl = G
      END IF

      G = (Gl + Gr)/2._dp
      nopt  = nopt + 1
    END DO

    IF (nopt .ge. 20) THEN
      @:ERROR_HANDLING(iError, -1, "Gmax could not be determined in getGRmax")
    END IF

    !       now find best value of R by refining the interval [Rl,Rr]
    R = (Rl + Rr)/2._dp
    nopt  = 0
    DO WHILE ( Rspace(R,alpha) .gt. tol .and. nopt .le. 20)

      IF ( Rspace(R,alpha) .lt. (tol - 1e-10) ) THEN
        Rr = R
      END IF

      IF ( Rspace(alpha,R) .gt. (tol + 1e-10) ) THEN
        Rl = R
      END IF

      R = (Rl + Rr)/2._dp
      nopt  = nopt + 1
    END DO

    IF (nopt .ge. 20) THEN
      @:ERROR_HANDLING(iError, -2, "Rmax could not be determined in getGRmax")
    END IF

    Gmax = G
    Rmax = R
    PRINT*,"GMAX:", GMAX, "REC:", Gspace(Gmax,alpha,vol)
    PRINT*,"RMAX:", RMAX, "real*8:", Rspace(Rmax,alpha)

  END subroutine getGRmax



  !==============================================================================
  !      get differnce between decline in reciprocal and real space
  !      this function is only used by function getalpha
  !
  !      INPUT:
  !      real(kind=dp) ::   alpha         convergence parameter
  !      real(kind=dp) ::   G             square of norm of smallest G
  !      real(kind=dp) ::   R             norm of smallest R
  !      real(kind=dp) ::   vol           cell volume
  !
  !      RETURNS:
  !      real(kind=dp) ::                 difference between decline in reciprocal 
  !                          space (rec(2G)-rec(3G)) and real space (real(2R) 
  !                          - real(3R))
  !==============================================================================

  function diffrecreal(alpha,G,R,vol)

    real(kind=dp) ::  alpha, G, R, vol
    real(kind=dp) ::  diffrecreal
    real(kind=dp) ::  diffrec, diffreal

    !       make differences between decline at 2G and 3G / 2R and 3R
    diffrec = Gspace(2._dp*G,alpha,vol) - Gspace(3._dp*G,alpha,vol)
    diffreal= Rspace(2._dp*R,alpha) - Rspace(3._dp*R,alpha)

    !       return difference between reciprocal and realspace decline
    diffrecreal        = diffrec - diffreal

  END function diffrecreal


  !==============================================================================
  !       returns the "r independent" G space part of the Ewald sum
  !      
  !       INPUT:
  !       real(kind=dp) ::  G       norm of G
  !       real(kind=dp) ::  alpha   chosen convergence parameter
  !       real(kind=dp) ::  vol     cell volume
  !
  !       RETURNS:     
  !       real(kind=dp) ::          "r independent" G space part of the Ewald sum
  !============================================================================== 

  real(kind=dp) function Gspace(G,alpha,vol)

    real(kind=dp) ::  G, alpha, vol
    
    !       evaluate reciprocal space term at G
    Gspace = exp(-G**2/(4._dp*alpha*alpha))/(G**2)
    Gspace = (4._dp*Pi*Gspace)/vol

  END function Gspace


  !==============================================================================
  !       returns the R space part of the Ewald sum
  !      
  !       INPUT:
  !       real(kind=dp) ::  R       norm of R
  !       real(kind=dp) ::  alpha   chosen convergence parameter
  !
  !       RETURNS:     
  !       real(kind=dp) ::          R space part of the Ewald sum
  !============================================================================== 

  real(kind=dp) function Rspace(R,alpha)

    real(kind=dp) ::  R, alpha
 
    !       evaluate real space term at R
    Rspace  = terfc(alpha*R)/R

  END function Rspace
 

end module gewald
