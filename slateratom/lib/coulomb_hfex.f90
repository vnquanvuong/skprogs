!> Module that builds up Coulomb and Hartree-Fock exchange supermatrix elements.
module coulomb_hfex

  use common_accuracy, only : dp
  use common_anglib, only : realGaunt
  use common_poisson, only : becke_grid_params, becke_integrator, integrator_init,&
      & integrator_set_kernel_param, integrator_precomp_fdmat, integrator_build_LU,&
      & integrator_get_coords, integrator_solve_helmholz
  use utilities, only : fak
  use core_overlap, only : v

  implicit none
  private

  public :: coulomb, hfex, hfex_lr


  !>
  type TSlaterBasisIndex

    real(dp) :: alpha
    integer  :: poly_order
    integer  :: ll

  end type TSlaterBasisIndex


contains

  !> Calculates Coulomb supermatrix elements,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 6 and eqn. 21
  pure subroutine coulomb(jj, max_l, num_alpha, alpha, poly_order, uu, ss)

    !> coulomb supermatrix
    real(dp), intent(out) :: jj(0:,:,:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> nucleus-electron supervector
    real(dp), intent(out) :: uu(0:,:,:)

    !> overlap supervector
    real(dp), intent(out) :: ss(0:,:,:)

    !> temporary storage
    real(dp) :: alpha1, alpha2

    !> auxiliary variables
    integer :: ii, jjj, kk, ll, mm, nn, oo, pp, qq, rr, sss, tt, uuu, vv
    integer :: nlpq, nmrs

    jj(:,:,:,:,:,:) = 0.0_dp

    do ii = 0, max_l
      sss = 0
      do jjj = 1, num_alpha(ii)
        do kk = 1, poly_order(ii)
          sss = sss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uuu = 0
                do oo = 1, num_alpha(nn)
                  do pp = 1, poly_order(nn)
                    uuu = uuu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        alpha1 = (alpha(ii, jjj) + alpha(ii, ll)) /&
                            &(alpha(nn, oo) + alpha(nn, qq))
                        alpha2 = (alpha(nn, oo) + alpha(nn, qq)) /&
                            &(alpha(ii, jjj) + alpha(ii, ll))
                        nlpq = kk + mm + 2 * ii
                        nmrs = pp + rr + 2 * nn

                        jj(ii, sss, tt, nn, uuu, vv) =&
                            & uu(ii, sss, tt) * ss(nn, uuu, vv) *&
                            & cc(nlpq - 1, nmrs, alpha1) +&
                            & uu(nn, uuu, vv) * ss(ii, sss, tt) *&
                            & cc(nmrs - 1, nlpq, alpha2)

                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine coulomb


  !> Builds HF exchange supermatrix,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 7/8 and eqn. 21
  pure subroutine hfex(kk, max_l, num_alpha, alpha, poly_order, problemsize)

    !> Hartree-Fock exchange supermatrix
    real(dp), intent(out) :: kk(0:,:,:,0:,:,:)

    !> maximum angular momentum
    integer, intent(in) :: max_l

    !> number of exponents in each shell
    integer, intent(in) :: num_alpha(0:)

    !> basis exponents
    real(dp), intent(in) :: alpha(0:,:)

    !> highest polynomial order + l in each shell
    integer, intent(in) :: poly_order(0:)

    !> maximum size of the eigenproblem
    integer, intent(in) :: problemsize

    !> Rev. Mod. Phys. 32, 186 (1960) eqn. 21
    real(dp), allocatable :: knu(:,:,:,:,:,:,:)

    !> temporary storage
    real(dp) :: alpha1, alpha2, alpha3, alpha4, beta1, beta2, beta3, beta4

    !> four terms of Rev. Mod. Phys. 32, 186 (1960) eqn. 21
    real(dp) :: pre, t1, t2, t3, t4

    !> auxiliary variables
    integer :: ii, jj, kkk, ll, mm, nn, oo, pp, qq, rr, ss, tt, uu, vv
    integer :: nu, nlp, nlq, nmr, nms

    allocate(knu(0:max_l, problemsize, problemsize, 0:max_l, problemsize, problemsize,&
        & 0:2 * max_l + 2))

    kk(:,:,:,:,:,:) = 0.0_dp
    knu(:,:,:,:,:,:,:) = 0.0_dp

    ! build knu according to Rev. Mod. Phys. 32, 186 (1960) eqn. 8

    do ii = 0, max_l
      ss = 0
      do jj = 1, num_alpha(ii)
        do kkk = 1, poly_order(ii)
          ss = ss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uu = 0
                do oo = 1, num_alpha(nn)
                  do pp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        alpha1 = 0.5_dp * (alpha(ii, jj) + alpha(nn, oo))
                        alpha2 = 0.5_dp * (alpha(ii, ll) + alpha(nn, qq))
                        alpha3 = 0.5_dp * (alpha(ii, jj) + alpha(nn, qq))
                        alpha4 = 0.5_dp * (alpha(ii, ll) + alpha(nn, oo))
                        beta1 = alpha1 / alpha2
                        beta2 = alpha2 / alpha1
                        beta3 = alpha3 / alpha4
                        beta4 = alpha4 / alpha3
                        nlp = kkk + ii
                        nlq = mm + ii
                        nmr = pp + nn
                        nms = rr + nn

                        pre = 1.0_dp / sqrt(&
                            & v(alpha(ii, jj), 2 * (kkk + ii))&
                            & * v(alpha(ii, ll), 2 * (mm + ii))&
                            & * v(alpha(nn, oo), 2 * (pp + nn))&
                            & * v(alpha(nn, qq), 2 * (rr + nn)))

                        do nu = abs(ii - nn), ii + nn, 2

                          t1 = v(alpha1, nlp + nmr - nu - 1) * v(alpha2, nlq + nms + nu) *&
                              & cc(nlp + nmr - nu - 1, nlq + nms + nu, beta1)
                          t2 = v(alpha2, nlq + nms - nu - 1) * v(alpha1, nlp + nmr + nu) *&
                              & cc(nlq + nms - nu - 1, nlp + nmr + nu, beta2)
                          t3 = v(alpha3, nlp + nms - nu - 1) * v(alpha4, nlq + nmr + nu) *&
                              & cc(nlp + nms - nu - 1, nlq + nmr + nu, beta3)
                          t4 = v(alpha4, nlq + nmr - nu - 1) * v(alpha3, nlp + nms + nu) *&
                              & cc(nlq + nmr - nu - 1, nlp + nms + nu, beta4)

                          knu(ii, ss, tt, nn, uu, vv, nu) = pre * (t1 + t2 + t3 + t4)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ! build kk according to Rev. Mod. Phys. 32, 186 (1960) eqn. 7

    do ii = 0, max_l
      ss = 0
      do jj = 1, num_alpha(ii)
        do kkk = 1, poly_order(ii)
          ss = ss + 1
          tt = 0
          do ll = 1, num_alpha(ii)
            do mm = 1, poly_order(ii)
              tt = tt + 1
              do nn = 0, max_l
                uu = 0
                do oo = 1, num_alpha(nn)
                  do pp = 1, poly_order(nn)
                    uu = uu + 1
                    vv = 0
                    do qq = 1, num_alpha(nn)
                      do rr = 1, poly_order(nn)
                        vv = vv + 1

                        do nu = abs(ii - nn), ii + nn, 2

                          kk(ii, ss, tt, nn, uu, vv) = kk(ii, ss, tt, nn, uu, vv) +&
                              & almn(ii, nn, nu) * knu(ii, ss, tt, nn, uu, vv, nu)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine hfex


  !>
  subroutine hfex_lr(k, max_l, num_alpha, alpha, poly_order, problemsize, kappa, grid_params)
    real(dp), intent(out) :: k(0:,:,:,0:,:,:)
    integer, intent(in) :: max_l
    integer, intent(in) :: num_alpha(0:)
    integer, intent(in) :: poly_order(0:)
    real(dp), intent(in) :: alpha(0:4,10)
    real(dp), intent(in) :: kappa
    type(becke_grid_params), intent(in) :: grid_params

    ! local variables:
    real(dp),allocatable :: knu(:,:,:,:,:,:,:)
    real(dp) :: alpha1,alpha2,alpha3,alpha4,beta1,beta2,beta3,beta4
    real(dp) :: pre,t1,t2,t3,t4
    integer :: ii,jj,kk,ll,mm,nn,oo,pp,qq,rr,ss,tt,uu,vv,ww,xx,yy,zz
    integer :: nu,problemsize
    integer :: nlp,nlq,nmr,nms

    ! instance of becke integrator
    type(becke_integrator) :: t_integ

    type(TSlaterBasisIndex) :: basis_index1, basis_index2
    type(TSlaterBasisIndex) :: basis_index3, basis_index4

    ! the inner integral
    real(dp), allocatable :: Vin(:,:,:,:)
    integer :: eta_max,sigma_max,ll_max,N_radial,numalpha
    integer :: porder, nalpha,ll1,ll2,alp,bet,sigma1,sigma2
    integer :: eta1,eta2,aa,bb,lambda,mu,maxllind,npl

    integer, allocatable :: llind(:,:), lambdamu(:,:)
    integer, allocatable :: aind(:,:), ps(:,:),ne_ind(:,:),ne(:,:)
    integer, allocatable :: sigmaind(:,:), ep_es(:,:)

    real(dp), allocatable :: Integrand(:,:,:)
    real(dp), pointer :: rr1(:)
    real(dp), allocatable :: rho_lm(:),expon(:),IK(:,:,:)

    real(dp), allocatable :: GauntFaktor(:,:,:), innerint(:), normfaktor(:,:,:)
    real(dp) :: gaunt_lm, gaunt2_lm,norm

    real(dp) :: res


    N_radial = grid_params%N_radial
    ll_max = grid_params%ll_max

    ! inititalize the becke_integrator
    call integrator_init(t_integ, grid_params)

    ! set the kernel parameter
    call integrator_set_kernel_param(t_integ, kappa)
    call integrator_precomp_fdmat(t_integ)
    call integrator_build_LU(t_integ)

    call integrator_get_coords(t_integ, [3,1,1], rr1)

    !! If poly order is not the same for all shells this will give a bug!
    porder = poly_order(0)
    nalpha = num_alpha(0)

    eta_max = 2*(porder + max_l)

    sigma_max = nalpha*(nalpha +1)/2

    allocate(Vin(ll_max,eta_max,sigma_max,N_radial))
    allocate(rho_lm(N_radial))
    allocate(expon(sigma_max))
    allocate(Integrand(eta_max,sigma_max,N_radial))

    maxllind = (max_l+1)*(max_l + 2)/2

    allocate(llind(max_l+1,max_l+1))
    allocate(lambdamu(maxllind,2))

    ss=1
    do ii=1, max_l+1
       do jj=1, ii
          llind(ii,jj)=ss
          llind(jj,ii)=ss
          lambdamu(ss,1)=ii
          lambdamu(ss,2)=jj
          ss=ss+1
       end do
    end do

    allocate(aind(porder*nalpha,porder*nalpha))
    allocate(ps(porder*nalpha*(porder*nalpha + 1)/2,2))

    ss=1
    do ii=1, porder*nalpha
       do jj=1, ii
          aind(ii,jj)=ss
          aind(jj,ii)=ss
          ps(ss,1)=ii
          ps(ss,2)=jj
          ss=ss+1
       end do
    end do

    allocate(ne_ind(nalpha,porder))
    allocate(ne(porder*nalpha,2))
    ss=1
    do ii=1, nalpha
       do jj=1, porder
          ne_ind(ii,jj)=ss
          ne(ss,1)=ii
          ne(ss,2)=jj
          ss=ss+1
       end do
    end do

    allocate(sigmaind(nalpha,nalpha))
    allocate(ep_es(sigma_max,2))
    ss=1
    do ii=1, nalpha
       do jj=1, ii
          sigmaind(ii,jj)=ss
          sigmaind(jj,ii)=ss
          ep_es(ss,1)=ii
          ep_es(ss,2)=jj
          ss=ss+1
       end do
    end do

    ! precompute the angular part
    allocate(GauntFaktor(max_l+1, max_l+1, ll_max))

    do ll1=0,max_l
       do ll2=0,max_l
          do ll=0,ll_max-1
             gaunt2_lm=0.0_dp
             do alp=-ll1, ll1
                do bet=-ll2, ll2
                   do mm=-ll,ll
                      gaunt_lm = realGaunt(ll1,alp,ll2,bet,ll,mm)
                      gaunt_lm=gaunt_lm*gaunt_lm
                      gaunt2_lm=gaunt2_lm + gaunt_lm
                   end do
                end do
             end do
             GauntFaktor(ll1+1,ll2+1,ll+1)=gaunt2_lm
          end do
       end do
    end do

    do ii=1,sigma_max
       expon(ii) = alpha(0,ep_es(ii,1))+alpha(0,ep_es(ii,2))
    end do

    ! fill in the integrands
    do ii=1, eta_max
       do kk=1,sigma_max
          Integrand(ii,kk,:)=rr1**(ii-1)*exp(-expon(kk)*rr1)
          do ll=0,ll_max-1
             rho_lm(:)=Integrand(ii,kk,:)
             call integrator_solve_helmholz(t_integ,ll,rho_lm)
             Vin(ll+1,ii,kk,:)=rho_lm(:)/rr1(:)*t_integ%integration_grid(3)%weight(:)
          end do
       end do
    end do

    allocate(innerint(N_radial))
    allocate(IK(2*(max_l+1)*(max_l+2)-1, porder*nalpha*(porder*nalpha + 1)/2,&
        & porder*nalpha*(porder*nalpha + 1)/2))
    ! perform the integrals and construct the exchange supermatrix
    do ll=0, maxllind - 1 ! 2 * (max_l + 1) * (max_l + 2) - 1
       do aa=1, porder * nalpha * (porder * nalpha + 1) / 2
          do bb = 1, aa

             lambda=lambdamu(ll+1,1)-1
             mu=lambdamu(ll+1,2)-1

             pp=ps(aa,1)
             ss=ps(aa,2)
             qq=ps(bb,1)
             rr=ps(bb,2)

             eta1=ne(pp,2)+ne(ss,2)+lambda+mu-1
             sigma1=sigmaind(ne(pp,1),ne(ss,1)) !ne(pp,1)+ne(ss,1)
             eta2=ne(qq,2)+ne(rr,2)+lambda+mu-1
             sigma2=sigmaind(ne(qq,1),ne(rr,1)) !ne(qq,1)+ne(rr,1)

             innerint=0.0_dp
             do ll1=abs(lambda-mu),lambda+mu,2
                if(GauntFaktor(lambda+1,mu+1,ll1+1) .ge. 1.0e-16_dp) then
                   innerint(:) = innerint(:) + Vin(ll1+1,eta1,sigma1,:)&
                       &*GauntFaktor(lambda+1,mu+1,ll1+1)
                else
                   print *, "G=0"
                end if
             end do
             IK(ll+1,aa,bb)=sum(innerint*Integrand(eta2,sigma2,:))
             IK(ll+1,bb,aa)=IK(ll+1,aa,bb)
          end do
       end do
    end do

    allocate(normfaktor(porder,nalpha,max_l+1))
    do ii=1, porder
       do kk=1, nalpha
          do jj=1, max_l+1
             npl=ii+jj-1
             normfaktor(ii,kk,jj)=sqrt(2.0_dp*alpha(0,kk)/fak(2*npl))*(2.0_dp&
                 &*alpha(0,kk))**npl
          end do
       end do
    end do

    allocate(knu(0:max_l,problemsize,problemsize,0:max_l,problemsize,&
        &problemsize,0:2*max_l+2))

    k=0.0d0
    knu=0.0d0

    ! build kk according to Rev. Mod. Phys. 32, 186 (1960) eqn. 8

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                        vv=vv+1

                        alpha1=0.5d0*(alpha(ii,jj)+alpha(nn,oo))
                        alpha2=0.5d0*(alpha(ii,ll)+alpha(nn,qq))
                        alpha3=0.5d0*(alpha(ii,jj)+alpha(nn,qq))
                        alpha4=0.5d0*(alpha(ii,ll)+alpha(nn,oo))
                        beta1=alpha1/alpha2
                        beta2=alpha2/alpha1
                        beta3=alpha3/alpha4
                        beta4=alpha4/alpha3
                        nlp=kk+ii
                        nlq=mm+ii
                        nmr=pp+nn
                        nms=rr+nn

                        pre=1.0d0/sqrt(v(alpha(ii,jj),2*(kk+ii))*&
                            &       v(alpha(ii,ll),2*(mm+ii))*&
                            &       v(alpha(nn,oo),2*(pp+nn))*&
                            &       v(alpha(nn,qq),2*(rr+nn)))

                        do nu=abs(ii-nn),ii+nn,2

                          t1=v(alpha1,nlp+nmr-nu-1)*v(alpha2,nlq+nms+nu)*&
                              &cc(nlp+nmr-nu-1,nlq+nms+nu,beta1)
                          t2=v(alpha2,nlq+nms-nu-1)*v(alpha1,nlp+nmr+nu)*&
                              &cc(nlq+nms-nu-1,nlp+nmr+nu,beta2)
                          t3=v(alpha3,nlp+nms-nu-1)*v(alpha4,nlq+nmr+nu)*&
                              &cc(nlp+nms-nu-1,nlq+nmr+nu,beta3)
                          t4=v(alpha4,nlq+nmr-nu-1)*v(alpha3,nlp+nms+nu)*&
                              &cc(nlq+nmr-nu-1,nlp+nms+nu,beta4)

                          knu(ii,ss,tt,nn,uu,vv,nu)=pre*(t1+t2+t3+t4)

                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

    ! build kk according to Rev. Mod. Phys. 32, 186 (1960) eqn. 7

    do ii=0,max_l
      ss=0
      do jj=1,num_alpha(ii)
        do kk=1,poly_order(ii)
          ss=ss+1
          tt=0
          do ll=1,num_alpha(ii)
            do mm=1,poly_order(ii)
              tt=tt+1
              do nn=0,max_l
                uu=0
                do oo=1,num_alpha(nn)
                  do pp=1,poly_order(nn)
                    uu=uu+1
                    vv=0
                    do qq=1,num_alpha(nn)
                      do rr=1,poly_order(nn)
                         vv=vv+1

                         do nu=abs(ii-nn),ii+nn,2
                            k(ii,ss,tt,nn,uu,vv)=k(ii,ss,tt,nn,uu,vv)+&
                                &almn(ii,nn,nu)*knu(ii,ss,tt,nn,uu,vv,nu)
                         end do

                         norm=normfaktor(kk,jj,ii+1)*normfaktor(mm,ll&
                             &,ii+1)
                         norm=norm*normfaktor(pp,oo,nn+1)*normfaktor(rr&
                             &,qq,nn+1)/(real((2*ii+1)*(2*nn+1),dp))

                         k(ii,ss,tt,nn,uu,vv)=k(ii,ss,tt,nn,uu,vv)-0.5_dp*(IK(llind(ii+1,nn+1),aind(ss&
                             &,vv),aind(tt,uu))+IK(llind(ii+1,nn+1),aind(ss&
                             &,uu),aind(tt,vv)))*norm

                       end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do

  end subroutine hfex_lr


  ! !>
  ! function get_lc_contrib(t_int, basis1, basis2, basis3, basis4, kappa) result(res)
  !   type(becke_integrator), intent(inout) :: t_int
  !   type(slater_basis_index), intent(in) :: basis1, basis2
  !   type(slater_basis_index), intent(in) :: basis3, basis4
  !   real(dp), intent(in) :: kappa
  !   !
  !   real(dp), pointer :: rr1(:), theta1(:), phi1(:)
  !   ! real(dp), pointer :: rr2(:), theta2(:), phi2(:)
  !   ! real(dp), pointer :: rr3(:), theta3(:), phi3(:)
  !   integer :: ngrid, ii, kk, mm1, mm2
  !   real(dp), allocatable, target :: center1(:), center2(:), center3(:), center4(:)
  !   real(dp) :: res

  !   integer :: np1,np2,np3,np4
  !   real(dp) :: alpha_x_2,alpha1,alpha2,alpha3,alpha4
  !   real(dp) :: norm1

  !   type(realtess) :: tes1, tes2, tes3, tes4
  !   type(separable_integrand) :: int1, int2

  !   res = 0.0_dp

  !   ! get the coordinates for the inner one-center integral
  !   call integrator_get_coords(t_int, [3,1,1], rr1)

  !   ! evaluate the orbital on the grid
  !   ngrid = size(rr1)

  !   allocate(int1%radial(size(rr1)))
  !   allocate(int2%radial(size(rr1)))

  !   int1%ll1=basis1%ll
  !   int1%mm1=0
  !   int1%ll2=basis2%ll
  !   int1%mm2=0

  !   int2%ll1=basis3%ll
  !   int2%mm1=0
  !   int2%ll2=basis4%ll
  !   int2%mm2=0

  !   res = 0.0_dp

  !   ! normalization const 1
  !   alpha1 = basis1%alpha
  !   alpha_x_2=basis1%alpha*2.0_dp
  !   np1 = basis1%poly_order + basis1%ll
  !   norm1=(alpha_x_2)**np1 * sqrt(alpha_x_2/fak(2*np1))

  !   ! normalization const 2
  !   alpha2 = basis2%alpha
  !   alpha_x_2=basis2%alpha*2.0_dp
  !   np2 = basis2%poly_order + basis2%ll
  !   norm1=norm1*(alpha_x_2)**np2 * sqrt(alpha_x_2/fak(2*np2))

  !   ! normalization const 3
  !   alpha3 = basis3%alpha
  !   alpha_x_2=basis3%alpha*2.0_dp
  !   np3 = basis3%poly_order + basis3%ll
  !   norm1=norm1*(alpha_x_2)**np3 * sqrt(alpha_x_2/fak(2*np3))

  !   ! normalization const 4
  !   alpha4 = basis4%alpha
  !   alpha_x_2=basis4%alpha*2.0_dp
  !   np4 = basis4%poly_order + basis4%ll
  !   norm1=norm1*(alpha_x_2)**np4 * sqrt(alpha_x_2/fak(2*np4))

  !   int1%radial = rr1**(np1+np2-2) * exp(-(alpha1+alpha2)*rr1)
  !   int2%radial = rr1**(np3+np4-2) * exp(-(alpha3+alpha4)*rr1)

  !   res=integrator_process(t_int, int1, int2)*norm1

  ! end function get_lc_contrib


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
  pure function cc(alpha, beta, tt)

    !> indices of Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
    integer, intent(in) :: alpha, beta

    !> argument of Rev. Mod. Phys. 32, 186 (1960) eqn. 22 and eqn. 23
    real(dp), intent(in) :: tt

    !> result
    real(dp) :: cc

    !> prefactor of Rev. Mod. Phys. 32, 186 (1960) eqn. 23
    real(dp) :: factor

    !> array that holds all C_{\alpha, \beta}(t)
    real(dp), allocatable :: carray(:,:)

    !> alpha and beta indices to carry out the recursion
    integer :: aa, bb

    ! early return if index smaller than zero

    if (alpha < 0) then
      cc = 0.0_dp
      return
    end if

    if (beta < 0) then
      cc = 0.0_dp
      return
    end if

    allocate(carray(0:alpha, 0:beta))

    factor = 1.0_dp / (1.0_dp + tt)

    ! Overall this is naive, the matrix could be reused to some extent ...
    ! OTOH, the matrices are relatively small.

    ! first handle Kronecker delta, three cases
    carray(0, 0) = factor
    do aa = 1, alpha
      carray(aa, 0) = factor * (tt * carray(aa - 1, 0) + 1.0_dp)
    end do
    do bb = 1, beta
      carray(0, bb) = factor * (carray(0, bb - 1))
    end do

    ! now build up from alpha, beta = 1
    do aa = 1, alpha
      do bb = 1, beta
        carray(aa, bb) = factor * (tt * carray(aa - 1, bb) + carray(aa, bb - 1))
      end do
    end do

    cc = carray(alpha, beta)

  end function cc


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 9
  pure function aa(rho)

    !> even integer, according to Rev. Mod. Phys. 32, 186 (1960) eqn. 9
    integer, intent(in) :: rho

    !> result
    real(dp) :: aa

    aa = fak(rho) / ((fak(rho / 2))**2)

  end function aa


  !> Auxiliary function,
  !! see Rev. Mod. Phys. 32, 186 (1960) eqn. 9
  pure function almn(lambda, mu, nu)

    !> indices of Rev. Mod. Phys. 32, 186 (1960) eqn. 9
    integer, intent(in) :: lambda, mu, nu

    !> result
    real(dp) :: almn

    almn = aa(lambda + mu - nu) * aa(lambda - mu + nu) * aa(mu - lambda + nu)&
        & / (real(lambda + mu + nu + 1, dp) * aa(lambda + mu + nu))

  end function almn

end module coulomb_hfex
