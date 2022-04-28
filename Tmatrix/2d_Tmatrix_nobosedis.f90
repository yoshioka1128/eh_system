      Program Tmatrix2d
      implicit none
! test: kmax=5.0d0/a2d, dw=12.5d0*E2d/Esc
! present: kmax=10.0d0/a2d, dw=37.5d0*E2d/Esc

    ! me and mh="electron and hole mass", mr="reduced mass"        
      real(8),parameter :: me=0.0665d0,mh=0.11d0,mr=me*mh/(me+mh)
    ! mass
      real(8),parameter :: mass(2)=(/me/mr,mh/mr/)
    ! sum of mass
      real(8),parameter :: M(3)=(/2.0d0*me/mr,(me+mh)/mr,2.0d0*mh/mr/)
    ! a2d="2d exciton Bohr radius per a3d", E2d="2d exciton binding energy per E3d"
      real(8),parameter :: a2d=1.0d0/2.0d0,E2d=4.0d0
    ! E2d*Esc="energy unit"
      real(8),parameter :: Esc=2.0d0

    ! initial condition iread=0: SHF result, read=1: read previdous result
      integer(4),parameter :: iread=0

    ! temp="temperature", num="carrier density"
      real(8),parameter :: temp=0.5d0*E2d/Esc,beta=1.0d0/temp,num=0.1d0/a2d**2
    ! const="SPP parameter"
      real(8),parameter :: const=4.0d0

    ! imaginary number and circular constant
      complex(8),parameter :: im=(0.0d0,1.0d0)
      real(8),parameter :: pai=dacos(-1.0d0)

! parameter --------------------------------------------------------
    ! alpha="mix parameter"
      real(8),parameter :: alpha=1.0d0
    ! gmm0="initial broadening", cnv="order of convergence"
      real(8),parameter :: gmm0=0.1d0*E2d/Esc,cnv=1.0d-4
    ! iend="iteration limit"
      integer(4),parameter :: iend=500
    ! initial shift of chemical potential
      real(8),parameter :: dmu0=0.1d0/Esc
! ------------------------------------------------------------------


! set momentum ------------------------------------------------------
  ! integral parameter in main loop
    ! kmax="k cutoff", divk="k step size"
      real(8),parameter :: kmax=5.0d0/a2d,divk=0.1d0/a2d
    ! nk="number of k point"
      integer(4),parameter :: nk=int(kmax/divk)
    ! xk and xyk="point and weight",xp and xq="relative and center of mass momentum"
      real(8) :: xk(0:nk),xyk(0:nk),xp(0:nk),xq(0:2*nk)

  ! integral parameter in Coulomb hole
    ! kCH="k cutoff", divCH="k step size"
      real(8),parameter :: kCH=500.0d0/a2d,divCH=1.0d-4/a2d
    ! nCH="number of k point"
      integer(4),parameter :: nCH=int(kCH/divCH)
    ! xCH and yCH="point and weight"
      real(8) :: xCH(0:nCH),yCH(0:nCH)

  ! integral parameter for 0 point correction
    ! n0ptk="number of fine k points"
      integer(4),parameter :: n0ptk=500
    ! x0ptk and y0ptk="point and weight"
      real(8) :: x0ptk(n0ptk),y0ptk(n0ptk)
! ----------------------------------------------------------------------


! set omega -------------------------------------------------------
    ! divw="step size of omega", dw="extarnal part"
      real(8),parameter :: divw=1.0d0/40.0d0*E2d/Esc,dw=12.5d0*E2d/Esc
    ! wemax and whmax="highest energy of free electron and hole energy"
      real(8), parameter :: wemax=0.5d0*kmax**2/mass(1),whmax=0.5d0*kmax**2/mass(2)

  ! single particle
    ! nwn1="low energy edge point of electron or hole energy"
      integer(4),parameter :: nwn1=-int(dw/divw)
    ! nwe2 and nwh2="high energy edge point of electron and hole energy"
      integer(4),parameter :: nwe2=int((dw+wemax)/divw),nwh2=int((dw+whmax)/divw)
    ! set
      integer(4),parameter :: nwg1(2)=(/nwe2,nwh2/)

  ! two particle
    ! nwnn1="low energy edge point of pair energy"
      integer(4),parameter :: nwnn1=nwn1*2
    ! nwee2, nwhh2, and nweh2 ="high energy edge point of pair energy"
      integer(4),parameter :: nwee2=nwe2*2,nwhh2=nwh2*2,nweh2=nwe2+nwh2
    ! set
      integer(4),parameter :: nwg2(3)=(/nwee2,nweh2,nwhh2/)

    ! omega points
      real(8) :: xw(nwnn1:nwee2)
! ------------------------------------------------------------------


! set angle k1,k2 <-> p,q Transformation ---------------------------------
    ! nt0="number of angle point"
      integer(4),parameter :: nt0=60
    ! divt="step size of angle"
      real(8),parameter :: divt=pai/dble(nt0)
    ! xt0 and yt0="point and weight"
      real(8) :: xt0(0:nt0),yt0(0:nt0)
! ------------------------------------------------------------------


! set angle partial wave expansion parameter (interaction)-----------------------
    ! mmax="expansion component"
      integer(4),parameter :: mmax=3
    ! nt="number of angle point"
      integer(4) :: nt(0:mmax)
    ! xtt and ytt="point and weight"
      real(8) :: xtt(nt0*(mmax+1)),ytt(nt0*(mmax+1))
    ! xt and yt="point and weight in gauleg"
      real(8) :: xt(nt0*(mmax+1),0:mmax),yt(nt0*(mmax+1),0:mmax)
! ---------------------------------------------------------------------


! pair green's function -----------------------------------------------------
    ! gpq="pair green's function"
    ! allocate (gpq(nwnn1:nwee2,0:nk,0:2*nk,3))
      complex(8),allocatable :: gpq(:,:,:,:)
    ! Imbgpq="bare pair green's function without bose distribution function"
    ! allocate (Imbgpq(nwnn1:nwee2,0:nk,0:2*nk,3))
      real(8),allocatable :: Imbgpq(:,:,:,:)
! ---------------------------------------------------------------------------


! effective interaction ----------------------------------------------------
    ! Wscexp and V0exp="screened and bare interaction"
      real(8) :: Wscexp(0:nk,0:nk,0:mmax),V0exp(0:nk,0:nk,0:mmax)
    ! WscV0exp="Wsc*V0"
      real(8) :: WscV0exp(0:nk,0:nk,0:mmax)
! --------------------------------------------------------------------------


! Tmatrix and Montroll Ward ----------------------------------------------
    ! npnum="number of p points", IPIV="use in ZGESV"
      integer(4) :: npnum(0:2*nk,3),IPIV(nk+1)
    ! ImTpq and ImTkk="Im part of Tmatrix+MW as a function of (p,q) and (k1,k2)"
    ! allocate (ImTpq(nwnn1:nwee2,0:nk,0:2*nk,3))
    ! allocate (ImTkk(nwnn1:nwee2,0:nk,0:nk,2,2))
      real(8),allocatable :: ImTpq(:,:,:,:),ImTkk(:,:,:,:,:)
    ! ImbTpq and ImbTkk="Im part of Tmatrix+MW without bose distribution function"
    ! allocate (ImbTpq(nwnn1:nwee2,0:nk,0:2*nk,3))
    ! allocate (ImbTkk(nwnn1:nwee2,0:nk,0:nk,2,2))
      real(8),allocatable :: ImbTpq(:,:,:,:),ImbTkk(:,:,:,:,:)
    ! Wsc2nd="Wsc*V0-Wsc**2", ImMW="Montroll-Word term", ImbMW=bfnc*ImMW
      real(8) :: Wsc2nd(0:nk,0:nk,0:mmax),ImMW(nk+1),ImbMW(nk+1)
    ! matrix and T2="use in ZGESV"
      complex(8) :: matrix(nk+1,nk+1),T2(nk+1,nk+1)
    ! dmatrix and dT2="use in DGESV"
      real(8) :: dmatrix(nk+1,nk+1),dT2(nk+1,nk+1)
    ! cffT and cffTex="coefficients used in TmatrixMW"
      real(8) :: cffT(0:mmax),cffTex(0:mmax,3)
! -------------------------------------------------------------------------


! self energy -------------------------------------------------------        
    ! "exchange term or SHF self energy"
      real(8) :: self0(0:nk,2)
    ! selfT and self="correlation term and total selfenergy"
      complex(8) :: selfT(nwn1:nwe2,0:nk,2),self(nwn1:nwe2,0:nk,2),selfold(nwn1:nwe2,0:nk,2)
! --------------------------------------------------------------------


! single particle -------------------------------------------------------
    ! Dos and Akw="density of states and spectral weight"
      real(8) :: Dos(nwn1:nwe2,2),Akw(nwn1:nwe2,0:nk,2),Akwold(nwn1:nwe2,0:nk,2)
    ! engkf, engkSHF, and engkQP="free, SHF, and QP energy"
      real(8) :: engkf(0:nk,2),engkSHF(0:nk,2),engkQP(0:nk,2)
    ! fermi and bose="fermi and bosedistribution function"
      real(8) :: fermi0(0:nk,2),fermi(nwn1:nwe2,2),bose(nwnn1:nwee2,3)
    ! "pauli blocking factor"
      real(8) :: pblock(nwn1:nwe2,nwn1:nwe2,3)
    ! numeh and numehold ="carrier density", numeh0="QP carrier density"
      real(8) :: numeh(2),numehold(2),numeh0(2),mu(2)
! -----------------------------------------------------------------------


! optical response -----------------------------------------------------------
    ! g2="q=0 pair green's function"
      complex(8) :: g2(nwnn1:nweh2,0:nk)
    ! opt and opt0="full and bare optical susceptabilitity"
      complex(8) :: opt(nwnn1:nweh2),opt0(nwnn1:nweh2)
    ! bopt="opt tiems bose distribution function"
      real(8) :: bopt(nwnn1:nweh2)
    ! vertex="work array"
      complex(8) :: vertex(nk+1)
! -----------------------------------------------------------------------


! Principal Integration ---------------------------------------------------
    ! log part of principal integration
      real(8) :: logrm0(nwnn1:nwee2,nwnn1:nwee2)
    ! Sw and fnc="input and output function"
      real(8) :: Sw(nwnn1:nwee2),fnc(nwnn1:nwee2)
! --------------------------------------------------------------------------

      
! k1,k2 <-> p,q transformation ---------------------------------------------
    ! pq="from k1,k2 to p,q", kk="from p,q to k1,k2"
      real(8) :: pq(0:nt0,0:nk,0:nk,3,2),kk(0:nt0,0:nk,0:2*nk,3,2)
    ! overlap angle in angle average
      real(8) :: angleave(0:nk,0:2*nk,3)
    ! shift momentum number 
      integer(4) :: iplusk(0:nk),iplusq(0:2*nk)
    ! ioutpq="read k1,k2 output p,q",ioutkk="read p,q output k1,k2"
      integer(4) :: ioutpq(0:nt0,0:nk,0:nk,3,2),ioutkk(0:nt0,0:nk,0:2*nk,3,2)
    ! nt0bd and npbd="boundary of angle and momentum integration"
      integer(4) :: nt0bd(0:nk,0:2*nk,3,2),npbd(0:2*nk,3,2)
! --------------------------------------------------------------------------------

      integer(4) :: i,j,iw,iw1,iw2,jm,l,ieh,ieh1,ieh2,&
           ik1,ik2,ik,ip,iq,it,it0,nwmax,icount,INFO
      real(8) :: k1,k2,selfCH,ksc,ffnc,bfnc,w,ww,p,q,t0,&
           CH,abc,def,ghi,jkl,aion,screening,cffInt






      write(6,*) 'omega points',nwee2-nwnn1+1
      write(6,*) 'gpq size',&
           dble(nwee2-nwnn1+1)*dble(nk+1)*dble(2*nk+1)*3.0d0*8.0d-9,'G byte'
      write(6,*) 'Akw size',&
           dble(nwe2-nwn1+1)*dble(nk+1)*2.0d0*8.0d-9,'G byte'

      write(6,*) 'T (E2d)',Esc/beta/E2d
      write(6,*) 'carrier density (1/a2d^2)',num*a2d**2
      write(6,*)
      write(6,*) 'setup'
      write(6,*) 'kmax=',kmax*a2d
      write(6,*) 'nk=',nk
      write(6,*) 'mmax=',mmax
      write(6,*) 'dw=',dw*Esc/E2d
      write(6,*) 'divw=',divw*Esc/E2d
      write(6,*) 'alpha=',alpha
      write(6,*) 'iread=',iread
      write(6,*)








! set momentum ----------------------------------------
  ! integral parameter in main loop
    ! xk="k points"
      do ik=0,nk ! k
         xk(ik)=divk*dble(ik)
      end do
    ! xyk="weight of radial integral of k"
      xyk(0)=0.0d0
      do ik=1,nk
         xyk(ik)=divk**2*dble(ik)
      end do
    ! xp="relative momentum"
      xp=xk
    ! xq="center of mass momentum"
      do iq=0,2*nk
         xq(iq)=divk*dble(iq)
      end do

  ! integral parameter in Coulomb hole
    ! xk="k points"
      do ik=0,nCH
         xCH(ik)=divCH*dble(ik)
      end do
    ! xyk="weight of radial integral of k"
      yCH(0)=0.5d0*divCH
      do ik=1,nCH
         yCH(ik)=divCH
      end do

  ! integral parameter in 0pt correction
      call gauleg(0.0d0,divk,x0ptk,y0ptk,n0ptk)
! -----------------------------------------------------------


! set omega -----------------------------------------------
    ! xw="omega points"
      do iw=nwnn1,nwee2
         xw(iw)=divw*dble(iw)
      end do
! -----------------------------------------------------------


! set angle average parameter --------------------------------
      do it0=0,nt0
        xt0(it0)=divt*dble(it0)
      end do
      
      yt0(0)=0.5d0*divt
      do it0=1,nt0-1
        yt0(it0)=divt
      end do
      yt0(nt0)=0.5d0*divt
! ------------------------------------------------------------


! set partial wave expansion paramter (interaction) ------------
    ! jm label
      do jm=0,mmax
        nt(jm)=nt0*(jm+1) ! nt(j) 
      end do

    ! Cos expansion
      xt=0.0d0
      yt=0.0d0
      do jm=0,mmax
        call gauleg(0.0d0,pai,xtt,ytt,nt(jm)) ! numerical recipes
        do it=1,nt(jm)
          xt(it,jm)=xtt(it)
          yt(it,jm)=ytt(it)
        end do
      end do
! -----------------------------------------------------------------

      
! principal integration ------------------------------------------------
      logrm0=0.0d0
      do iw=nwnn1,nwee2 ! w
        w=xw(iw)
        do iw2=nwnn1,nwee2-1 ! iw2
          if(iw2.eq.iw-1.or.iw2.eq.iw) then
            abc=1.0d0
          else
            abc=(xw(iw2+1)-w)/(xw(iw2)-w)
          end if
          logrm0(iw2,iw)=dlog(abc)
        end do
      end do
! ------------------------------------------------------------------------


! set coefficients used in TmatrixMW -------------------------------------
    ! coefficient of Tmatrix direct term 
      cffT=2.0d0
      cffT(0)=1.0d0
    ! coefficient of Tmatrix exchange term
      cffTex=0.0d0
      do jm=1,mmax
        cffTex(jm,1)=2.0d0*(-1.0d0)**jm
        cffTex(jm,3)=2.0d0*(-1.0d0)**jm
      end do
      cffTex(0,1)=1.0d0
      cffTex(0,3)=1.0d0
! -----------------------------------------------------------------------


      
! k1,k2 <-> p,q transformation -----------------------------------------------
    ! read: k1,k2,theta_k1k2 out: p,q
      call outpq(ioutpq,pq,nt0,nk,xk,xt0,mass,M,xp,xq)
    ! read: p,q,theta_pq out: k1,k2
      call outkk(kk,nt0,nk,ioutkk,xk,xq,xp,xt0,yt0,npbd,nt0bd,npnum,angleave,mass,M)
    ! shift momentum number
      do ik=0,nk ! k ! p or k
        iplusk(ik)=ik+1
      end do
      iplusk(nk)=nk
      do iq=0,2*nk ! q
        iplusq(iq)=iq+1
      end do
      iplusq(2*nk)=2*nk
! --------------------------------------------------------------------------------

      
! free particle energy ---------------------------------------------
      do ieh=1,2 
        do ik1=0,nk ! k1
          engkf(ik1,ieh)=0.5d0*xk(ik1)**2/mass(ieh)
        end do
      end do
! --------------------------------------------------------------


! set initial condition -----------------------------------------------
      if(iread.eq.1) then 
      ! previous mu and selfenergy read write(89)
        read(5,*) mu(1),mu(2)
        write(6,*) 'previous mu Ry'
        write(6,*) Esc*mu(1)/E2d,Esc*mu(2)/E2d
        do ieh1=1,2 
          do ik1=0,nk ! k1
            do iw=nwn1,nwg1(ieh1)
              read(5,*) self(iw,ik1,ieh1)
            end do
          end do
        end do
        
      else 
      ! SHF mu and selfenergy 
      ! set initial mu
        do i=1,2
          mu(i)=dlog(dexp(beta*pai*num/mass(i))-1.0d0)/beta ! free mu
        end do
      ! SHF
        call SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
           xyk,mmax,nt,nt0,yt,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt,&
           nCH,xCH,yCH,self0,cnv,const,iend)
      ! set complex selfenergy
        do ieh1=1,2 
          do ik1=0,nk ! k1
            do iw=nwn1,nwg1(ieh1) ! w
              self(iw,ik1,ieh1)=self0(ik1,ieh1)-im*gmm0
            end do
          end do
        end do

      end if
! --------------------------------------------------------------------------


! set Akw Dos and calculate carrier density(numeh) -------------------------
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,xk,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,xyk,numeh,numeh0,&
      mu,beta)
      write(6,*) 'e-h particle density (1/a2d^2)'
      write(6,*) 0,numeh(1)*a2d**2,numeh(2)*a2d**2      
! --------------------------------------------------------------------------      


! resutlts in progress ------------------------------------------------------------
      open(29,file='T0.5n0.1_E_num_ite.txt',status='unknown')
      open(30,file='T0.5n0.1_H_num_ite.txt',status='unknown')
      open(129,file='T0.5n0.1_ksc_ite.txt',status='unknown')
      open(130,file='T0.5n0.1_mu_ite.txt',status='unknown')
      open(133,file='T0.5n0.1_alpha_ite.txt',status='unknown')
      write(29,*) 0,numeh(1)*a2d**2,num*a2d**2
      write(30,*) 0,numeh(2)*a2d**2,num*a2d**2
      write(129,*) 0,ksc*a2d
      write(130,*) 0,Esc*(mu(1)+mu(2))/E2d
      write(133,*) 0,1.0d0
! --------------------------------------------------------------------------
      

! set initial condition ---------------------------------
      selfold=self
      numehold=numeh
      Akwold=Akw
! -------------------------------------------------------


! set bare interaction -------------------------------------------
      call Intexp(mmax,mmax,nt,nk,V0exp,yt,0.0d0,0.0d0,nt0,&
      n0ptk,x0ptk,y0ptk,divk,xk,xt)
! ----------------------------------------------------------------      

    ! main loop 700
      icount=0
      do 700
    ! icount="iteration number"
      icount=icount+1
      
! Mix selfenergy by the rate, alpha -----------------------------
      do ieh=1,2 
        do ik=0,nk ! k
          do iw=nwn1,nwg1(ieh)
            self(iw,ik,ieh)=alpha*self(iw,ik,ieh)+&
            (1.0d0-alpha)*selfold(iw,ik,ieh)
          end do
        end do
      end do
      selfold=self
! --------------------------------------------------------------


! set chemical potential ------------------------------------------
    ! set Akw Dos and calculate carrier density(numeh)
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,xk,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,xyk,numeh,numeh0,&
      mu,beta)

    ! adjust mu (numeh=num)
      do ieh=1,2 
        nwmax=nwg1(ieh)
        call shiftmu(num,ieh,mu,dmu0,numeh,&
        beta,Dos,xw,divw,nwn1,nwe2,nwnn1,nwee2,nwmax)
        do iw=nwn1,nwe2
          fermi(iw,ieh)=ffnc(beta,xw(iw),mu(ieh)) ! fermi distribution
        end do
      end do
    ! pauli blocking factor
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iw2=nwn1,nwg1(ieh2)
            do iw1=nwn1,nwg1(ieh1)
               pblock(iw1,iw2,ieh)=1.0d0-fermi(iw1,ieh1)-fermi(iw2,ieh2)
            end do
          end do
        end do
      end do
    ! QP fermi distribution
      do ieh=1,2 
        do ik=0,nk ! k
          fermi0(ik,ieh)=ffnc(beta,engkQP(ik,ieh),mu(ieh))
        end do
      end do

    ! bose distribution (0point correction)
!      call bosedis(bose,mu,nwnn1,nwee2,xw,&
!      beta,divw)
! --------------------------------------------------------------

      
! update interaction -------------------------------------------
    ! screening wave number
      ksc=screening(beta,engkQP,nk,xyk,fermi0)
    ! coefficient used in Intexp and Int2nd
      cffInt=const*ksc/32.0d0/num/pai
    ! effective interaction 
      call Intexp(mmax,mmax,nt,nk,Wscexp,yt,ksc,cffInt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,xk,xt)
    ! set WscV0exp="Wsc*V0" (for Montroll Word term)
      call Int2nd(mmax,mmax,nt,nk,WscV0exp,&
        yt,ksc,cffInt,nt0,n0ptk,x0ptk,y0ptk,divk,xk,xt)
! --------------------------------------------------------------


! HF selfenergy (bare interaction) ------------------------
      self0=0.0d0
      do ieh=1,2 
        do ik2=0,nk ! k2
          abc=0.0d0
          do iw=nwn1,nwg1(ieh)
            abc=abc+fermi(iw,ieh)*Akw(iw,ik2,ieh)
          end do
          abc=xyk(ik2)*divw*abc/(2.0d0*pai)
          do ik1=0,nk ! k1
            self0(ik1,ieh)=self0(ik1,ieh)+V0exp(ik1,ik2,0)*abc
          end do
        end do
      end do
      self0=-self0/(2.0d0*pai)
! ---------------------------------------------------------------------------


! set pair Green's fnc. -----------------------------------------------
      allocate (gpq(nwnn1:nwee2,0:nk,0:2*nk,3))
      allocate (Imbgpq(nwnn1:nwee2,0:nk,0:2*nk,3))
      gpq=0.0d0
    ! set imaginary part of gpq
      call Impair(gpq,Imbgpq,nweh2,nk,nt0,nwn1,nwe2,nwnn1,nwee2,&
        nwg1,Akw,pblock,fermi,xk,xw,yt0,npbd,nt0bd,divk,iplusk,kk,ioutkk,&
        angleave,divw,nwg2)

    ! Kramers-Kronig translation
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
              do iw=nwnn1,nwg2(ieh) ! w
                fnc(iw)=dimag(gpq(iw,ip,iq,ieh))
              end do
              call principal(nwnn1,nwg2(ieh),nwnn1,nwee2,&
              xw,fnc,logrm0,Sw)
              do iw=nwnn1,nwg2(ieh)
                gpq(iw,ip,iq,ieh)=gpq(iw,ip,iq,ieh)+Sw(iw)/pai
              end do
            end do
          end do
        end do
      end do
! --------------------------------------------------------------      

        
! set Tmatrix and Montroll Word -----------------------------------------       
    ! set WscV0exp="Wsc*V0-Wsc**2"
      do jm=0,mmax
        do ik2=0,nk ! k2
          do ik1=0,nk ! k1
            Wsc2nd(ik1,ik2,jm)=&
            WscV0exp(ik1,ik2,jm)-Wscexp(ik1,ik2,jm)**2
          end do
        end do
      end do
      allocate (ImTpq(nwnn1:nwee2,0:nk,0:2*nk,3))
      allocate (ImbTpq(nwnn1:nwee2,0:nk,0:2*nk,3))
    ! imaginary part of correlation term
      call TmatrixMW(ImTpq,ImbTpq,gpq,Imbgpq,nwee2,nwnn1,&
           nk,xyk,npbd,npnum,cffT,cffTex,Wscexp,Wsc2nd,&
        mmax,nwg2,dT2,T2,dmatrix,matrix,ImMW,ImbMW,IPIV,xk,xw)
    ! set g0(p,q=0,w)
      do iw=nwnn1,nweh2
        do ik=0,nk ! k
          g2(iw,ik)=-gpq(iw,ik,0,2)
        end do
      end do
      deallocate(gpq,Imbgpq)
! --------------------------------------------------------------------------


! theta_pq integration and (p,q)->(k1,k2) -----------------------------------
      allocate (ImTkk(nwnn1:nwee2,0:nk,0:nk,2,2))
      allocate (ImbTkk(nwnn1:nwee2,0:nk,0:nk,2,2))
      call trnsTpq(ImTpq,ImbTpq,ImTkk,ImbTkk,nwnn1,nwee2,nt0,nk,nwg2,&
                   ioutpq,pq,xk,xq,yt0,iplusk,iplusq,divk)
      deallocate(ImTpq,ImbTpq)
    ! eh pair <- he pair
      do ik2=0,nk ! k2
        do ik1=0,nk ! k1
          do iw=nwnn1,nwg2(2)
            ImTkk(iw,ik1,ik2,1,2)=ImTkk(iw,ik2,ik1,2,1)
            ImbTkk(iw,ik1,ik2,1,2)=ImbTkk(iw,ik2,ik1,2,1)
          end do
        end do
      end do
! --------------------------------------------------------------------------


! Tmatrix selfenergy ------------------------------------------------------
    ! imaginary part of self energy
      call ImselfTmatrixMW(selfT,ImTkk,ImbTkk,nwn1,nwe2,xyk,&
           nwnn1,nwee2,nk,nwg1,fermi,bose,Akw,divw)
      deallocate(ImTkk,ImbTkk)
    ! check sign of Im part selfenergy 
      do ieh=1,2 
        do ik1=0,nk ! k1
          do iw1=nwn1,nwg1(ieh)
            if(dimag(selfT(iw1,ik1,ieh)).gt.0.0d0) then
              write(6,*) 'im_self is negative'
              selfT(iw1,ik1,ieh)=0.0d0
            end if
          end do
        end do
      end do
    ! Kramers-Kronig translation
      do ieh=1,2 
        do ik=0,nk ! k
          do iw=nwn1,nwg1(ieh)
            fnc(iw)=dimag(selfT(iw,ik,ieh))
          end do
          call principal(nwn1,nwg1(ieh),nwnn1,nwee2,xw,fnc,logrm0,Sw)
          do iw=nwn1,nwg1(ieh)
            selfT(iw,ik,ieh)=selfT(iw,ik,ieh)+Sw(iw)/pai
          end do
        end do
      end do
! ---------------------------------------------------------------------------

      
! total selfenergy ---------------------------------------------------------
      do ieh=1,2 
        do ik=0,nk ! k
          do iw=nwn1,nwg1(ieh)
            self(iw,ik,ieh)=selfT(iw,ik,ieh)+self0(ik,ieh)
          end do
        end do
      end do
! ---------------------------------------------------------------------------

      
! total carrier density and set Akw ------------------------------------
      call Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,xk,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,xyk,numeh,numeh0,&
      mu,beta)
! ----------------------------------------------------------------------      


! results in progress -------------------------------------------------
      call output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,ksc,Akw,numeh,numeh0,mu,self)
      aion=(numeh0(1)/numeh(1)+numeh0(2)/numeh(2))/2.0d0
      write(6,*) icount,numeh(1)*a2d**2,numeh(2)*a2d**2
      write(29,*) icount,numeh(1)*a2d**2,num*a2d**2
      write(30,*) icount,numeh(2)*a2d**2,num*a2d**2
      write(129,*) icount,ksc*a2d
      write(130,*) icount,Esc*(mu(1)+mu(2))/E2d
      write(133,*) icount,aion
! --------------------------------------------------------------

! optical response -----------------------------------------------------------
      call optical(nk,nwnn1,nwee2,nweh2,Wscexp,mmax,opt,opt0,bopt,&
        xw,xyk,g2,vertex,matrix,IPIV,beta,mu,num,divw)
! --------------------------------------------------------------------------


! convergence check --------------------------------------------------
      if(abs(numeh(1)-numehold(1))/numeh(1).lt.cnv.and.&
      abs(numeh(2)-numehold(2))/numeh(2).lt.cnv.and.&
      abs(numeh(1)-numeh(2))/numeh(1).lt.cnv.and.&
      abs(numeh(1)-num)/num.lt.cnv) then
        write(6,*) 'converged !'
        exit
      else if(icount.ge.iend) then
        write(6,*) 'convergence error'
        exit
      end if
! -------------------------------------------------------------------------

    ! set previous density
      numehold=numeh
      Akwold=Akw

700   end do
      close(29)
      close(30)
      close(129)
      close(130)
      close(133)
      
! check difference of Akw ------------------------------------------------      
      open(70,file='T0.5n0.1_convergence_check.txt',status='unknown')
      write(70,*) 'E-H num=',numeh(1),numeh(2)
      do ieh1=1,2 
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          do ik1=0,nk ! k1
            k1=xk(ik1)
            abc=Akw(iw,ik1,ieh1)-Akwold(iw,ik1,ieh1)
            if(abc/2.0d0*pai.gt.1.0d-3) then
              write(70,*) k1*a2d,abc
            end if
          end do
        end do
      end do
      close(70)
! ---------------------------------------------------------------------------    


! output optical response and results ---------------------------------------
      call optical2(nk,nwnn1,nwee2,nweh2,Wscexp,mmax,opt,opt0,bopt,&
        xw,xyk,g2,vertex,matrix,IPIV,beta,mu,ksc,cffInt,yt,xt,nt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,mr,divw,num)
      call output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,ksc,Akw,numeh,numeh0,mu,self)
! --------------------------------------------------------------------------

      
      end program Tmatrix2d
! ************************************************************************













      
      
! ***********************************************************************
      subroutine SHF(mu,beta,num,nk,engkSHF,engkf,fermi0,dmu0,&
           xyk,mmax,nt,nt0,yt,Wscexp,n0ptk,x0ptk,y0ptk,divk,xk,xt,&
           nCH,xCH,yCH,self0,cnv,const,iend)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: i,ieh,ik,ik1,ik2,nk,icount,mmax,nt,nt0,n0ptk,nCH,iend
      real(8) :: mu(2),engkSHF(0:nk,2),engkf(0:nk,2),fermi0(0:nk,2),&
           xyk(0:nk),Wscexp(0:nk,0:nk,0:mmax),yt(nt0*(mmax+1),0:mmax),&
           x0ptk(n0ptk),y0ptk(n0ptk),xk(0:nk),xt(nt0*(mmax+1),0:mmax),&
           xCH(0:nCH),yCH(0:nCH),self0(0:nk,2),numeh(2),numehold(2)
      real(8) :: beta,num,const,dmu0,f,ksc,cffInt,divk,selfCH,CH,&
           abc,k1,cnv,screening,ffnc

    ! initial energy and fermi distribution
      do ieh=1,2 
        do ik=0,nk ! k
          engkSHF(ik,ieh)=engkf(ik,ieh)
          fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
        end do
      end do

      write(6,*) 'SHF e-h particle density (1/a2d^2)'
      icount=0
! SHF iteration -------------------------------------------------
      do 600
        icount=icount+1
      ! set mueh
        do ieh=1,2 
          call shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk,engkSHF)
          do ik=0,nk ! k
            fermi0(ik,ieh)=ffnc(beta,engkSHF(ik,ieh),mu(ieh))
          end do
        end do

      ! effective interaction
        ksc=screening(beta,engkSHF,nk,xyk,fermi0)
        cffInt=const*ksc/32.0d0/num/pai
        call Intexp(0,mmax,nt,nk,Wscexp,yt,ksc,cffInt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,xk,xt)        

      ! Coulomb hole enery
        selfCH=CH(ksc,nCH,xCH,yCH,cffInt)

      ! exchange energy
        self0=0.0d0
        do ieh=1,2 
           do ik2=0,nk ! k2
              abc=xyk(ik2)*fermi0(ik2,ieh)
              do ik1=0,nk ! k1
                 self0(ik1,ieh)=self0(ik1,ieh)+Wscexp(ik1,ik2,0)*abc
              end do
           end do
        end do
        self0=-self0/(2.0d0*pai)+selfCH

      ! total energy
        do ieh=1,2 
          do ik=0,nk ! k
            k1=xk(ik)
            engkSHF(ik,ieh)=engkf(ik,ieh)+self0(ik,ieh)
          end do
        end do

      ! total density
        numeh=0.0d0
        do ieh=1,2 
           abc=0.0d0
           do ik1=0,nk ! k1
              abc=abc+xyk(ik1)*ffnc(beta,engkSHF(ik1,ieh),mu(ieh))
           end do
           numeh(ieh)=abc
        end do
        numeh=2.0d0*numeh/(2.0d0*pai)

        write(6,*) icount,numeh(1)*a2d**2,numeh(2)*a2d**2

      ! convergence check
        if(abs(numeh(1)-numehold(1))/numeh(1).lt.cnv.and.&
             abs(numeh(2)-numehold(2))/numeh(2).lt.cnv.and.&
             abs(numeh(1)-numeh(2))/numeh(1).lt.cnv.and.&
             abs(numeh(1)-num)/num.lt.cnv) exit

        
        if(icount.eq.iend) then
          write(6,*) 'convergence error'
          stop
        end if
        
        numehold=numeh

600   end do
!-------------------------------------------------------------------

! output SHF results -----------------------------------------------------
      write(6,*) 'SHF results'
      write(6,*) 'total mu (E2d)'
      write(6,*) Esc*(mu(1)+mu(2))/E2d
      write(6,*) 'BGR (E2d)'
      write(6,*) (self0(0,1)+self0(0,2))*Esc/E2d
      write(6,*) 'Coulomb Hole (E2d)'
      write(6,*) selfCH*Esc/E2d

      open(27,file='T0.5n0.1_SHF_results.txt',status='unknown')
      write(27,*) 'screening wave number',ksc*a2d
      write(27,*) 'Coulomb hole self energy',selfCH*Esc/E2d
      write(27,*) 'E-H screened exchange self energy',&
           (self0(0,1)-selfCH)*Esc/E2d,(self0(0,2)-selfCH)*Esc/E2d
      write(27,*) 'BGR (E2d)',&
           (self0(0,1)+self0(0,2))*Esc/E2d
! -----------------------------------------------------------------

      end subroutine SHF
! ***********************************************************************
      
      
      
      
      
      
      
      
      
      
      
! ***********************************************************************
      subroutine shiftmuSHF(mu,ieh,dmu0,num,beta,nk,xyk,engk)
      implicit none
      integer(4) :: i,j,sign,nk,ik1,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-7
      real(8) :: mu(2),mu2,num,dmu0,beta,mumax,mumin,dense,ffnc
      real(8) :: xyk(0:nk),engk(0:nk,2)
      
! electron ------------------------------------------------------------
! check upper limit
    ! SHF density
      dense=0.0d0
      do ik1=0,nk ! k1
        dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu(ieh))
      end do
      dense=2.0d0*dense/(2.0d0*pai)

      if(dense.gt.num) then
        sign=-1
        mumax=mu(ieh)
      else
        sign=1
        mumin=mu(ieh)
      end if
! narrow down the chemical potential
      i=0
      do
        i=i+sign
        mu2=mu(ieh)+dble(i)*dmu0
        
      ! SHF density
        dense=0.0d0
        do ik1=0,nk ! k1
          dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai)
        if(dense.lt.num) then
          mumin=mu2
          if(sign.eq.-1) exit
        else
          mumax=mu2
          if(sign.eq.1) exit
        end if
      end do
      
      i=0
      do
        i=i+1
        mu2=(mumax+mumin)/2.0d0

      ! SHF density  
        dense=0.0d0
        do ik1=0,nk ! k1
          dense=dense+xyk(ik1)*ffnc(beta,engk(ik1,ieh),mu2)
        end do
        dense=2.0d0*dense/(2.0d0*pai)

        if(abs(dense-num)/num.lt.hdig) then
!          write(6,*) (mumax+mumin)/2.0d0
          exit
        else if(dense.lt.num) then
          mumin=mu2
        else
          mumax=mu2
        end if
        if(i.gt.100) then
! chemical potential located in the energy gap of finite size mesh pt.
          write(6,*) 'step size of wavenumber (1/nk) is too large'
          stop
        end if
      end do
      mu(ieh)=mu2
!      write(6,*) 'finish',mu,dense

      return
      end subroutine shiftmuSHF
! *********************************************************************







! *********************************************************************
      subroutine shiftmu(num,ieh,mu,dmu0,numeh,&
           beta,Dos,xw,divw,nwn1,nwe2,nwnn1,nwee2,nwmax)
      implicit none
      integer(4) :: i,j,sign,iw,nwmax,nwn1,nwe2,nwnn1,nwee2,ieh
      real(8),parameter :: pai=dacos(-1.0d0),hdig=1.0d-7
      real(8) :: &
        mu0,mu2,num,dmu0,beta,mumax,mumin,dense
      real(8) :: Dos(nwn1:nwe2,2),xw(nwnn1:nwee2),divw,&
      ffnc,w,mu(2),numeh(2)
      
! electron ------------------------------------------------------------
! check upper limit
      mu0=mu(ieh)
    ! total density
      dense=0.0d0
      do iw=nwn1,nwmax
        w=xw(iw)
        dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu0)
      end do
      dense=dense*divw

      if(dense.gt.num) then
        sign=-1
        mumax=mu0
      else
        sign=1
        mumin=mu0
      end if
! narrow down the chemical potential
      i=0
      do
        i=i+sign
        mu2=mu0+dble(i)*dmu0
    ! total density        
        dense=0.0d0
        do iw=nwn1,nwmax
          w=xw(iw)
          dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu2)
        end do
        dense=dense*divw
        
        if(dense.lt.num) then
          mumin=mu2
          if(sign.eq.-1) exit
        else
          mumax=mu2
          if(sign.eq.1) exit
        end if
        
      end do
      
      i=0
      do
        i=i+1
        mu2=(mumax+mumin)/2.0d0
    ! total density        
        dense=0.0d0
        do iw=nwn1,nwmax
          w=xw(iw)
          dense=dense+Dos(iw,ieh)*ffnc(beta,w,mu2)
        end do
        dense=dense*divw
        
!        write(6,*) (mumax+mumin)/2.0d0,dense
        if(abs(dense-num)/num.lt.hdig) then
!          write(6,*) (mumax+mumin)/2.0d0
          exit
        else if(dense.lt.num) then
          mumin=mu2
        else
          mumax=mu2
        end if
        if(i.gt.100) then
! chemical potential located in the energy gap of finite size mesh pt.
          write(6,*) 'step size of wavenumber (1/nk) is too large'
          stop
        end if
      end do
      mu(ieh)=mu2
      numeh(ieh)=dense

      return
      end subroutine shiftmu
! **********************************************************************



! **********************************************************************
      function screening(beta,engk,nk,xyk,fermi0)
      implicit none
      integer(4) :: nk,ik,ieh
      real(8) :: beta,screening,engk(0:nk,2),xyk(0:nk),fermi0(0:nk,2)
      
      screening=0.0d0
      do ieh=1,2 
        do ik=0,nk ! k
          screening=screening+&
          xyk(ik)*fermi0(ik,ieh)*(1.0d0-fermi0(ik,ieh))
        end do
      end do
      close(24)
      screening=2.0d0*beta*screening

      return
      end function screening
! **********************************************************************



! **********************************************************************
      subroutine Intexp(mlim,mmax,nt,nk,Wscexp,yt,ksc,cffInt,nt0,&
      n0ptk,x0ptk,y0ptk,divk,xk,xt)
      implicit none
      integer(4) :: nk,ik1,ik2,it,mmax,mlim,jm,ieh1,ieh2,nt0
      integer(4) :: nt(0:mmax),n0ptk
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: &
        k1,k2,t0,q,const,Opl2,ksc,Wscexp(0:nk,0:nk,0:mmax),&
        yt(nt0*(mmax+1),0:mmax),S1,Stot,cffInt,x0ptk(n0ptk),y0ptk(n0ptk),&
        St(500,0:500),xint,divk,xk(0:nk),xt(nt0*(mmax+1),0:mmax)
    
    ! interaction (no 0poit correction)
      Wscexp=0.0d0
      do jm=0,mlim ! m
        do ik2=0,nk ! k2
          k2=xk(ik2)
          do ik1=0,nk ! k1
            k1=xk(ik1) 
            S1=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm) 
              q=dsqrt(k1**2+k2**2-2.0d0*k1*k2*dcos(t0)) ! |q|
              S1=S1+&
              dcos(dble(jm*xt(it,jm)))*yt(it,jm)/(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Wscexp(ik1,ik2,jm)=S1/pai
          end do
        end do
      end do
 
    ! interaction with 0point correction on 0<ik1=ik2<nk
      do jm=0,mlim
        do ik1=1,nk-1 ! k1
          k1=xk(ik1)
          Stot=0.0d0
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm)
              q=dsqrt(k1**2+(k1-x0ptk(ik2))**2-&
              2.0d0*k1*(k1-x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+dcos(dble(jm*xt(it,jm)))*yt(it,jm)&
                /(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/pai
          end do
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm)
              q=dsqrt(k1**2+(k1+x0ptk(ik2))**2-&
              2.0d0*k1*(k1+x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+dcos(dble(jm*xt(it,jm)))*yt(it,jm)&
                /(q+ksc/(1.0d0+q**2*cffInt))
            end do
            Stot=Stot+y0ptk(ik2)*S1/pai
          end do
          
        ! remove adjacent contribution
          Wscexp(ik1,ik1,jm)=&
          (Stot-divk*(Wscexp(ik1,ik1-1,jm)+Wscexp(ik1,ik1+1,jm))/2.0d0)/divk
        end do
      end do

    ! unused interaction at ik1=ik2=0
      do jm=0,mlim
         Wscexp(0,0,jm)=0.0d0
      end do
      Wscexp=2.0d0*pai*Wscexp

      return
      end subroutine Intexp
! **********************************************************************




! **********************************************************************
      subroutine Int2nd(mlim,mmax,nt,nk,WscV0exp,&
        yt,ksc,cffInt,nt0,n0ptk,x0ptk,y0ptk,divk,xk,xt)
      implicit none
      integer(4) :: nk,ik1,ik2,it,mmax,mlim,jm,ieh1,ieh2,nt0
      integer(4) :: nt(0:mmax),n0ptk
      real(8),parameter :: pai=dacos(-1.0d0)
      real(8) :: &
        k1,k2,t0,q,const,Opl2,ksc,WscV0exp(0:nk,0:nk,0:mmax),&
        yt(nt0*(mmax+1),0:mmax),S1,S2,Stot,cffInt,x0ptk(n0ptk),y0ptk(n0ptk),&
        St(500,500),xint,divk,xk(0:nk),xt(nt0*(mmax+1),0:mmax)

    ! Wsc*V (no 0poit correction)    
      WscV0exp=0.0d0
      do jm=0,mlim ! m
        do ik2=0,nk ! k2
          k2=xk(ik2)
          do ik1=0,nk ! k1
            k1=xk(ik1)
            S1=0.0d0
            S2=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm) 
              q=dsqrt(k1**2+k2**2-2.0d0*k1*k2*dcos(t0)) ! |q|
              S1=S1+dcos(dble(jm*xt(it,jm)))*yt(it,jm)&
              /(q+ksc/(1.0d0+q**2*cffInt)) ! Wsc
              S2=S2+dcos(dble(jm*xt(it,jm)))*yt(it,jm)/q ! V0
            end do
            WscV0exp(ik1,ik2,jm)=S1*S2/pai**2
          end do
        end do
      end do
      
    ! ik1=ik2
      do jm=0,mlim ! m
        do ik1=1,nk-1 ! k
          k1=xk(ik1)
          Stot=0.0d0
          do ik2=1,n0ptk ! dk
            S1=0.0d0
            S2=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm)
              q=dsqrt(k1**2+(k1-x0ptk(ik2))**2-&
              2.0d0*k1*(k1-x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+dcos(dble(jm*xt(it,jm)))*yt(it,jm)&
                /(q+ksc/(1.0d0+q**2*cffInt)) ! Wsc
              S2=S2+dcos(dble(jm*xt(it,jm)))*yt(it,jm)/q ! V0
            end do
            Stot=Stot+y0ptk(ik2)*S1*S2/pai**2
          end do

          do ik2=1,n0ptk ! dk
            S1=0.0d0
            S2=0.0d0
            do it=1,nt(jm) ! theta_k1k2
              t0=xt(it,jm)
              q=dsqrt(k1**2+(k1+x0ptk(ik2))**2-&
              2.0d0*k1*(k1+x0ptk(ik2))*dcos(t0)) ! |q|
              S1=S1+dcos(dble(jm*xt(it,jm)))*yt(it,jm)&
                /(q+ksc/(1.0d0+q**2*cffInt)) ! Wsc
              S2=S2+dcos(dble(jm*xt(it,jm)))*yt(it,jm)/q ! V0
            end do
            Stot=Stot+y0ptk(ik2)*S1*S2/pai**2
          end do
          
        ! remove adjacent contribution
          WscV0exp(ik1,ik1,jm)=&
          (Stot-divk*(WscV0exp(ik1,ik1-1,jm)&
          +WscV0exp(ik1,ik1+1,jm))/2.0d0)/divk
        end do
      end do

    ! unused interaction at ik1=ik2=0
      do jm=0,mlim
        WscV0exp(0,0,jm)=0.0d0
      end do
      WscV0exp=(2.0d0*pai)**2*WscV0exp

      return
      end subroutine Int2nd
! **********************************************************************









! ***********************************************************************
      function bfnc(beta,ek)
      implicit none
      real(8) :: bfnc,beta,ek
         bfnc=1.0d0/(dexp(beta*ek)-1.0d0)
      end function bfnc
! ***********************************************************************
! ***********************************************************************
      function ffnc(beta,ek,mu)
      implicit none
      real(8) :: ffnc,beta,ek,mu
         ffnc=1.0d0/(dexp(beta*(ek-mu))+1.0d0)
      end function ffnc
! ***********************************************************************





! ***********************************************************************
      function CH(ksc,nCH,xCH,yCH,cffInt)
      implicit none
      integer(4) :: nCH,iq
      real(8) :: cffInt,ksc,CH,xCH(0:nCH),yCH(0:nCH),q
      real(8),parameter :: pai=dacos(-1.0d0)
      
      CH=0.0d0
      do iq=0,nCH
        q=xCH(iq)
        
        CH=CH+yCH(iq)*ksc/(ksc+q+cffInt*q**3)

        if(ksc/(ksc+q+cffInt*q**3)/CH.le.1.0d-10) exit
      end do
      CH=-pai*CH/(2.0d0*pai)
      
      return
      end function CH
! ***********************************************************************



! ***********************************************************************
      SUBROUTINE gauleg(x1,x2,x,w,n)
      real(8) :: pai=dacos(-1.0d0)
      INTEGER(4) :: n
      REAL(8) :: x1,x2,x(n),w(n)
      REAL(8),PARAMETER :: EPS=3.d-14
      INTEGER(4) :: i,j,m
      REAL(8) :: p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(pai*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END SUBROUTINE gauleg
! ***********************************************************************












! ***********************************************************************
      subroutine Dosnumber(Akw,engkf,engkQP,nk,nwg1,xw,xk,divw,&
      nwn1,nwe2,nwnn1,nwee2,self,Dos,xyk,numeh,numeh0,&
      mu,beta)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      complex(8),parameter :: im=(0.0d0,1.0d0)
      integer(4) :: &
      nk,nwg1(2),nwn1,nwe2,nwnn1,nwee2,ieh1,ik1,iw,iw1,iw2,&
      iph,iwph
      real(8) :: &
      Akw(nwn1:nwe2,0:nk,2),xk(0:nk),&
      engkf(0:nk,2),engkQP(0:nk,2),xw(nwnn1:nwee2),&
      divw,Dos(nwn1:nwe2,2),&
      xyk(0:nk),numeh(2),numeh0(2),mu(2),&
      beta,a,w,ww,wph,abc,def,ghi,jkl,ffnc,k1
      complex(8) :: self(nwn1:nwe2,0:nk,2)

    ! Akw
      Akw=0.0d0
      do ieh1=1,2 
         do ik1=0,nk
            do iw=nwn1,nwg1(ieh1)
               w=xw(iw)
               Akw(iw,ik1,ieh1)=&
                    dimag(1.0d0/(w-engkf(ik1,ieh1)-self(iw,ik1,ieh1)))
            end do
         end do
      end do
      Akw=-2.0d0*Akw


    ! Dos
      Dos=0.0d0
      do ieh1=1,2 
        do iw=nwn1,nwg1(ieh1)
          do ik1=0,nk
            Dos(iw,ieh1)=Dos(iw,ieh1)+xyk(ik1)*Akw(iw,ik1,ieh1)
          end do
        end do
      end do
      Dos=2.0d0*Dos/(2.0d0*pai)**2

      
    ! total density
      numeh=0.0d0
      do ieh1=1,2 
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          numeh(ieh1)=numeh(ieh1)+Dos(iw,ieh1)*ffnc(beta,w,mu(ieh1))
        end do
      end do
      numeh=divw*numeh
! ------------------------------------------------------------------
      
      
    ! quasi-particle
      do ieh1=1,2 
        do ik1=0,nk
          do iw=nwg1(ieh1)-1,nwn1,-1
            w=xw(iw)
            abc=(w-(engkf(ik1,ieh1)+dreal(self(iw,ik1,ieh1))))
            if(abc.lt.0.0d0) then
              ww=xw(iw+1)
              def=(ww-(engkf(ik1,ieh1)+dreal(self(iw+1,ik1,ieh1))))
              ghi=-abc/(def-abc)*(ww-w)
              exit
            end if
          end do
          engkQP(ik1,ieh1)=w+ghi
        end do
      end do
      
      open(44,file='T0.5n0.1_qp_pole_check.txt',status='unknown')
      do ieh1=1,2 
        do iw=nwn1,nwg1(ieh1)
          ww=xw(iw)
          def=(ww-(engkf(0,ieh1)+dreal(self(iw,0,ieh1))))
          write(44,*) ww,def
        end do
      end do
      close(44)
! -------------------------------------------------------------------


! quasi-particle density --------------------------------------------
      numeh0=0.0d0
      do ieh1=1,2 
        do ik1=0,nk
          numeh0(ieh1)=numeh0(ieh1)+&
          xyk(ik1)*ffnc(beta,engkQP(ik1,ieh1),mu(ieh1))
        end do
      end do
      numeh0=2.0d0*numeh0/(2.0d0*pai)
! -------------------------------------------------------------------






! Akw check -------------------------------------------------------
      open(20,file='T0.5n0.1_Akwshift_e.txt',status='unknown')
      open(21,file='T0.5n0.1_Akwshift_h.txt',status='unknown')

      do ik1=0,nk,nk/25
        k1=xk(ik1)
        write(20,*) 'k=',k1*a2d
        write(21,*) 'k=',k1*a2d
        ieh1=1
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          write(20,*) w*Esc/E2d,Akw(iw,ik1,ieh1)*E2d/Esc+k1*a2d
        end do
        ieh1=2
        do iw=nwn1,nwg1(ieh1)
          w=xw(iw)
          write(21,*) w*Esc/E2d,Akw(iw,ik1,ieh1)*E2d/Esc+k1*a2d
        end do
      end do

      close(20)
      close(21)
! ---------------------------------------------------------------


! Akw check2 ----------------------------------------------------------
      open(33,file='T0.5n0.1_Akwk0_e.txt',status='unknown')
      open(63,file='T0.5n0.1_Akwk0_h.txt',status='unknown')
      ieh1=1
      ik1=0
      k1=xk(0)
      write(33,*) 'k=',k1*a2d
      write(63,*) 'k=',k1*a2d
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(33,*) w*Esc/E2d,Akw(iw,ik1,ieh1)*E2d/Esc
      end do
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(63,*) w*Esc,Akw(iw,ik1,ieh1)*E2d/Esc
      end do
      close(33)
      close(63)
! --------------------------------------------------------------
      

! Dos check -------------------------------------------------------
      open(45,file='T0.5n0.1_Dos_e.txt',status='unknown')
      open(46,file='T0.5n0.1_Dos_h.txt',status='unknown')
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(45,*) w*Esc/E2d,Dos(iw,ieh1)/Esc*E2d
      end do
      ieh1=2
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        write(46,*) w*Esc/E2d,Dos(iw,ieh1)/Esc*E2d
      end do
      close(45)
      close(46)
! ------------------------------------------------------------------


! qp check ------------------------------------------------------
      open(43,file='T0.5n0.1_qp_dispersion_e.txt',status='unknown')
      open(44,file='T0.5n0.1_qp_dispersion_h.txt',status='unknown')
      ieh1=1
      do ik1=0,nk
        k1=xk(ik1)
        write(43,*) k1*a2d,engkQP(ik1,ieh1)*Esc/E2d
      end do
      ieh1=2
      do ik1=0,nk
        k1=xk(ik1)
        write(44,*) k1*a2d,engkQP(ik1,ieh1)*Esc/E2d
      end do
      close(43)
      close(44)
! ----------------------------------------------------------------

      
! fermi distribution check -------------------------------------------      
      open(176,file='T0.5n0.1_E_distribution.txt',status='unknown')
      open(177,file='T0.5n0.1_H_distribution.txt',status='unknown')
      do iw=nwn1,nwg1(1)
        w=xw(iw)
        write(176,*) w*Esc/E2d,ffnc(beta,w,mu(1))
      end do
      do iw=nwn1,nwg1(2)
        w=xw(iw)
        write(177,*) w*Esc/E2d,ffnc(beta,w,mu(2))
      end do
      close(176)
      close(177)
! ----------------------------------------------------------------

      return
      end subroutine Dosnumber
! **********************************************************************










! ************************************************************************
      subroutine output(nk,nwg1,xk,xw,nwn1,nwe2,nwnn1,nwee2,&
        engkQP,ksc,Akw,numeh,numeh0,mu,self)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: &
      ieh1,ieh2,ik1,nk,iw,nwg1(2),nwn1,nwe2,nwnn1,nwee2,&
      nave,i,j
      real(8) :: xk(0:nk),xw(nwnn1:nwee2),&
      engkQP(0:nk,2),Akw(nwn1:nwe2,0:nk,2),&
      numeh(2),numeh0(2),mu(2),k1,w,ksc
      complex(8) :: self(nwn1:nwe2,0:nk,2)
      
! selfenergy output for read-in ------------------------------------------
      open(89,file='T0.5n0.1_results.txt',status='unknown')
      write(89,*) mu(1),mu(2)
      
      do ieh1=1,2 
        do ik1=0,nk
          do iw=nwn1,nwg1(ieh1)
            write(89,*) self(iw,ik1,ieh1)
          end do
        end do
      end do
      close(89)
! ----------------------------------------------------------------------------


! output Akw (contour plot) -------------------------------------------------
      open(25,file='T0.5n0.1_Akwcontour_e.txt',status='unknown')
      open(26,file='T0.5n0.1_Akwcontour_h.txt',status='unknown')
      ieh1=1
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        if(abs(w).gt.10.0d0) goto 20
        do ik1=0,nk
          k1=xk(ik1)
          if(abs(k1).gt.10.0d0) exit
          write(25,'(3F20.15)') k1*a2d,w*Esc/E2d,&
          Akw(iw,ik1,ieh1)/Esc*E2d
          
        end do
        write(25,*) 
20    end do
      ieh1=2
      do iw=nwn1,nwg1(ieh1)
        w=xw(iw)
        if(abs(w).gt.10.0d0) goto 21
        do ik1=0,nk
          k1=xk(ik1)
          if(abs(k1).gt.10.0d0) exit
          write(26,'(3F20.15)') k1*a2d,w*Esc/E2d,&
          Akw(iw,ik1,ieh1)/Esc*E2d
          
        end do
        write(26,*) 
21    end do
      close(25)
      close(26)
! -------------------------------------------------------------------
      
      
! thermodynamic results -----------------------------------------------------      
      open(47,file='T0.5n0.1_thermodynamic_results.txt',status='unknown')
      write(47,*) 'E chemical potential (E2d)'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,mu(1)*Esc/E2d
      write(47,*) 'H chemical potential (E2d)'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,mu(2)*Esc/E2d
      write(47,*) 'total chemical potential (E2d)'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,(mu(1)+mu(2))*Esc/E2d
      write(47,*)
      
      write(47,*) 'E-H carrier density (1/a2d^-2)'
      write(47,*) numeh(1)*a2d**2,numeh(2)*a2d**2,&
      (numeh(1)+numeh(2))*a2d**2/2.0d0
      write(47,*) 'E-H rs parameter'
      write(47,*) &
        1.0d0/dsqrt(pai*numeh(1))/a2d,1.0d0/dsqrt(pai*numeh(2))/a2d,&
        (1.0d0/dsqrt(pai*numeh(1))/a2d+1.0d0/dsqrt(pai*numeh(1))/a2d)/2.0d0
      write(47,*)

      write(47,*) 'electron inonization ratio'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,&
      numeh0(1)/numeh(1)
      write(47,*) 'hole inonization ratio'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,&
      numeh0(2)/numeh(2)
      write(47,*) 'mean inonization ratio'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,&
      (numeh0(1)/numeh(1)+numeh0(2)/numeh(2))/2.0d0
      write(47,*)
      
      write(47,*) 'screening wave number (1/a2d)'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,ksc*a2d
      
      write(47,*) 'BGR (E2d)'
      write(47,*) (numeh(1)+numeh(2))*a2d**2/2.0d0,&
        (engkQP(0,1)+engkQP(0,2))*Esc/E2d
      
      close(47)
! --------------------------------------------------------------------------      
      
      return
      end subroutine output
! ****************************************************************************






! ***************************************************************************
      subroutine optical(nk,nwnn1,nwee2,nweh2,Wscexp,mmax,opt,opt0,bopt,&
        xw,xyk,g2,vertex,matrix,IPIV,beta,mu,num,divw)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: &
        nk,IPIV(nk+1),nwnn1,nwee2,nweh2,i,&
        mmax,iw0
      real(8) :: xyk(0:nk),Wscexp(0:nk,0:nk,0:mmax),mu0,bopt(nwnn1:nweh2),&
        xw(nwnn1:nwee2),mu(2),w,bfnc,beta,num,divw,def,ghi,ainter,&
        minw,maxw,peakw
      complex(8) :: &
        g2(nwnn1:nweh2,0:nk),vertex(nk+1),matrix(nk+1,nk+1)
      
      integer(4) :: iq,ik1,ik2,iw,INFO
      complex(8) :: zabc,zdef,opt(nwnn1:nweh2),opt0(nwnn1:nweh2)

      do iw=nwnn1,nweh2
        w=xw(iw)
        
        matrix=0.0d0
        do ik2=0,nk ! k2
          matrix(ik2+1,ik2+1)=1.0d0
          zabc=xyk(ik2)*g2(iw,ik2)/(2.0d0*pai)
          do ik1=0,nk ! k1
            matrix(ik1+1,ik2+1)=matrix(ik1+1,ik2+1)&
            -Wscexp(ik1,ik2,0)*zabc
          end do
        end do
        
        vertex=1.0d0
        call ZGESV(nk+1,1,matrix,nk+1,IPIV,vertex,nk+1,INFO)
        if(INFO.ne.0) then
          write(6,*) 'ZGESV error INFO=',INFO
          stop
        end if

        zabc=0.0d0
        zdef=0.0d0
        do ik1=0,nk
          zabc=zabc+xyk(ik1)*vertex(ik1+1)*g2(iw,ik1)
          zdef=zdef+xyk(ik1)*g2(iw,ik1)
        end do
        opt(iw)=zabc
        opt0(iw)=zdef
        bopt(iw)=bfnc(beta,w-mu(1)-mu(2))*dimag(zabc)
      end do
      opt=2.0d0*opt/(2.0d0*pai)
      opt0=2.0d0*opt0/(2.0d0*pai)
      bopt=2.0d0*bopt/(2.0d0*pai)

! output optical response ----------------------------------------------
      do iw=nwnn1,nweh2
        w=xw(iw)
        if(abs(w*Esc/E2d).lt.20.0d0) then
        open(222,file='T0.5n0.1_Imx.txt',status='unknown')
        open(125,file='T0.5n0.1_Rex.txt',status='unknown')
        open(123,file='T0.5n0.1_PL.txt',status='unknown')
        open(32,file='T0.5n0.1_Imx0.txt',status='unknown')
        open(132,file='T0.5n0.1_Rex0.txt',status='unknown')
        open(35,file='T0.5n0.1_PL0.txt',status='unknown')
        open(42,file='T0.5n0.1_Imx_per_Imx0.txt',status='unknown')
        open(142,file='T0.5n0.1_Rex_per_Rex0.txt',status='unknown')
        write(222,*) w*Esc/E2d,dimag(opt(iw))*a2d**2*E2d/Esc
        write(125,*) w*Esc/E2d,dreal(opt(iw))*a2d**2*E2d/Esc
        write(123,*) w*Esc/E2d,bopt(iw)*a2d**2*E2d/Esc
        write(32,*) w*Esc/E2d,dimag(opt0(iw))*a2d**2*E2d/Esc
        write(132,*) w*Esc/E2d,dreal(opt0(iw))*a2d**2*E2d/Esc
        write(35,*) &
          w*Esc/E2d,&
          bfnc(beta,w-mu(1)-mu(2))*dimag(opt0(iw))*a2d**2*E2d/Esc
        write(42,*) w*Esc/E2d,&
          dimag(opt(iw))/dimag(opt0(iw))*a2d**2*E2d/Esc
        write(142,*) w*Esc/E2d,&
          dreal(opt(iw))/dreal(opt0(iw))*a2d**2*E2d/Esc
        end if
      end do
      close(222)
      close(125)
      close(123)
      close(32)
      close(132)
      close(35)
      close(42)
      close(142)
! ---------------------------------------------------------------------

     
! output pair susceptibility ------------------------------------------
      open(56,file='T0.5n0.1_pair_susceptibility.txt',status='unknown') 
      mu0=(mu(1)+mu(2))
      do iw=nwnn1,nweh2
        w=xw(iw)
        if(w-mu0.gt.0) then
          iw0=iw-1
          exit
        end if
      end do
      def=ainter(mu0-xw(iw0),dreal(opt(iw0)),dreal(opt(iw0+1)),divw)
      ghi=ainter(mu0-xw(iw0),dreal(opt0(iw0)),dreal(opt0(iw0+1)),divw)
      write(56,*) "Rex at w=mu (1/(a2d^2*E2d))" 
      write(56,*) num*a2d**2,def*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      write(56,*) "Rex0 at w=mu (1/(a2d^2*E2d))"
      write(56,*) num*a2d**2,ghi*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num*a2d**2,def/ghi,1.0d0/beta*Esc/E2d
      write(56,*) "Rex-Rex0 at w=mu (1/(a2d^2*E2d))"
      write(56,*) num*a2d**2,(def-ghi)*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      close(56)
! ---------------------------------------------------------------------


! peak position -------------------------------------------------------
      open(30,file='T0.5n0.1_peak_position.txt',status='unknown')
      call peakposition(dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "absorption_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
      call peakposition(-dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "gain_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
      call peakposition(bopt,xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "PL_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
! ---------------------------------------------------------------------

      end subroutine optical
! **************************************************************************







! ***************************************************************************
      subroutine optical2(nk,nwnn1,nwee2,nweh2,Wscexp,mmax,opt,opt0,bopt,&
        xw,xyk,g2,vertex,matrix,IPIV,beta,mu,ksc,cffInt,yt,xt,nt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,mr,divw,num)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      integer(4) :: &
        nk,IPIV(nk+1),nwnn1,nwee2,nweh2,i,&
        mmax
      real(8) :: xyk(0:nk),Wscexp(0:nk,0:nk,0:mmax),&
        xw(nwnn1:nwee2),mu(2),w,bfnc,beta,bopt(nwnn1:nweh2)
      complex(8) :: &
        g2(nwnn1:nweh2,0:nk),vertex(nk+1),matrix(nk+1,nk+1)
      
      integer(4) :: iq,ik,ik1,ik2,iw,INFO
      complex(8) :: zabc,opt(nwnn1:nweh2),opt0(nwnn1:nweh2)

! added part to optical ----------------------------------------------
    ! me0,mh0,mr="3d GaAs e,h,reduced mass"
      real(8),parameter :: me0=0.0665d0,mh0=0.457d0,mr0=me0*mh0/(me0+mh0)
    ! E3d0="3d GaAs binding energy (eV)"
      real(8),parameter :: E3d0=4.2d-3
    ! Eg and E3d="band gap (1.5eV) and 3d exciton Rydberg",mr="reduced mass"
      real(8) :: Eg=1.5d0,E3d,mr
    ! kmax2="k cutoff"
      real(8),parameter :: kmax2=100.0d0
    ! nt="number of angle point"
      integer(4) :: nt(0:mmax),nt0
    ! xt and yt="point and weight in gauleg"
      real(8) :: xt(nt0*(mmax+1),0:mmax),yt(nt0*(mmax+1),0:mmax)
  ! integral parameter for 0 point correction
    ! n0ptk="number of fine k points", nk2="number of k point"
      integer(4) :: n0ptk,nk2
    ! x0ptk and y0ptk="point and weight"
      real(8) :: x0ptk(n0ptk),y0ptk(n0ptk)

      integer(4) :: iw0
      real(8) :: ksc,cffInt,divk,abc,def,ghi,mu0,ainter,divw,num,&
           minw,maxw,peakw
      complex(8) :: zdef
      real(8),allocatable :: Wscexp2(:,:,:),xk2(:),xyk2(:),dipole(:),&
           engkf2(:)
      integer(4),allocatable :: IPIV2(:)
      complex(8),allocatable :: g22(:),matrix2(:,:),vertex2(:)
! --------------------------------------------------------------------

    ! set energy gap
      Eg=Eg/(E3d0*mr/mr0)/Esc
    ! nk2="number of k point"
      nk2=int(kmax2/divk)

      allocate(Wscexp2(0:nk2,0:nk2,0:mmax),xk2(0:nk2),xyk2(0:nk2),&
           dipole(0:nk2),engkf2(0:nk2),g22(0:nk2),matrix2(nk2+1,nk2+1),&
           vertex2(nk2+1),IPIV2(nk2+1))
    ! xk="k points"
      do ik=0,nk2 ! k
         xk2(ik)=divk*dble(ik)
      end do
    ! xyk="weight of radial integral of k"
      do ik=0,nk2
         xyk2(ik)=divk**2*dble(ik)
      end do
    ! free q=0 energy
      do ik1=0,nk2
        engkf2(ik1)=0.5d0*xk2(ik1)**2
      end do
    ! dipole matrix element
      do ik1=0,nk2
        dipole(ik1)=Eg/(Eg+engkf2(ik1))
      end do

    ! effective interaction 
      call Intexp(mmax,mmax,nt,nk2,Wscexp2,yt,ksc,cffInt,nt0,&
        n0ptk,x0ptk,y0ptk,divk,xk2,xt)

      do iw=nwnn1,nweh2
        w=xw(iw)
        if(w*Esc/E2d.lt.20.0d0) then

        do ik=0,nk
          g22(ik)=g2(iw,ik)
        end do
        do ik=nk+1,nk2
          g22(ik)=-1.0d0/(w-engkf2(ik))
        end do

        matrix2=0.0d0
        do ik2=0,nk2
          matrix2(ik2+1,ik2+1)=1.0d0
          zabc=xyk2(ik2)*g22(ik2)/(2.0d0*pai)
          do ik1=0,nk2
            matrix2(ik1+1,ik2+1)=matrix2(ik1+1,ik2+1)&
            -Wscexp2(ik1,ik2,0)*zabc
          end do
        end do

        do ik=0,nk2
           vertex2(ik+1)=dipole(ik)
        end do
        call ZGESV(nk2+1,1,matrix2,nk2+1,IPIV2,vertex2,nk2+1,INFO)
        if(INFO.ne.0) then
          write(6,*) 'ZGESV error INFO=',INFO,iw,xw(iw)*Esc/E2d
        end if

        zabc=0.0d0
        zdef=0.0d0
        do ik1=0,nk2
          zabc=zabc+dipole(ik1)*xyk2(ik1)*vertex2(ik1+1)*g22(ik1)
          zdef=zdef+dipole(ik1)*xyk2(ik1)*g22(ik1)
        end do
        opt(iw)=zabc
        opt0(iw)=zdef
        bopt(iw)=bfnc(beta,w-mu(1)-mu(2))*dimag(zabc)
        end if
      end do
      opt=2.0d0*opt/(2.0d0*pai)
      opt0=2.0d0*opt0/(2.0d0*pai)
      bopt=2.0d0*bopt/(2.0d0*pai)

! output optical response ----------------------------------------------
      do iw=nwnn1,nweh2
        w=xw(iw)
        if(abs(w*Esc/E2d).lt.20.0d0) then
        open(222,file='T0.5n0.1_Imx_mod.txt',status='unknown')
        open(125,file='T0.5n0.1_Rex_mod.txt',status='unknown')
        open(123,file='T0.5n0.1_PL_mod.txt',status='unknown')
        open(32,file='T0.5n0.1_Imx0_mod.txt',status='unknown')
        open(132,file='T0.5n0.1_Rex0_mod.txt',status='unknown')
        open(35,file='T0.5n0.1_PL0_mod.txt',status='unknown')
        open(42,file='T0.5n0.1_Imx_per_Imx0_mod.txt',status='unknown')
        open(142,file='T0.5n0.1_Rex_per_Rex0_mod.txt',status='unknown')
        write(222,*) w*Esc/E2d,dimag(opt(iw))*a2d**2*E2d/Esc
        write(125,*) w*Esc/E2d,dreal(opt(iw))*a2d**2*E2d/Esc
        write(123,*) w*Esc/E2d,bopt(iw)*a2d**2*E2d/Esc
        write(32,*) w*Esc/E2d,dimag(opt0(iw))*a2d**2*E2d/Esc
        write(132,*) w*Esc/E2d,dreal(opt0(iw))*a2d**2*E2d/Esc
        write(35,*) w*Esc/E2d,&
          bfnc(beta,w-mu(1)-mu(2))*dimag(opt0(iw))*a2d**2*E2d/Esc
        write(42,*) w*Esc/E2d,&
          dimag(opt(iw))/dimag(opt0(iw))*a2d**2*E2d/Esc
        write(142,*) w*Esc/E2d,&
          dreal(opt(iw))/dreal(opt0(iw))*a2d**2*E2d/Esc
        end if
      end do
      close(222)
      close(125)
      close(123)
      close(32)
      close(132)
      close(35)
      close(42)
      close(142)
! ---------------------------------------------------------------------

     
! output pair susceptibility -----------------------------------------
      open(56,file='T0.5n0.1_pair_susceptibility_mod.txt',status='unknown') 
      mu0=(mu(1)+mu(2))
      do iw=nwnn1,nweh2
        w=xw(iw)
        if(w-mu0.gt.0) then
          iw0=iw-1
          exit
        end if
      end do
      def=ainter(mu0-xw(iw0),dreal(opt(iw0)),dreal(opt(iw0+1)),divw)
      ghi=ainter(mu0-xw(iw0),dreal(opt0(iw0)),dreal(opt0(iw0+1)),divw)
      write(56,*) "Rex at w=mu (1/(a2d^2*E2d))" 
      write(56,*) num*a2d**2,def*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      write(56,*) "Rex0 at w=mu (1/(a2d^2*E2d))"
      write(56,*) num*a2d**2,ghi*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      write(56,*) "Rex/Rex0 at w=mu" 
      write(56,*) num*a2d**2,def/ghi,1.0d0/beta*Esc/E2d
      write(56,*) "Rex-Rex0 at w=mu (1/(a2d^2*E2d))"
      write(56,*) num*a2d**2,(def-ghi)*a2d**2*E2d/Esc,1.0d0/beta*Esc/E2d
      close(56)
! ---------------------------------------------------------------------


! peak position -------------------------------------------------------
      open(30,file='T0.5n0.1_peak_position_mod.txt',status='unknown')
      call peakposition(dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "absorption_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
      call peakposition(-dimag(opt),xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "gain_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
      call peakposition(bopt,xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      write(30,*) "PL_peak_position"
      write(30,*) num*a2d**2,peakw*Esc/E2d,minw*Esc/E2d,maxw*Esc/E2d
! ---------------------------------------------------------------------


      end subroutine optical2
! **************************************************************************




! **************************************************************************
      subroutine peakposition(spec,xw,nwnn1,nweh2,nwee2,peakw,minw,maxw)
      implicit none
      real(8) :: spec(nwnn1:nweh2),xw(nwnn1:nwee2),maxspec,peakw,minw,maxw,w
      integer(4) :: nwnn1,nweh2,nwee2,imin,i
      maxspec=0.0d0
      do i=nwnn1,nweh2
        if(spec(i).gt.maxspec) then
          maxspec=spec(i)
          peakw=xw(i)
        end if
      end do
      
      do i=nwnn1,nweh2
        if(spec(i).gt.0.5d0*maxspec) then
          minw=(xw(i)-xw(i-1))/(spec(i)-spec(i-1))*&
          (0.5d0*maxspec-spec(i-1))+xw(i-1)
          imin=i
          exit
        end if
      end do
      
      do i=imin,nweh2
        if(spec(i).lt.0.5d0*maxspec) then
          maxw=(xw(i)-xw(i-1))/(spec(i)-spec(i-1))*&
          (0.5d0*maxspec-spec(i-1))+xw(i-1)
          exit
        end if
      end do

      
      end subroutine peakposition
! **************************************************************************



! ************************************************************************
      subroutine Impair(gpq,Imbgpq,nweh2,nk,nt0,nwn1,nwe2,nwnn1,nwee2,&
        nwg1,Akw,pblock,fermi,xk,xw,yt0,npbd,nt0bd,divk,iplusk,kk,ioutkk,&
        angleave,divw,nwg2)
      implicit none
      real(8),parameter :: Esc=2.0d0,E2d=4.0d0,a2d=0.5d0,pai=dacos(-1.0d0)
      complex(8),parameter :: im=(0.0d0,1.0d0)
      integer(4) :: &
        nk,iw,iw1,iw2,nwn1,nwe2,nwg1(2),nwg2(3),ieh,ieh1,ieh2,iplusk(0:nk),&
        nwnn1,nwee2,nweh2,iq,ip,it0,nt0,npbd(0:2*nk,3,2),&
        ioutkk(0:nt0,0:nk,0:2*nk,3,2),ik1,ik2,&
        nt0bd(0:nk,0:2*nk,3,2)
      real(8) :: &
        w,f,S1,Akw(nwn1:nwe2,0:nk,2),fermi(nwn1:nwe2,2),&
        yt0(0:nt0),xk(0:nk),pblock(nwn1:nwe2,nwn1:nwe2,3),xw(nwnn1:nwee2),&
        dk,kk(0:nt0,0:nk,0:2*nk,3,2),ainter,divk,divw,&
        angleave(0:nk,0:2*nk,3),p,Imbgpq(nwnn1:nwee2,0:nk,0:2*nk,3)
      complex(8) :: gpq(nwnn1:nwee2,0:nk,0:2*nk,3)
    ! Akw as a function of (p,q)
      real(8),allocatable :: Apqw(:,:,:)


! angle average -------------------------------------------------------------
      allocate(Apqw(0:nt0,nwn1:nwe2,2))
      Apqw=0.0d0
      gpq=0.0d0
      Imbgpq=0.0d0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
            
            ! k1,k2->p,q transformation
!$omp parallel do private(it0,ik1,ik2,dk,iw)
              do it0=nt0bd(ip,iq,ieh,1),nt0bd(ip,iq,ieh,2) ! theta_pq
                ik1=ioutkk(it0,ip,iq,ieh,1)
                dk=kk(it0,ip,iq,ieh,1)-xk(ik1)
                do iw=nwn1,nwg1(ieh1)
                  Apqw(it0,iw,1)=&
                    ainter(dk,Akw(iw,ik1,ieh1),Akw(iw,iplusk(ik1),ieh1),divk)
                end do

                ik2=ioutkk(it0,ip,iq,ieh,2)
                dk=kk(it0,ip,iq,ieh,2)-xk(ik2)
                do iw=nwn1,nwg1(ieh2)
                  Apqw(it0,iw,2)=&
                    ainter(dk,Akw(iw,ik2,ieh2),Akw(iw,iplusk(ik2),ieh2),divk)
                end do
              end do
!$omp end parallel do

              do iw2=nwn1,nwg1(ieh2)
!$omp parallel do private(iw1,iw,S1)
                do iw1=nwn1,nwg1(ieh1)
                  iw=iw1+iw2
                  
                  S1=0.0d0
                  do it0=nt0bd(ip,iq,ieh,1),nt0bd(ip,iq,ieh,2)
                    S1=S1+yt0(it0)*Apqw(it0,iw1,1)*Apqw(it0,iw2,2)
                  end do
                  
                  gpq(iw,ip,iq,ieh)=gpq(iw,ip,iq,ieh)+pblock(iw1,iw2,ieh)*S1
                  Imbgpq(iw,ip,iq,ieh)=Imbgpq(iw,ip,iq,ieh)+&
                       fermi(iw1,ieh1)*fermi(iw2,ieh2)*S1
                end do
!$omp end parallel do
              end do
!              do iw=nwnn1,nwg2(ieh) ! local average
!                gpq(iw,ip,iq,ieh)=gpq(iw,ip,iq,ieh)/angleave(ip,iq,ieh)
!              end do

            end do
          end do
        end do
      end do
      deallocate(Apqw)
      gpq=-im*divw*gpq/pai/(4.0d0*pai) ! (/pai) pai average  2012/11/2
      Imbgpq=-divw*Imbgpq/pai/(4.0d0*pai) ! (/pai) pai average  2012/11/2
! -------------------------------------------------------------------



! check pair Green's function -------------------------------------
      open(75,file='T0.5n0.1_Imgpq_ee.txt',status='unknown')
      open(76,file='T0.5n0.1_Imgpq_eh.txt',status='unknown')
      open(78,file='T0.5n0.1_Imgpq_hh.txt',status='unknown')
      do ip=0,nk,nk/25 ! p
        p=xk(ip)
        write(75,*) 'p=',p*a2d
        write(76,*) 'p=',p*a2d
        write(78,*) 'p=',p*a2d
        do iw=nwnn1,nwg2(1)
          w=xw(iw)
          write(75,*) w*Esc/E2d,&
            -dimag(gpq(iw,ip,0,1))*E2d/Esc+p*a2d*0.25d0
        end do
        do iw=nwnn1,nwg2(2)
          w=xw(iw)
          write(76,*) w*Esc/E2d,&
            -dimag(gpq(iw,ip,0,2))*E2d/Esc+p*a2d*0.25d0
        end do
        do iw=nwnn1,nwg2(3)
          w=xw(iw)
          write(78,*) w*Esc/E2d,&
            -dimag(gpq(iw,ip,0,3))*E2d/Esc+p*a2d*0.25d0
        end do
      end do
      close(75)
      close(76)
      close(78)
! --------------------------------------------------------

      return
      end subroutine Impair
! *********************************************************************





! **************************************************************************
      subroutine TmatrixMW(ImTpq,ImbTpq,gpq,Imbgpq,nwee2,nwnn1,&
           nk,xyk,npbd,npnum,cffT,cffTex,Wscexp,Wsc2nd,&
        mmax,nwg2,dT2,T2,dmatrix,matrix,ImMW,ImbMW,IPIV,xk,xw)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0),Esc=2.0d0,E2d=4.0d0,a2d=0.5d0
      real(8),parameter :: sgn(3)=(/1.0d0,-1.0d0,1.0d0/)
      integer(4) :: &
        nwnn1,nwee2,nk,iw,iq,ip1,ip3,ip2,ip,jm,mmax,&
        npnum(0:2*nk,3),npbd(0:2*nk,3,2),nwg2(3),&
        INFO,IPIV(nk+1),ieh,ieh1,ieh2,i1,i2,i3
      real(8) :: &
        ImbTpq(nwnn1:nwee2,0:nk,0:2*nk,3),Imbgpq(nwnn1:nwee2,0:nk,0:2*nk,3),&
        ImTpq(nwnn1:nwee2,0:nk,0:2*nk,3),xyk(0:nk),ImMW(nk+1),ImbMW(nk+1),&
        Wscexp(0:nk,0:nk,0:mmax),Wsc2nd(0:nk,0:nk,0:mmax),abc,def,&
        cffT(0:mmax),cffTex(0:mmax,3),xk(0:nk),xw(nwnn1:nwee2),p,w,&
        dmatrix(nk+1,nk+1),dT2(nk+1,nk+1)
      complex(8) :: &
        matrix(nk+1,nk+1),T2(nk+1,nk+1),&
        gpq(nwnn1:nwee2,0:nk,0:2*nk,3),zabc
      
      
      IPIV=0
      ImTpq=0.0d0
      T2=0.0d0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do jm=0,mmax

!$omp parallel do private(iw,iq,matrix,dmatrix,ip1,ip3,ip,ip2,zabc,&
             T2,dT2,IPIV,INFO,abc,i2,i3,i1)
      do iw=nwnn1,nwg2(ieh)
        do iq=0,2*nk ! q

! Tmatrix ---------------------------------------------------
        ! set right hand side
          do ip2=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p2
            i2=ip2-npbd(iq,ieh,1)+1
            do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p1
              i1=ip1-npbd(iq,ieh,1)+1
              T2(i1,i2)=sgn(ieh)*Wscexp(ip1,ip2,jm)
            end do
          end do
          
        ! set left hand side
          matrix=0.0d0
          do ip3=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p3
            i3=ip3-npbd(iq,ieh,1)+1
            zabc=-sgn(ieh)*xyk(ip3)*gpq(iw,ip3,iq,ieh)/(2.0d0*pai)
            do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p1
              i1=ip1-npbd(iq,ieh,1)+1
              matrix(i1,i3)=Wscexp(ip1,ip3,jm)*zabc
            end do
            matrix(i3,i3)=matrix(i3,i3)+1.0d0
          end do
          
          call ZGESV(npnum(iq,ieh),npnum(iq,ieh),&
          matrix,nk+1,IPIV,T2,nk+1,INFO)
          
          do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
            i1=ip-npbd(iq,ieh,1)+1
            ImTpq(iw,ip,iq,ieh)=ImTpq(iw,ip,iq,ieh)+&
            2.0d0*cffT(jm)*dimag(T2(i1,i1))&
            -cffTex(jm,ieh)*dimag(T2(i1,i1))
          end do

! Tmatrix multipled by bose distribution function --------------------
        ! set right hand side
          dT2=0.0d0
          do ip2=npbd(iq,ieh,1),npbd(iq,ieh,2)
            i2=ip2-npbd(iq,ieh,1)+1
            do ip3=npbd(iq,ieh,1),npbd(iq,ieh,2)
              i3=ip3-npbd(iq,ieh,1)+1
              abc=sgn(ieh)*xyk(ip3)*dreal(T2(i3,i2))*&
                   Imbgpq(iw,ip3,iq,ieh)
              do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2)
                i1=ip1-npbd(iq,ieh,1)+1
                dT2(i1,i2)=dT2(i1,i2)+Wscexp(ip1,ip3,jm)*abc
              end do
            end do
          end do
          dT2=dT2/(2.0d0*pai)

        ! set left hand side
          dmatrix=0.0d0
          do ip3=npbd(iq,ieh,1),npbd(iq,ieh,2)
            i3=ip3-npbd(iq,ieh,1)+1
            abc=-sgn(ieh)*xyk(ip3)*&
              dreal(gpq(iw,ip3,iq,ieh))/(2.0d0*pai)
            do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2)
              i1=ip1-npbd(iq,ieh,1)+1
              dmatrix(i1,i3)=Wscexp(ip1,ip3,jm)*abc
            end do
            dmatrix(i3,i3)=dmatrix(i3,i3)+1.0d0
          end do

          call DGESV(npnum(iq,ieh),npnum(iq,ieh),dmatrix,nk+1,&
          IPIV,dT2,nk+1,INFO)
          
          do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2)
            i1=ip1-npbd(iq,ieh,1)+1
            ImbTpq(iw,ip1,iq,ieh)=ImbTpq(iw,ip1,iq,ieh)+&
              2.0d0*cffT(jm)*dT2(i1,i1)-cffTex(jm,ieh)*dT2(i1,i1)
          end do

        end do
      end do
!$omp end parallel do
          end do
        end do
      end do
      
    ! for Montroll-Word
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do jm=0,mmax
!$omp parallel do private(iw,iq,ImMW,ImbMW,ip1,ip3,ip,abc,def)
      do iw=nwnn1,nwg2(ieh)
        do iq=0,2*nk ! q
        
          ImMW=0.0d0
          ImbMW=0.0d0
          do ip3=npbd(iq,ieh,1),npbd(iq,ieh,2)
            abc=xyk(ip3)*dimag(gpq(iw,ip3,iq,ieh))/(2.0d0*pai)
            def=xyk(ip3)*Imbgpq(iw,ip3,iq,ieh)/(2.0d0*pai)
            do ip1=npbd(iq,ieh,1),npbd(iq,ieh,2)
              ImMW(ip1+1)=ImMW(ip1+1)+Wsc2nd(ip1,ip3,jm)*abc
              ImbMW(ip1+1)=ImbMW(ip1+1)+Wsc2nd(ip1,ip3,jm)*def
            end do
          end do
          
          do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
            ImTpq(iw,ip,iq,ieh)=&
            ImTpq(iw,ip,iq,ieh)+2.0d0*cffT(jm)*ImMW(ip+1)
            ImbTpq(iw,ip,iq,ieh)=&
            ImbTpq(iw,ip,iq,ieh)+2.0d0*cffT(jm)*ImbMW(ip+1)
          end do
          
        end do
      end do
!$omp end parallel do
          end do
        end do
      end do
      
! check imaginary part of Tmatrix -------------------------------------
      open(65,file='T0.5n0.1_ImTpq_ee.txt',status='unknown')
      open(66,file='T0.5n0.1_ImTpq_eh.txt',status='unknown')
      open(68,file='T0.5n0.1_ImTpq_hh.txt',status='unknown')
      do ip=0,nk,nk/25 ! p
        p=xk(ip)
        write(65,*) 'p=',p*a2d
        write(66,*) 'p=',p*a2d
        write(68,*) 'p=',p*a2d
        do iw=nwnn1,nwg2(1)
          w=xw(iw)
          write(65,*) w*Esc/E2d,&
            -ImTpq(iw,ip,0,1)*Esc/E2d*a2d**2+p*a2d*15.0d0
        end do
        do iw=nwnn1,nwg2(2)
          w=xw(iw)
          write(66,*) w*Esc/E2d,&
            -ImTpq(iw,ip,0,2)*Esc/E2d*a2d**2+p*a2d*15.0d0
        end do
        do iw=nwnn1,nwg2(3)
          w=xw(iw)
          write(68,*) w*Esc/E2d,&
            -ImTpq(iw,ip,0,3)*Esc/E2d*a2d**2+p*a2d*15.0d0
        end do
      end do
      close(65)
      close(66)
      close(68)
! ---------------------------------------------------------------------

      return
      end subroutine TmatrixMW
! **************************************************************************






! **************************************************************************
      subroutine &
      trnsTpq(ImTpq,ImbTpq,ImTkk,ImbTkk,nwnn1,nwee2,nt0,nk,nwg2,&
        ioutpq,pq,xp,xq,yt0,iplusk,iplusq,divk)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0)
      integer(4) :: &
        ieh,ieh1,ieh2,ik1,ik2,nk,iw,nwnn1,nwg2(3),nwee2,it0,nt0,ip,iq,&
        ioutpq(0:nt0,0:nk,0:nk,3,2),iplusk(0:nk),iplusq(0:2*nk)
      real(8) :: &
        ImTkk(nwnn1:nwee2,0:nk,0:nk,2,2),&
        ImbTkk(nwnn1:nwee2,0:nk,0:nk,2,2),pq(0:nt0,0:nk,0:nk,3,2),&
        dp,dq,yt0(0:nt0),divk,&
        ainter2,xp(0:nk),xq(0:2.0d0*nk),&
        ImTpq(nwnn1:nwee2,0:nk,0:2*nk,3),&
        ImbTpq(nwnn1:nwee2,0:nk,0:2*nk,3),S1,S2
        

      ImTkk=0.0d0
      ImbTkk=0.0d0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
           ieh=ieh+1

          do ik2=0,nk ! k2
            do ik1=0,nk ! k1
!$omp parallel do private(it0,ip,iq,dp,dq,iw,S1,S2)
              do iw=nwnn1,nwg2(ieh)
                S1=0.0d0
                S2=0.0d0
                do it0=0,nt0 ! theta_k1k2
                  ip=ioutpq(it0,ik1,ik2,ieh,1) ! p
                  iq=ioutpq(it0,ik1,ik2,ieh,2) ! q
                  dp=pq(it0,ik1,ik2,ieh,1)-xp(ip)
                  dq=pq(it0,ik1,ik2,ieh,2)-xq(iq)
                ! normal
                  S1=&
                  S1+&
                  yt0(it0)*ainter2(dp,dq,&
                  ImTpq(iw,ip,iq,ieh),&
                  ImTpq(iw,ip,iplusq(iq),ieh),&
                  ImTpq(iw,iplusk(ip),iq,ieh),&
                  ImTpq(iw,iplusk(ip),iplusq(iq),ieh),divk)
                ! multipled by bose distribution function
                  S2=&
                  S2+&
                  yt0(it0)*ainter2(dp,dq,&
                  ImbTpq(iw,ip,iq,ieh),&
                  ImbTpq(iw,ip,iplusq(iq),ieh),&
                  ImbTpq(iw,iplusk(ip),iq,ieh),&
                  ImbTpq(iw,iplusk(ip),iplusq(iq),ieh),divk)
                end do
                ImTkk(iw,ik1,ik2,ieh1,ieh2)=S1
                ImbTkk(iw,ik1,ik2,ieh1,ieh2)=S2
              end do
!$omp end parallel do
            end do
          end do
        end do
      end do
      ImTkk=ImTkk/pai
      ImbTkk=ImbTkk/pai
      
      return
      end subroutine trnsTpq
! **************************************************************************


! *******************************************************************
      subroutine ImselfTmatrixMW(selfT,ImTkk,ImbTkk,nwn1,nwe2,xyk,&
        nwnn1,nwee2,nk,nwg1,fermi,bose,Akw,divw)
      implicit none
      real(8),parameter :: pai=dacos(-1.0d0)
      complex(8),parameter :: im=(0.0d0,1.0d0)
      integer(4) :: &
        ik1,ik2,nk,nwg1(2),nwn1,nwe2,nwnn1,nwee2,&
        iw1,iw2,iw,ieh,ieh1,ieh2
      real(8) :: &
        xyk(0:nk),ImTkk(nwnn1:nwee2,0:nk,0:nk,2,2),&
        ImbTkk(nwnn1:nwee2,0:nk,0:nk,2,2),S1,S2,divw,&
        Akw(nwn1:nwe2,0:nk,2),fermi(nwn1:nwe2,2),bose(nwnn1:nwee2,3)
      complex(8) :: selfT(nwn1:nwe2,0:nk,2)
      
      selfT=0.0d0
! Tmatrix self energy
      do ieh2=1,2 
        do ieh1=1,2 
          if(ieh1.eq.ieh2) then
             if(ieh1.eq.1) ieh=1
             if(ieh1.eq.2) ieh=3
          else
             ieh=2
          end if

          do ik1=0,nk
!$omp parallel do private(iw1,iw2,ik2,iw,S1,S2)
        do iw1=nwn1,nwg1(ieh1) ! w1
          do iw2=nwn1,nwg1(ieh2) ! w2
            iw=iw1+iw2
            
            S1=0.0d0
            S2=0.0d0
            do ik2=0,nk ! k2
              S1=S1+&
              xyk(ik2)*ImTkk(iw,ik1,ik2,ieh1,ieh2)*Akw(iw2,ik2,ieh2)
              S2=S2+&
              xyk(ik2)*ImbTkk(iw,ik1,ik2,ieh1,ieh2)*Akw(iw2,ik2,ieh2)
            end do
            
            selfT(iw1,ik1,ieh1)=selfT(iw1,ik1,ieh1)&
                 +S1*fermi(iw2,ieh2)+S2
          end do
        end do
!$omp end parallel do
          end do
        end do
      end do
    ! 1.0d0/(2.0d0*pai)**2 from w,k integration
      selfT=im*divw*selfT/(2.0d0*pai)**2
      
      return
      end subroutine ImselfTmatrixMW
! ************************************************************************






! **************************************************************************
      function ainter2(xx1,yy1,y11,y12,y21,y22,divw)
      implicit none
      real(8) :: xx1,yy1,y11,y12,y21,y22,ainter2,divw
      
        ainter2=(y22-y12-y21+y11)*(xx1)*(yy1)/divw**2&
        +(y21-y11)*(xx1)/divw+(y12-y11)*(yy1)/divw+y11
      end function ainter2
! **************************************************************************




! ************************************************************************
      function ainter(xx1,y1,y2,divk)
      implicit none
      real(8) :: xx1,divk,y1,y2,ainter
        ainter=(y2-y1)/divk*xx1+y1
      end function ainter
! ******************************************************************



! ***********************************************************************
      subroutine principal(nw1,nw2,nwnn1,nwee2,xw,fnc,logrm0,Sw)
      implicit none
      integer(4) :: iw,iw2,nw1,nw2,nwnn1,nwee2
      real(8) :: Sw(nwnn1:nwee2),fnc(nwnn1:nwee2),xw(nwnn1:nwee2),&
        logrm0(nwnn1:nwee2,nwnn1:nwee2),w,S1,abc
        
      Sw=0.0d0
      do iw2=nw1,nw2-1
        abc=(fnc(iw2+1)-fnc(iw2))/(xw(iw2+1)-xw(iw2))
!$omp parallel do private(iw,w)
        do iw=nw1,nw2
          w=xw(iw)
          Sw(iw)=Sw(iw)+&
          logrm0(iw2,iw)*((w-xw(iw2))*abc+fnc(iw2))
        end do
!$omp end parallel do
      end do


      end subroutine principal
! **********************************************************************









! ********************************************************************
      subroutine &
      bosedis(bose,mu,nwnn1,nwee2,xw,beta,divw)
      implicit none
      integer(4) :: nwnn1,nwee2,ieh,ieh1,ieh2,iw,iw2,iw0
      real(8) :: &
        xw(nwnn1:nwee2),mu(2),bfnc,bose(nwnn1:nwee2,3),beta,mu0,&
        abc,S1,S2,S3,divw,w,a1,a2,a3,a4
      
    ! bose distribution function (no 0pt correction)
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iw=nwnn1,nwee2 ! w
            bose(iw,ieh)=bfnc(beta,xw(iw)-mu(ieh1)-mu(ieh2))
          end do
        end do
      end do


    ! bose distribution function (no 0pt correction)
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          mu0=mu(ieh1)+mu(ieh2) ! sum of mu

        ! set sign inversion point
          iw0=nwnn1
          do iw=nwnn1,nwee2 ! w
            if(xw(iw)-mu0.gt.0.0d0) then
              iw0=iw
              exit
            end if
          end do

          if(iw0.ne.nwnn1) then
          ! set close to mu point
            if(xw(iw0)-mu0.gt.divw*0.5d0) iw0=iw0-1

            a1=dexp(-beta*(xw(iw0-1)-mu0))-1.0d0
            a4=dexp(-beta*(xw(iw0+1)-mu0))-1.0d0
            S2=dlog(abs(a4/a1))/beta
            
          ! remove adjacent contribution
            bose(iw0,ieh)=&
            (S2-divw*(bose(iw0-1,ieh)+bose(iw0+1,ieh))/2.0d0)/divw
          end if

        end do
      end do
      
      end subroutine bosedis
! *************************************************************************************



! *************************************************************************************
      subroutine outpq(ioutpq,pq,nt0,nk,xk,xt0,mass,M,xp,xq)
      implicit none
      integer(4) :: ieh,ieh1,ieh2,ik1,ik2,nt0,nk,it0,ip,iq
      integer(4) :: ioutpq(0:nt0,0:nk,0:nk,3,2)
      real(8) :: pq(0:nt0,0:nk,0:nk,3,2),xk(0:nk),xt0(0:nt0),mass(2),M(3),&
           xp(0:nk),xq(0:2*nk)
      real(8) :: k1,k2,p,q,t0,abc

    ! set p and q
      pq=0.0d0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do ik2=0,nk ! k2
            k2=xk(ik2)
            do ik1=0,nk ! k1
              k1=xk(ik1)
              do it0=0,nt0 ! theta_k1k2
                t0=xt0(it0)
              ! set p
                abc=(mass(ieh2)*k1)**2+(mass(ieh1)*k2)**2&
                -2.0d0*mass(ieh1)*mass(ieh2)*k1*k2*dcos(t0)
                if(abc.lt.0.0d0) abc=0.0d0
                pq(it0,ik1,ik2,ieh,1)=dsqrt(abc)/M(ieh)
              ! set q
                abc=k1**2+k2**2+2.0d0*k1*k2*dcos(t0)
                if(abc.lt.0.0d0) abc=0.0d0
                pq(it0,ik1,ik2,ieh,2)=dsqrt(abc)
                
              end do
            end do
          end do
        end do
      end do
      
    ! set integer ip and iq (xk(ip)<p<xk(ip+1) and xk(iq)<q<xk(iq))
      ioutpq=0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do ik2=0,nk ! k2
            k2=xk(ik2)
            do ik1=0,nk ! k1
              k1=xk(ik1)
              do it0=0,nt0 ! theta_k1k2
                t0=xt0(it0)
                
                p=pq(it0,ik1,ik2,ieh,1)
                q=pq(it0,ik1,ik2,ieh,2)
              ! set ip
                do ip=0,nk-1 ! p
                  if(p.ge.xp(ip).and.p.lt.xp(ip+1)) ioutpq(it0,ik1,ik2,ieh,1)=ip
                end do
                if(p.eq.xp(nk)) ioutpq(it0,ik1,ik2,ieh,1)=nk
              ! set iq
                do iq=0,2*nk-1 ! q
                  if(q.ge.xq(iq).and.q.lt.xq(iq+1)) ioutpq(it0,ik1,ik2,ieh,2)=iq
                end do
                if(q.eq.xq(2*nk)) ioutpq(it0,ik1,ik2,ieh,2)=2*nk
                
              end do
            end do
          end do
        end do
      end do
      end subroutine outpq
! ***************************************************************************



! ***************************************************************************
      subroutine outkk(kk,nt0,nk,ioutkk,xk,xq,xp,xt0,yt0,npbd,nt0bd,npnum,angleave,mass,M)
      implicit none
      integer(4) :: ieh,ieh1,ieh2,iq,ip,nt0,nk,it0,ik1,ik2,icount
      integer(4) :: ioutkk(0:nt0,0:nk,0:2*nk,3,2),npbd(0:2*nk,3,2),nt0bd(0:nk,0:2*nk,3,2),&
           npnum(0:2*nk,3)
      real(8) :: kk(0:nt0,0:nk,0:2*nk,3,2),xq(0:2*nk),xp(0:nk),mass(2),M(3),&
           xk(0:nk),angleave(0:nk,0:2*nk,3),xt0(0:nt0),yt0(0:nt0)
      real(8) :: k1,k2,p,q,t0,abc

    ! set k1 and k2
      kk=0.0d0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            q=xq(iq)
            do ip=0,nk ! p
              p=xp(ip)
              do it0=0,nt0 ! theta_pq
                t0=xt0(it0)
              ! set k1
                abc=(mass(ieh1)/M(ieh)*q)**2+p**2+&
                2.0d0*mass(ieh1)/M(ieh)*q*p*dcos(t0)
                if(abc.le.0.0d0) abc=0.0d0
                kk(it0,ip,iq,ieh,1)=dsqrt(abc)
              ! set k2
                abc=(mass(ieh2)/M(ieh)*q)**2+p**2-&
                2.0d0*mass(ieh2)/M(ieh)*q*p*dcos(t0)
                if(abc.le.0.0d0) abc=0.0d0
                kk(it0,ip,iq,ieh,2)=dsqrt(abc)
                
              end do
            end do
          end do
        end do
      end do
      
    ! set integer ik1 and ik2 (xk(ik1)<k1<xk(ik1+1) and xk(ik2)<k2<xk(ik2+1))
      ioutkk=10000
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            q=xq(iq)
            do ip=0,nk ! p
              p=xp(ip)
              do it0=0,nt0 ! theta_pq
                t0=xt0(it0)
                
                k1=kk(it0,ip,iq,ieh,1)
                k2=kk(it0,ip,iq,ieh,2)
              ! set ik1
                do ik1=0,nk-1 ! k1
                  if(k1.ge.xk(ik1).and.k1.lt.xk(ik1+1)) ioutkk(it0,ip,iq,ieh,1)=ik1
                end do
                if(k1.eq.xk(nk)) ioutkk(it0,ip,iq,ieh,1)=nk
              ! set ik2
                do ik2=0,nk-1 ! k2
                  if(k2.ge.xk(ik2).and.k2.lt.xk(ik2+1)) ioutkk(it0,ip,iq,ieh,2)=ik2
                end do
                if(k2.eq.xk(nk)) ioutkk(it0,ip,iq,ieh,2)=nk
                
              end do
            end do
          end do
        end do
      end do
      
    ! theta_pq boundary
      npbd=0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=0,nk ! p
              
            ! search p lower limit
              icount=0
              do it0=0,nt0 ! theta_pq
                if(ioutkk(it0,ip,iq,ieh,1).ne.10000.and.&
                   ioutkk(it0,ip,iq,ieh,2).ne.10000) then
                   icount=icount+1
                end if
              end do
              
              if(icount.ne.0) then
                npbd(iq,ieh,1)=ip ! set p lower limit
                exit
              end if
              if(ip.eq.nk) then
                npbd(iq,ieh,1)=10000
                write(6,*) 'no ip points',iq,ieh1,ieh2
              end if
              
            end do
          end do
        end do
      end do
      
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            icount=0
            do ip=nk,npbd(iq,ieh,1),-1 ! p

            ! search p upper limit
              icount=0
              do it0=0,nt0 ! theta_pq
                if(ioutkk(it0,ip,iq,ieh,1).ne.10000.and.&
                   ioutkk(it0,ip,iq,ieh,2).ne.10000) then
                   icount=icount+1
                end if
              end do
              
              if(icount.ne.0) then
                npbd(iq,ieh,2)=ip ! set p upper limit
                exit
              end if
              
            end do
          end do
        end do
      end do
      
    ! nt0 boundary
      nt0bd=0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
              icount=0
              do it0=0,nt0 ! theta_pq
              
              ! search theta_pq lower limit
                if(ioutkk(it0,ip,iq,ieh,1).ne.10000.and.&
                   ioutkk(it0,ip,iq,ieh,2).ne.10000.and.icount.eq.0) then
                  nt0bd(ip,iq,ieh,1)=it0
                  icount=1
                end if
              
                if(icount.ne.0) then
                  nt0bd(ip,iq,ieh,1)=it0 ! set theta_pq lower limit
                  exit
                end if
                if(it0.eq.nt0) then
                  write(6,*) 'no it0 points',it0,ip,iq,ieh1,ieh2
                  stop
                end if
              
              end do
            end do
          end do
        end do
      end do
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
              icount=0
              do it0=nt0,nt0bd(ip,iq,ieh,1),-1 ! theta_pq
              
              ! search theta_pq upper limit
                if(ioutkk(it0,ip,iq,ieh,1).ne.10000.and.&
                   ioutkk(it0,ip,iq,ieh,2).ne.10000.and.icount.eq.0) then
                  nt0bd(ip,iq,ieh,2)=it0
                  icount=1
                end if
              
                if(icount.ne.0) then
                  nt0bd(ip,iq,ieh,2)=it0 ! set theta_pq upper limit
                  exit
                end if
                
              end do
            end do
          end do
        end do
      end do

    ! set number of p points
      npnum=0
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            icount=0
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
              icount=icount+1
            end do
            npnum(iq,ieh)=icount ! number of p points
          end do
        end do
      end do
      
    ! set overlap angle integration
      ieh=0
      do ieh2=1,2 
        do ieh1=ieh2,2 
          ieh=ieh+1
          do iq=0,2*nk ! q
            do ip=npbd(iq,ieh,1),npbd(iq,ieh,2) ! p
              abc=0
              do it0=nt0bd(ip,iq,ieh,1),nt0bd(ip,iq,ieh,2) ! theta_pq
                abc=abc+yt0(it0)
              end do
              angleave(ip,iq,ieh)=abc ! overlap angle integration
            end do
          end do
        end do
      end do
      end subroutine outkk
! ***********************************************************************
