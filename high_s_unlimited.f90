      implicit none
      integer n_layers,layer,ilayer
      integer i,n
      integer seed_size
      complex(kind=16),allocatable::tau(:),h(:),s(:)   
      complex(kind=16),allocatable::w(:)
      complex(kind=16) recursive_heaviside,heavi,tt
      real,allocatable:: sr(:),si(:) 
      real dt
      integer,allocatable:: seed(:)
      character*5 ilayer5
      call random_seed(size=seed_size)
      allocate(seed(seed_size))
      seed=20260730
      call random_seed(put=seed)
      n_layers=14*60
      n=n_layers
      dt=float(n)/(60.*365.*24.0)
      print *,"dt=",dt
      allocate(sr(n+1),si(n+1))
      allocate(tau(0:n+1),s(n+1),h(n+1))
      allocate(w(n+1))
           call random_number(sr)
           sr=sr*11.q0+0.1q0
           s=qcmplx(sr,0.0) 
           s(n+1)=qcmplx(5.0,0.0)
           do i=1,n+1
           tau(i)=qcmplx(float(i)*dt,0.q0)
           enddo
           tau=tau/tau(n+1)
            tau(0)=qcmplx(1.q0,0.q0)
      open(11,file="real_high_s.txt")
      do i=n+1,1,-1
      write(11,111) qreal(tau(i)),qreal(s(i))
      enddo
      close(11) 
      open(11,file="tau.txt")
      do i=n+1,1,-1
      write(11,111) float(i)*dt,qreal(tau(i))
      enddo
      close(11) 
111 format(2(1x,e14.7))
      w=qcmplx(0.q0,0.q0)
      do ilayer=n,1,-1
      tt=tau(ilayer)
      h(ilayer)=recursive_heaviside(tt,tau,s,ilayer,n+1)
      w(ilayer)=1.q0-h(ilayer)
      enddo
      open(11,file="real_high_s_u.asc") 
      do ilayer=n,1,-1
      write(11,11) qreal(tau(ilayer)),qreal(w(ilayer)) 
      enddo
      close(11)
11    format(2(1x,e14.7))
      stop
      end 
!
      recursive function recursive_heaviside(t,tau,s,level,n_layers) result(u)
      implicit none
      complex(kind=16),intent(in)::t,tau(0:n_layers),s(*)
      integer, intent(in) :: level,n_layers
      complex(kind=16) :: u,h_input,inner,arg,heavi
      if(level.eq.n_layers) then
      h_input=-t+tau(level)
      arg=-2.q0/s(level)*h_input
      if(abs(arg).gt.1.d+10) then
      u=0.0
      else
      u=1.q0/(1.q0+exp(arg))
      endif
      else
      inner=recursive_heaviside(t,tau,s,level+1,n_layers)
      h_input=-t+tau(level)*inner
      arg=-2.q0/s(level)*h_input
      u=1.q0/(1.q0+exp(arg))
      endif
      end  function recursive_heaviside      
!
