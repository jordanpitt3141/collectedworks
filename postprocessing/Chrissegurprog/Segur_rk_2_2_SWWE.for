      program Segur_rk_2_2_SWWE

c============================================================
c Combine evolution of w = h+z and h for shallow flows
c Using Heun second-order Runga-Kutta
c Solving the sudden drop in bed experimental data from Hammack and Segur
c J. Fluid Mechanics, 84(2), 337-358, 
c Channel length is L = 50m, Q = 0.0m^3/s
c h0 = 10cm, u = 0m/s, b = 61cm and dh = 1.0, 3.0 and 5.0cm.
c============================================================

      implicit none

      integer md, nxmax, nxd, mn, i

      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)
      real*8 h(nxd), uh(nxd), z(nxd), z2(nxd), x(nxd),
     .       b, dh, L

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      real*8 toll, theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau
      common /rvariables/ toll, theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      call util_open_files
c
c Parameters
c
      nx = 1201
      tf = 50.d0
      nxx = nx + 2*md
      if (nx.gt.nxmax) write(*,*)'nx too large'

      toll  = 1.0d-10
      theta = 1.d0
      g     = 9.81d0
      h0    = 0.1d0
      u0    = 0.d0
      b     = 1.22d0/2.d0
      dh    = 0.01d0
      L     = 50.d0

      dx  = 2.d0*L/dble(nx-1)
      dt = dx/(dsqrt(g*h0))*.2d0
      tau = 0.0d0
      tc = 0.d0
c
c ---- Setup Terrain, Boundary Values and Initial Values
      do i = 1,nxx
        x(i)  = dble(i-md-1)*dx
        z(i)  = 0.d0
        z2(i) = 0.d0
      end do
c
c ---- Initial condition
      do i = 1, nxx
        h(i)  = h0
        uh(i) = h(i)*u0
        if (x(i).gt.x(nx+3)/2.d0-b.and.x(i).le.x(nx+3)/2.d0+b)then
          h(i) = h0 - dh
        end if
      end do
c
c ---- Boundary condition
      call boundary_apply(h,uh,x,tc)
c
c ---- Setup for Time looping
      isteps = 0
      
      call util_dump_solution(h,uh,z,x)
c
c ---- Evolve
      call evolver_rk_2(h,uh,z,z2,x)

      call util_dump_solution(h,uh,z,x)
c
c ---- Close down
      do i=1,nxx
c         write(16,*) x(i)
c         write(17,*) z(i)
      enddo
c      
c
      call util_close_files

      stop
      end

c===========================================================
c Shallow Water Flux
c===========================================================
      subroutine flux_huh(f1, f2, h0, uh0, em_x, l1, l2, g, toll)
      implicit none

      real*8 f1, f2, h0, uh0, l1, l2, em_x, g, toll, u, cvel, h, uh

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

c-----------------------------------

      uh = uh0
      h = h0

      if(h.le.toll)then
         u = 0.d0
         h = 0.0d0
         uh = 0.0d0
      else
         u = uh/h
      endif

      cvel = dsqrt( g*h )
      l1 = u - cvel
      l2 = u + cvel
      em_x = dabs(u) + cvel
      f1 = uh
      f2 = uh*u + 0.5d0*g*h**2

      return
      end

c===============================================================
c array limiter
c===============================================================
      subroutine array_limit(u, u_l, u_r)

      implicit none
      integer md, nxmax, nxd, mn
      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau
      common /rvariables/ toll,  theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      integer i

      real*8 a, b, t, xmin, xmic, u_x, d1, d2, u(nxd), u_l(nxd),
     .       u_r(nxd)

      xmin(a,b) = 0.5d0*(dsign(1.d0,a)+dsign(1.d0,b))*
     .            dmin1(dabs(a),dabs(b))
      xmic(t,a,b) = xmin(t*xmin(a,b), 0.5d0*(a+b) )

c------------------------------------

      do i=md,nx+md
          d1 = u(i)   - u(i-1)
          d2 = u(i+1) - u(i)
          u_x = xmic( theta, d1, d2 )
          u_l(i) = u(i) - 0.5d0*u_x
          u_r(i) = u(i) + 0.5d0*u_x
      enddo

      return
      end

c===========================================================
c Output Routine
c===========================================================
      subroutine util_dump_solution(h,uh,z,x)

      implicit none
      integer md,nxmax,nxd, mn
      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau
      common /rvariables/ toll,  theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      real*8 h(nxd), uh(nxd), z(nxd), x(nxd)

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      integer i

c      write(12,*)tc

      isteps = isteps + 1
      write(*,*)isteps, tc

      do i = 1, nx+2*md
c------ Dump h
         write(21,200) sngl(x(i)),sngl(h(i)),sngl(uh(i)/h(i))
  200 format(3f12.6)

      end do

      return
      end

c===============================================================
c The evolver using TVD Second-order Runga-Kutta method
c===============================================================
      subroutine evolver_rk_2(h,uh,z,z2,x)

      implicit none
      integer md,nxmax,nxd, mn
      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)
      real*8 h(nxd), uh(nxd), z(nxd), x(nxd)

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau, em_x
      common /rvariables/ toll,  theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      real*8 h_new(nxd),uh_new(nxd)
      real*8 h_star(nxd),uh_star(nxd)
      real*8 f_h(nxd), f_uh(nxd)
      real*8 h_l(nxd), h_r(nxd), uh_l(nxd), uh_r(nxd), z2(nxd)

      integer i, nt, i1, i2, i3, i4, i5

c----------------------------------------------------------
c
      nt =0

      do i=1,nxx
         f_h(i)  = 0.0d0
         f_uh(i) = 0.0d0
         h_l(i)  = h(i)
         h_r(i)  = h(i)
         uh_l(i) = uh(i)
         uh_r(i) = uh(i)
      enddo

      i1 = 50.61d0/dx + 5
      i2 = 55.61d0/dx + 5
      i3 = 60.61d0/dx + 5
      i4 = 65.61d0/dx + 5
      i5 = 70.61d0/dx + 5

      write(22,200) tc, h(i1), h(i2), h(i3), h(i4), h(i5)
        
      call boundary_apply(h,uh,x,tc)
      call array_limit  (uh, uh_l, uh_r)
      call array_limit  (h ,  h_l,  h_r)

      do i=1,nxx
         if ( h(i) .lt. 0.0d0 ) then
           write(*,*)'i=',i,' h(i) = ', h(i), tc
         endif
      enddo

c------------------------------------------------------
c
      do while (.true.)
c          write(*,*)nt
          nt = nt+1

c-----------------------------
c Interpolation and flux
c of first stage
c-----------------------------
          call boundary_apply(h,uh,x,tc)

          call array_limit  (uh, uh_l, uh_r)
          call array_limit  (h ,  h_l,  h_r)

          call numerical_flux_central(uh_l, uh_r, h_l, h_r, f_h, f_uh, 
     .                                em_x)

c-------------------------------------------------
c Compute the values of 'u' at the next time level
c-------------------------------------------------

          do i = md, md+nx

             h_star(i)  =  h(i)  - dt/dx*( f_h(i)  - f_h(i-1) )
             uh_star(i) =  uh(i) - dt/dx*( f_uh(i) - f_uh(i-1))

c--------------------------------------------------
c Bed Slope
c--------------------------------------------------
             uh_star(i) = uh_star(i)
     .              - dt/dx*g*(z2(i)-z2(i-1))*(h(i))


c--------------------------------------------------
c Friction
c--------------------------------------------------
             if(h(i).gt.0.d0)then
               uh_star(i) = uh_star(i) - dt*tau*uh(i)
             end if

             if (h_star(i).lt.-1.0d-5) then
                write(*,*)'h_star(',i,') = ',h_star(i)
             endif

             if (h_star(i).le.0.0d0) then
                h_star(i) = 0.d0
                uh_star(i) = 0.d0
             endif

          end do

c-----------------------------
c Interpolation and flux
c of second stage
c-----------------------------
          call boundary_apply(h_star,uh_star,x,tc+dt)

          call array_limit  (uh_star,uh_l, uh_r)
          call array_limit  (h_star , h_l,  h_r)

          call numerical_flux_central(uh_l, uh_r, h_l, h_r, f_h, f_uh, 
     .                                em_x)

c----------------------------
c Update second order
c----------------------------
          do i = md, md+nx

             h_new(i)  =  0.5d0*h_star(i) + 0.5d0*h(i)
     .                  - 0.5d0*dt/dx*( f_h(i)  - f_h(i-1) )
             uh_new(i) =  0.5d0*uh_star(i)+ 0.5d0*uh(i)
     .                  - 0.5d0*dt/dx*( f_uh(i) - f_uh(i-1))

c--------------------------------------------------
c Bed Slope
c--------------------------------------------------
             uh_new(i) = uh_new(i)
     .                  - 0.5d0*dt/dx*g*(z2(i)-z2(i-1))*(h_star(i))

c--------------------------------------------------
c Friction
c--------------------------------------------------
             if(h(i).gt.0.d0)then
               uh_new(i) = uh_new(i) - 0.5d0*dt*tau*uh_star(i)
             end if

             if (h_new(i).le.0.0d0) then
                h_new(i)  = 0.d0
                uh_new(i) = 0.d0
             endif

          end do

          do i = md,nx+md
             h(i) = h_new(i)
             uh(i) = uh_new(i)
          enddo

          write(22,200) tc, h(i1), h(i2), h(i3), h(i4), h(i5)
  200     format(11f12.6)

          tc = tc + dt

          write(*,*) tc

c          call util_dump_solution(h,uh,z,x)

          call boundary_apply(h,uh,x,tc)
          call array_limit  (uh, uh_l, uh_r)
          call array_limit  (h ,  h_l,  h_r)
      if(tc.ge.tf)go to 1001
c
c End timestep
c
      end do
c
 1001 write(*,*)isteps

      return

      end

c===========================================================
c Boundary Conditions
c===========================================================
      subroutine boundary_apply(h,uh,x,t)

      implicit none
      integer md, nxmax, nxd, mn
      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)
      real*8 h(nxd), uh(nxd), x(nxd), t

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau
      common /rvariables/ toll,  theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      integer i

c--------------------------------------------------


c==================
c Upstream Boundary
c==================
      do i = 1,md
        h(i)  = h0
        uh(i) = u0*h0
      end do

c==================
c Downstream Boundary
c==================
      do i = nx+md, nx+2*md
        h(i)  = h0
        uh(i) = u0*h0
      end do

      return
      end

c===========================================================
c Numerical Central Flux
c===========================================================
      subroutine numerical_flux_central(uh_l, uh_r, h_l, h_r,
     .                                  f_h, f_uh, em_x)

      implicit none
      integer md, nxmax, nxd, mn
      parameter (md=3, nxmax=15000, nxd=nxmax+2*md, mn=2)

      integer nx, nxx, version, isteps
      common /ivariables/ nx, nxx, version, isteps

      real*8 toll,  theta, g, dx, dt,
     .       tf, tc, hl, hr, uhl, uhr, h0, x0,
     .       u0, tau
      common /rvariables/ toll,  theta, g, dx, dt, 
     .                    tf, tc, hl, hr, uhl, uhr,
     .                    h0, x0, u0, tau

      real*8 f_h(nxd), f_uh(nxd), h_l(nxd), em_x,
     .       h_r(nxd), uh_l(nxd), uh_r(nxd), fh_l, fuh_l, fh_r, fuh_r,
     .       l1_l, l2_l, l1_r, l2_r, a_max, a_min, aa, axa

      integer i

c-----------------------------------

      em_x = 1.0d-15

      do i = 1,nxx-1
c
c ---- Calculate flux
c
         call flux_huh(fh_l,fuh_l,h_r(i),uh_r(i),
     .                 em_x,l1_l,l2_l,g,toll)
         call flux_huh(fh_r,fuh_r,h_l(i+1),uh_l(i+1),
     .                 em_x,l1_r,l2_r,g,toll)
c
c --- Calculate max wave speed
c
         a_max = dmax1(l2_l,l2_r,0.0d0)
         a_min = dmin1(l1_l,l1_r,0.0d0)
c
c --- Compute the numerical flux
c
         aa  = a_max - a_min
         axa = a_max*a_min
         if ( aa. gt. dx*dt/1000.0d0) then
            f_h(i)  =(a_max*fh_l  - a_min*fh_r  +
     .                axa*(h_l(i+1)-h_r(i)))/aa
            f_uh(i) =(a_max*fuh_l - a_min*fuh_r +
     .                axa*(uh_l(i+1)-uh_r(i)))/aa
         else
            f_h(i)  = 0.0d0
            f_uh(i) = 0.0d0
         endif

      enddo

      return
      end

c===========================================================
c Open Files
c===========================================================
      subroutine util_open_files

      open(21,file='Segur_rk_2_2_SWWE.out')
      open(22,file='Segur_rk_2_2_SWWE.r')

      return
      end

c===========================================================
c Close Files
c===========================================================
      subroutine util_close_files

      close(21)

      return
      end
