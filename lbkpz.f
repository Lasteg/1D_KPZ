c ================================
	program kpzD1Q3
c ================================
	use omp_lib
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c -------------------------------------------
        open(26,file='tip.out')
c --------------------------------------
        call input
        iseed = 95315571
        do ir=1,nreal
         iseed = iseed + 753 
         call init

         visco = 1./3.*(1./omega-0.5)
         Dramp = ramp*ramp/3.
         g     = lamb*lamb*Dramp/(visco**3)
         write(6,*) 'g coupling',g
         !pause

         call hydrovar
         call equili
    	 do istep = 1,nsteps
          call pbc
          call move
          call hydrovar
          call equili
          call colli
          call force(iseed)

          if (mod(istep,100000).eq.0) then
              write(6,*) 'calling diagno at realization ', ir
	      call diagno(istep)
          endif
	   if (mod(istep,ndiag).eq.0) then
            call profil(istep)
	   endif
        enddo

       end do

        stop
	end
c ==================================
	subroutine input
c ==================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c---------------------------------------------------
        iunit=7
        open(unit=iunit,file='lbkpz.inp') 

        read(iunit,*) omega
        read(iunit,*) rho0,u0
        read(iunit,*) lamb
        read(iunit,*) ramp
        read(iunit,*) nsteps,ndiag
        read(iunit,*) nreal
  	
	return
	end
c ===============================
	subroutine init
c ===============================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c--------------------------------------------------
        sigma = 10.
        w(0) = 2./3.
        w(1) = 1./6.
        w(2) = w(1)
	do i = 1, nx
          x=float(i-nx/2)/sigma
c         form = (x/sigma*sigma)*exp(-0.5*x*x)
c initialize with triangle
ctriangle          form=1.
ctriangle         if(i.gt.nx/2) form=-1.
c initialize zero h
          form = 0.

          rho(i) = rho0
	  u(i)   = u0
          f(0,i)=rho0*w(0)*form
          f(1,i)=rho0*w(1)*form
          f(2,i)=rho0*w(2)*form
        enddo
        
	return	
	end
c ==================================
	subroutine pbc
c periodic bc
c ==================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c --------------------------------------
        f(1,0)    = f(1,nx)
        f(2,nx+1) = f(2,1)

	return	
	end
c ==================================
	subroutine move
c ==================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c---------------------------------------------
        do i = nx,1,-1
	 f(1,i) = f(1,i-1)
        end do
        do i = 1,nx
	 f(2,i) = f(2,i+1)
        end do

	return	
	end
c =====================================
	subroutine hydrovar
c =====================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c----------------------------------------
	do i = 1, nx
          rho(i)=f(0,i)+f(1,i)+f(2,i) 
          u(i)  =(f(1,i)-f(2,i))/rho(i) 
	enddo

	return
	end
c ========================================
	subroutine equili
c ========================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c-------------------------------------------------
c rho is u in burgers
        cs2  = 1./3.
        do i = 1, nx
	 feq(0,i) = w(0)*rho(i)
	 feq(1,i) = w(1)*rho(i)*(1.0d0 + 0.5*lamb*rho(i)/cs2)
	 feq(2,i) = w(2)*rho(i)*(1.0d0 - 0.5*lamb*rho(i)/cs2)
        enddo
	return
	end
c =======================================================
	subroutine colli
c =======================================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
c----------------------------------------------------------
        do i = 1, nx
	  f(0,i) = f(0,i)-omega*(f(0,i)-feq(0,i))
	  f(1,i) = f(1,i)-omega*(f(1,i)-feq(1,i))
	  f(2,i) = f(2,i)-omega*(f(2,i)-feq(2,i))
	enddo

	return 
	end  
c =======================================================
	subroutine force(iseed)
c =======================================================
	implicit double precision(a-h,o-z)
	include'lbkpz.par'
        dimension rfrce(0:nx+1)
c----------------------------------------------------------
        frce = 0.
        do i=0,nx+1
         rfrce(i) = ramp*(2.*ranpang(iseed)-1)
        end do
c external force and langevin drag
        do i = 1, nx
         rf = 0.5*(rfrce(i+1)-rfrce(i-1))
         f(1,i) = f(1,i)+w(1)*(frce+rf)
         f(2,i) = f(2,i)-w(2)*(frce+rf)
        end do

	return 
	end
c ===================================
	subroutine diagno(istep)
c ===================================
	implicit double precision(a-h,o-z)
        include'lbkpz.par'
c----------------------------------------------------------
        rhotot=0.
        htot  =0.
        do i=1,nx
         rhotot = rhotot+rho(i)
         htot   = htot+h(i)
        end do
        write(6,*) 'total rho and h',istep, rhotot,htot
c roughness
        h(0)    = h(nx)
        h(nx+1) = h(1)
        hvar=0.
        do i=1,nx
         hbar=h(i+1)-2.*h(i)+h(i-1)
         hf = h(i)-hbar
         hvar = hvar + hf*hf
        end do
c tip motion and roughness
        write(26,*) istep,h(nx/2)-h(nx),sqrt(hvar/float(nx))         

	return
	end
c =======================
	subroutine profil(istep)
c =======================
	implicit double precision(a-h,o-z)
        include 'lbkpz.par'
c----------------------------------------------------------
c from burgers to kpz h=-grad u
        h(0)=0.
        do i=1,nx
        h(i+1) = h(i)-rho(i)
        end do    
        do i = 1,nx
         write(66,*) i,rho(i),h(i)
        end do
        write(66,'(bn)')
        write(66,'(bn)')

        return
        end
