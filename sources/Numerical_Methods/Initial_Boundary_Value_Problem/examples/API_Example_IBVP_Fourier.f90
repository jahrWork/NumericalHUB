module API_Example_IBVP_Fourier 


use Linear_systems
use Initial_Boundary_Value_Problems
use Collocation_methods
use Temporal_Schemes
use plots 
implicit none 


contains 
 
subroutine Fourier_IBVP_examples 

       call Fourier_Burgers_equation_1D
       call Fourier_Advection_diffusion_2D 
           
end subroutine 
 
!********************************************************************************************
subroutine Fourier_Burgers_equation_1D

       integer, parameter :: Nx = 64 !255 ! 
       integer, parameter ::  Nt = 2000
       real ::  x(0:Nx-1)
       real :: Time(0:Nt), U(0:Nt,0:Nx-1)  
       
       real, parameter :: PI = 4 * atan(1.) 
       real ::  x0 = 0, xf = 2*PI*Nx/(Nx-1), t0 = 0, tf =  4
       integer :: i, j, k, q = 10  
       integer, parameter :: Nl = 16 
       character(len=10) :: legends(0:Nl) 
    
       
     write (*, '(A50)') 'Time solution of the 1D Burgers equation'
     write(*,*) "press enter "; read(*,*)
     
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
         
     x(0) = x0; x(Nx-1) = xf
     
  ! call Grid_Initialization( "nonuniform", "x", x, q )
    call Grid_Initialization( "Fourier", "x", x )
       
     U(0, :)  =  sin(x) 
    
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Burgers_equation1D,       & 
                       Boundary_conditions   =  Burgers_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt:Nt/Nl,:)), legends(1:Nl), "$x$", "$u(x,t)$") 
           
contains 

real function Burgers_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u, ux, uxx
            
         F =   -u * ux  + 0.03*  uxx
    !     F =   - ux  
           
end function 
!-------------------------------------------------------
real function Burgers_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u, ux 
  
         ! No need to impose BCs. They are periodic.
           BC = 0 
       
end function

end subroutine 


character(len=5) function float_to_str(x) result(c) 
      real, intent(in) :: x

    write(c,'(f5.1)')  x
    
end function






!*************************************************************************
subroutine Fourier_Advection_diffusion_2D

      integer, parameter :: Nx = 32, Ny = 32, Nt = 100
      real :: x(0:Nx-1), y(0:Ny-1), Time(0:Nt), U(0:Nt, 0:Nx-1, 0:Ny-1) 
       
       real, parameter :: PI = 4 * atan(1.) 
       real :: x0 = 0, xf = Nx/real(Nx-1), y0 = 0, yf = Ny/real(Ny-1)
       real :: t0 = 0, tf = 1.0
       integer :: i, j
       real :: levels(10)
                 
     write(*,*) "Time solution of the 2D advection diffusion equation "
     write(*,*) "with periodic bounday conditions( Fourier expansions)"
     write(*,*)  "press enter "; read(*,*)  
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0)= x0; x(Nx-1) = xf; y(0)=y0; y(Ny-1) = yf 
     do i=1, 10; levels(i) = 0.05*i; end do; 
     
     call Grid_Initialization( "Fourier", "x", x )
     call Grid_Initialization( "Fourier", "y", y )
  
    U(0, :, :) = Tensor_product( exp(-200*(x-0.5)**2), exp(-200*(y-0.5)**2) ) 
    call plot_contour(x, y, U(0,:,:), "x","y", levels, graph_type ="isolines") 
 !  stop 
 
!    Advection diffusion 2D      
     call Initial_Boundary_Value_Problem(                             & 
         Time_Domain = Time, x_nodes = x, y_nodes = y,                &
         Differential_operator =  Advection_equation2D,               & 
         Boundary_conditions   =  Advection_BC2D,  Solution = U ) 
     
     !write(*,*) " maxval U =", maxval(U(1,:,:)) 
     !stop 
     
     do i=0, Nt, Nt/10 + 1 
       call plot_contour(x, y, U(i,:,:), "x","y", levels, graph_type ="isolines") 
     end do
     
contains
!----------------------------------------------------
function Advection_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
           real,intent(in) :: x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy
           real :: F 
        
        real :: nu = 0.02

        F = - Ux - Uy  ! +  nu * ( Uxx + Uyy )
       

end function
!-------------------------------------------------------
real function Advection_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
          real, intent(in) :: x, y, t, U, Ux, Uy

          ! No need to immpodr BCs. They are periodic 
            BC =  0                       
     
end function

end subroutine  



end module 