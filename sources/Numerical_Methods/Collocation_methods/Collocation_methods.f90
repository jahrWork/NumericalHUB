module Collocation_methods 
    use Dependencies_BC
    use Finite_differences
    use Fourier_interpolation
    use Chebyshev_interpolation
    implicit none 
    
    
     private 
     public ::                    &
      Grid_Initialization,        &     ! Coefficients of derivatives
      Derivative,                 &     ! k-th derivative of u(:)
      FREE_BOUNDARY_CONDITION,    &     ! to impose free BCs
      PERIODIC_BOUNDARY_CONDITION,&     ! to impose periodic BCs
      method                     ,&     ! Fourier, FiniteDifferences
      CPU_time_BC,                & 
      Derivative_3D,              & 
      Derivative_2D,              & 
      Spectral_transform,         &
      iSpectral_transform    
    
     interface Derivative
        module procedure Derivative3D, Derivative2D, Derivative1D
     end interface 
     
     abstract interface 
        subroutine nderivative( grid, order, U, dU )
           character(len=*), intent(in) :: grid
           integer, intent(in) :: order 
           real, intent(in) ::   U(0:)
           real, intent(out)::   dU(0:) 
        end subroutine 
     end interface 
     
     
     
     character(len=40), save :: method
     real, save :: CPU_time_BC = 0
     
     type domain1D
         character(len=40) :: method
         character(len=1) :: grid_name
         integer :: coordinate 
     end type 
     
     type (domain1D), save ::domain(3)
     
    contains 
!***********************************************************************************
!*  It initializes grid coefficients of different collocation methods: 
!*    inputs:
!*        -grid_spacing:
!*          1) Finite differences: uniform or nonuniform 
!*          2) Fourier collocation methods 
!*          3) Chebyshev collocation methods 
!*        -direction: "x", "y", "z" 
!*        -nodes: input are only extrema x0 and xn 
!*        -q (optional): order of the interpolant only for finite differences 
!*     
!*    outputs:   
!*        -nodes: grid points  and internal objects containing coefficients 
!*                of the high order derivatives calculated at initialization
!*
! Authors : Juan A Hernandez (juanantonio.hernandez@upm.es)  Jan 2021
!***********************************************************************************
subroutine Grid_Initialization( grid_spacing, direction, nodes, q ) 
    character(len=*),  intent(in) :: grid_spacing, direction 
    real, intent(inout) :: nodes(0:)
    integer, optional, intent(in) ::  q
          
       integer :: i
       
       if (direction == "x") then 
                                  i = 1 
       elseif (direction == "y") then 
                                  i = 2                                    
       elseif (direction == "z") then 
                                  i = 3
       end if 
       
       domain(i) % method = grid_spacing 
       domain(i) % grid_name = direction 
       domain(i) % coordinate = i
   !    domain(i) % direction = direction 
                                  
       if (grid_spacing == "Fourier" ) then 
           
           method = grid_spacing 
           call Fourier_Grid_Initialization( direction, nodes ) 
           write(*,'(A)') " Fourier expansions to calculate spatial derivatives"
           
       elseif (grid_spacing == "Chebyshev" ) then 
           
           method = grid_spacing 
           call Chebyshev_Grid_Initialization( direction, nodes ) 
           write(*,'(A)') " Chebyshev expansions to calculate spatial derivatives"
           
       else 
             method = "FiniteDifferences" 
             call FD_Grid_Initialization( grid_spacing , direction, nodes, q ) 
             write(*,'(A13, i4,A)') " High order=", q, & 
                     " Finite Differences to calculate spatial derivatives"
       end if 
          
end subroutine    
 

!****************************************************************************************
!* It calculates derivative in one dimensional grids. 
!* Grid_initialization must be called previously in order to set the method, grid and order 
!*
!*   Inputs: 
!*      direction        : " x", ...
!*      derivative_order : 1, 2, ..
!*      W(0:Nx)          : values of the function W(x) where derivatives are calculated
!*      j optional       : index of the mesh. if present(j) then, derivative is only evaluated at x_j 
!*
!*   Output: 
!*   
!*       Wxi  = (d^k W / dx^k)_i      i=0,.... Nx  
!*   
!* Author: Juan A Hernandez (juanantonio.hernandez@upm.es) Jan 2021
!****************************************************************************************
subroutine Derivative1D( direction, derivative_order, W, Wxi, j ) 

   character(len=*), intent(in) :: direction
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   W(0:)
   real, intent(out)::   Wxi(0:) 
   !real, intent(in) ::   W(:)
   !real, intent(out)::   Wxi(:) 
   integer, optional, intent(in) :: j 
   
      integer :: d 
      character(len=40) :: m
   
      d = findloc( domain(:) % grid_name, direction, dim = 1) 
      m = domain(d) % method
      
      if (m == "Fourier" ) then 
          call Fourier_Derivative1D( direction, derivative_order, W, Wxi )
          
      else if (m == "Chebyshev" ) then 
          call Chebyshev_Derivative1D( direction, derivative_order, W, Wxi )   
          
      else 
          call FD_Derivative1D( direction, derivative_order, W, Wxi, j  ) 
      end if 
      
            
      !if (method == "Fourier" ) then 
      !    call Fourier_Derivative1D( direction, derivative_order, W, Wxi )
      !    
      !else if (method == "Chebyshev" ) then 
      !    call Chebyshev_Derivative1D( direction, derivative_order, W, Wxi )   
      !    
      !else 
      !    call FD_Derivative1D( direction, derivative_order, W, Wxi, j  ) 
      !end if 
      
    
 end subroutine     

!****************************************************************************************
!*  Derivative 2D 
!****************************************************************************************
subroutine Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 

   character(len=*), intent(in) :: direction(1:2)
   integer, intent(in) :: coordinate, derivative_order
   real, intent(in) ::   W(0:, 0:)
   real, intent(out)::   Wxi(0:, 0:) 
   
      integer :: d 
      character(len=40) :: m
   
      d = findloc( domain(:) % grid_name, direction(coordinate), dim = 1) 
      m = domain(d) % method
      
      if (m == "Fourier" ) then 
         
          call Fourier_Derivative2D( direction, coordinate, derivative_order, W, Wxi )
          
          
      !else if (m == "Chebyshev" ) then 
      !    call Chebyshev_Derivative1D( direction, derivative_order, W, Wxi )   
      !    
      else 
          call FD_Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 
      end if 
   
   
   
      !if (method == "Fourier" ) then 
      !    call Fourier_Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 
      !else 
      !       call FD_Derivative2D( direction,  coordinate, derivative_order, W, Wxi ) 
      !end if   
      !
      
    
end subroutine
    
!****************************************************************************************
!* Derivative 3D 
!****************************************************************************************
subroutine Derivative3D( direction, coordinate, derivative_order, W, Wxi ) 

   character(len=*), intent(in) :: direction(1:3)
   integer, intent(in) :: coordinate
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   W(0:, 0:, 0:)
   real, intent(out)::   Wxi(0:, 0:, 0:) 
    
    
       call FD_Derivative3D( direction, coordinate, derivative_order, W, Wxi ) 
   

end subroutine









!****************************************************************************************
!*  Derivative 2D 
!****************************************************************************************
function Derivative_2D( grid, order, U) result (dU) 

   character(len=*), intent(in) :: grid
   integer, intent(in) ::  order
   real, intent(in) ::   U(0:, 0:)
   real :: dU( 0:size(U,dim=1)-1, 0:size(U, dim=2)-1 ) 
   
      integer ::  i, coordinate
      character(len=40) :: basis
         
      i = findloc( domain(:) % grid_name, grid, dim = 1) 
      basis = domain(i) % method
      coordinate = domain(i) % coordinate 
   !   write(*,*) "basis = ", basis, " coordinate =", coordinate 
      
      if (basis == "Fourier" ) then 
          
        dU = Directional_derivative2D( Fourier_Derivative1D, grid, coordinate, order, U)  
         
      else if (basis == "Chebyshev" ) then 
          
        dU = Directional_derivative2D( Chebyshev_Derivative1D, grid, coordinate, order, U)  
        
      else if (basis == "nonuniform" ) then 
            
        dU = Directional_derivative2D( FiniteD_Derivative1D, grid, coordinate, order, U)  
         
       else 
           write(*,*) "******* Error********* "
           stop 
           
      end if 
   
   
   
end function

function Directional_derivative2D(numerical_derivative, grid, coordinate, order, U ) result(dU)
   procedure(nderivative) numerical_derivative
   character(len=*), intent(in) :: grid
   integer, intent(in) :: coordinate, order
   real, intent(in) ::   U(0:, 0:)
   real :: dU( 0:size(U,dim=1)-1, 0:size(U,dim=2)-1) 

      integer :: i, j, N, M
     
      N = size( U, dim=1) 
      M = size( U, dim=2) 
   
      if (coordinate==1 ) then 
         
          do j=0, M-1 
            call numerical_derivative( grid, order, U(:,j), dU(:,j) )
          end do
          
      elseif (coordinate==2 ) then 
         
          do i=0, N-1
            call numerical_derivative( grid, order, U(i,:), dU(i,:) )
          end do
        
      end if 
   
   

end function    





!****************************************************************************************
!*  Derivative 3D 
!****************************************************************************************
function Derivative_3D( grid, order, U) result (dU) 

   character(len=*), intent(in) :: grid
   integer, intent(in) ::  order
   real, intent(in) ::   U(0:, 0:, 0:)
   real :: dU( 0:size(U,dim=1)-1, 0:size(U, dim=2)-1, 0:size(U, dim=3)-1 ) 
   
      integer ::  i, coordinate
      character(len=40) :: basis
         
      i = findloc( domain(:) % grid_name, grid, dim = 1) 
      basis = domain(i) % method
      coordinate = domain(i) % coordinate 
   !   write(*,*) "basis = ", basis, " coordinate =", coordinate 
      
      if (basis == "Fourier" ) then 
          
        dU = Directional_derivative( Fourier_Derivative1D, grid, coordinate, order, U)  
         
      else if (basis == "Chebyshev" ) then 
          
        dU = Directional_derivative( Chebyshev_Derivative1D, grid, coordinate, order, U)  
        
      else if (basis == "nonuniform" ) then 
            
        dU = Directional_derivative( FiniteD_Derivative1D, grid, coordinate, order, U)  
         
       else 
           write(*,*) "******* Error********* "
           stop 
           
      end if 
   
   
   
end function


function Directional_derivative(numerical_derivative, grid, coordinate, order, U ) result(dU)
   procedure(nderivative) numerical_derivative
   character(len=*), intent(in) :: grid
   integer, intent(in) :: coordinate, order
   real, intent(in) ::   U(0:, 0:, 0:)
   real :: dU( 0:size(U,dim=1)-1, 0:size(U,dim=2)-1, 0:size(U,dim=3)-1 ) 

      integer :: i, j, k,  N, M, L 
     
      N = size( U, dim=1) 
      M = size( U, dim=2)  
      L = size( U, dim=3)
   
      if (coordinate==1 ) then 
         
          do j=0, M-1; do k=0, L-1 
            call numerical_derivative( grid, order, U(:,j,k), dU(:,j,k) )
          end do; end do 
          
      elseif (coordinate==2 ) then 
         
          do i=0, N-1; do k=0, L-1 
            call numerical_derivative( grid, order, U(i,:,k), dU(i,:,k) )
          end do; end do  
        
      elseif (coordinate==3 ) then 
         
          do i=0, N-1; do j=0, M-1 
            call numerical_derivative( grid, order, U(i,j,:), dU(i,j,:) )
          end do; end do  
          
      end if 
   
   

end function      


!*************************************************************************
! Spectral transform 3D: Fourier x Fourier x Chebyshev  
!    
!        
!    1) Interpolant I(x,y,z) = sum_nml( c_nml exp(  I n x + I m y ) T_l( z )  )
!     
!    2) Conditions to obtain the interpolant     
!       u_ijk = sum_nml( c_nml exp(  2 PI n i/N I ) exp(  2 PI m j/M I ) T_l( z_k )    )
!
!       for i=0,...N-1, j=0,..M-1, k=0,...L  
!    
!    3) Coefficients determination     
!
!     c_nml = 1/(N M) 1/gamma_l 
!     sum_i{ sum_j[ sum_k( u_ijk T_l(z_k) ) exp( -2 PI n i/N I )] exp( -2 PI m j/M I )  }
!       
!  Author: Juan A. Hernandez, Jan., 2023   
!*************************************************************************
function Spectral_transform(N, M, L, U) result(C) 
     integer, intent(in) :: N, M, L 
     real, intent(in) :: U(0:N-1, 0:M-1, 0:L)
     complex :: C( -N/2:N/2-1, -M/2:M/2-1, 0:L ) 
        
     
     integer ::  i,j,k 
     complex :: Us(0:N-1, 0:M-1, 0:L)
     
     do i=0, N-1
         do j=0, M-1
           Us(i,j,:) = Chebyshev_transform( U(i,j,:) )
         end do 
     end do  
     
     do k=0, L 
            C(:,:,k) = FFT_2D(N, M, Us(:,:,k) )
     end do 
    
         
end function 
    
function iSpectral_transform(N, M, L, C) result(U) 
     integer, intent(in) :: N, M, L 
     complex, intent(in) :: C( -N/2:N/2-1, -M/2:M/2-1, 0:L )
     complex :: U(0:N-1, 0:M-1, 0:L) 
        
     
     integer ::  i,j,k 
     complex :: Cs( -N/2:N/2-1, -M/2:M/2-1, 0:L )
     
     
     do i=-N/2, N/2-1
         do j=-M/2, M/2-1
           Cs(i,j,:) = iChebyshev_transform( real(C(i,j,:)) )
         end do 
     end do  
     
     do k=0, L 
            U(:,:,k) = iFFT_2D(N, M, Cs(:,:,k) )
     end do 
    
         
end function 

    
    
end module 
    