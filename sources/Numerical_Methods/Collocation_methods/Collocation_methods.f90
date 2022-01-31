module Collocation_methods 
    use Dependencies_BC
    use Finite_differences
    use Fourier_interpolation
    implicit none 
    
     private 
     public ::                    &
      Grid_Initialization,        &     ! Coefficients of derivatives
      Derivative,                 &     ! k-th derivative of u(:)
      FREE_BOUNDARY_CONDITION,    &     ! to impose free BCs
      PERIODIC_BOUNDARY_CONDITION,&     ! to impose periodic BCs
      method                     ,&     ! Fourier, FiniteDifferences
      CPU_time_BC
    
     interface Derivative
        module procedure Derivative3D, Derivative2D, Derivative1D
     end interface 
     
     character(len=40), save :: method
     real, save :: CPU_time_BC = 0
     
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
          
       if (grid_spacing == "Fourier" ) then 
           
           method = grid_spacing 
           call Fourier_Grid_Initialization( direction, nodes ) 
           write(*,'(A)') " Fourier expansions to calculate spatial derivatives"
           
       elseif (grid_spacing == "Chebyshev" ) then 
           
      !     method = grid_spacing 
      !     call Chebyshev_Grid_Initialization( direction, nodes ) 
           write(*,'(A)') " Chebyshev expansions to calculate spatial derivatives (TODO)"
           stop
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
   real, intent(in) ::   W(:)
   real, intent(out)::   Wxi(:) 
   integer, optional, intent(in) :: j 
   
            
      if (method == "Fourier" ) then 
          call Fourier_Derivative1D( direction, derivative_order, W, Wxi ) 
      else 
          call FD_Derivative1D( direction, derivative_order, W, Wxi, j  ) 
      end if 
      
    
 end subroutine     

!****************************************************************************************
!*  Derivative 2D 
!****************************************************************************************
subroutine Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 

   character(len=*), intent(in) :: direction(1:2)
   integer, intent(in) :: coordinate, derivative_order
   real, intent(in) ::   W(0:, 0:)
   real, intent(out)::   Wxi(0:, 0:) 
   
      if (method == "Fourier" ) then 
          call Fourier_Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 
      else 
             call FD_Derivative2D( direction,  coordinate, derivative_order, W, Wxi ) 
      end if   
   
      
    
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



    
    
    
end module 
    