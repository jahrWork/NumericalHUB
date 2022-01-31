module Finite_differences
 use Lagrange_interpolation
 use Non_uniform_grids 

  implicit none 
 
  
  !interface Derivative
  !    module procedure FD_Derivative3D, FD_Derivative2D, FD_Derivative1D
  !end interface 
  
  !interface FD_Derivative
  !    module procedure FD_Derivative3D, FD_Derivative2D, FD_Derivative1D
  !end interface 
 
 
  type Grid
 
    character(len=30) :: name
    real, allocatable :: Derivatives( :, :, :)
    integer :: N
    real, allocatable :: nodes(:) 
    
  end type 
 
  integer, save :: Order
  integer, parameter :: Nmax = 20 
  type (Grid), save :: Grids(1:Nmax)
  integer, save :: ind = 0

 

  

    contains 
 

!***********************************************************************************
!*  It computes the coefficients of the high order derivatives once in a life time 
!*           Derivatives(:,:,:) 
!*
! Authors : Juan A Hernandez (juanantonio.hernandez@upm.es)  
!           Pablo Sierra Heras (pablo.sierra@hotmail.com)
!***********************************************************************************
subroutine FD_Grid_Initialization( grid_spacing , direction, nodes, q ) 

          character(len=*),  intent(in) :: grid_spacing, direction 
          integer, intent(in) ::  q
          real, intent(inout) :: nodes(:)
          
          integer d, df
                 
     Order = q 
     
     if (grid_spacing == "uniform") then 
         
                call  Uniform_grid( nodes ) 
                
     elseif  (grid_spacing == "nonuniform") then
                call  Non_uniform_grid( nodes, Order )
     endif 
         
     d = 0  
     d = findloc( Grids(:) % name, direction, dim=1   ) 
     
     if (d == 0) then
     
       ind = ind + 1                
       Grids(ind) % N = size(nodes) - 1
       Grids(ind) % name = direction
      
       allocate(Grids(ind)%nodes(0:Grids(ind) % N ))      
       allocate(Grids(ind)%Derivatives(-1:Order, 0:Order, 0:Grids(ind)%N)) 
     
       Grids(ind) % nodes = nodes
     
       call High_order_derivatives( Grids(ind) % nodes, Order,  & 
                                    Grids(ind) % Derivatives    )   
       write(*,*) " Grid name = ", Grids(ind) % name 
     
     elseif (d > 0) then
     
       Grids(d) % N = size(nodes) - 1
       Grids(d) % name = direction
       
       deallocate(Grids(d) % nodes, Grids(d) % Derivatives)
     
       allocate(Grids(d) % nodes(0:Grids(d) % N ))      
       allocate(Grids(d) % Derivatives(-1:Order, 0:Order, 0:Grids(d)%N)) 
     
       Grids(d) % nodes = nodes
     
       call High_order_derivatives( Grids(d) % nodes, Order,  & 
                                    Grids(d) % Derivatives    ) 
     endif
end subroutine 



!******************************************************************************
!* Lagrange coefficients and their derivatives at x_nodes 
!* Order:  maximum order of derivative and width of the stencil 
!******************************************************************************
subroutine High_order_derivatives( z_nodes, Order, Derivatives) 

  real, intent(in) ::  z_nodes(0:)
  integer, intent(in) :: Order  
  real, intent(out) :: Derivatives(-1:Order, 0:Order, 0:size(z_nodes)-1)  
  
   integer :: N, j, s 
   real :: xp 
   
   N = size(z_nodes) - 1; 
   
   do j=0, N
   
    if (mod(Order,2)==0) then 
                            s = max( 0, min(j-Order/2, N-Order) )   
    else 
                            s = max( 0, min(j-(Order-1)/2, N-Order) )
    endif 
  
    xp = z_nodes(j) 
    Derivatives(-1:Order, 0:Order, j) =  & 
    Lagrange_polynomials( x = z_nodes(s:s+Order), xp = xp ) 
   
   enddo 
end subroutine 


!****************************************************************************************
!* Derivative 3D 
!****************************************************************************************
subroutine FD_Derivative3D( direction, coordinate, derivative_order, W, Wxi ) 

   character(len=*), intent(in) :: direction(1:3)
   integer, intent(in) :: coordinate
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   W(0:, 0:, 0:)
   real, intent(out)::   Wxi(0:, 0:, 0:) 
   
    
    integer :: i, j, k, d1=0, d2=0, d3=0, Nx, Ny, Nz 
    integer, allocatable :: sx(:), sy(:), sz(:)
    integer :: m 
   
    d1 = findloc( Grids(:) % name, direction(1), dim=1   ) 
    d2 = findloc( Grids(:) % name, direction(2), dim=1   ) 
    d3 = findloc( Grids(:) % name, direction(3), dim=1   ) 
    m = derivative_order
    
    if (d1 > 0 .and. d2 > 0 .and. d3 > 0) then
    Nx = Grids(d1) % N 
    Ny = Grids(d2) % N
    Nz = Grids(d2) % N  
    
    sx = Stencilv( Order, Nx ) 
    sy = Stencilv( Order, Ny ) 
    sz = Stencilv( Order, Nz ) 
   
    do i=0, Nx
      do j=0, Ny  
         do k=0, Nz    

          if     (coordinate == 1) then 
                
           Wxi(i,j,k) = dot_product( Grids(d1) % Derivatives(m, 0:Order, i), &
                                      W(sx(i):sx(i)+Order, j, k) );

          elseif (coordinate == 2) then 
                 
           Wxi(i,j,k)  = dot_product( Grids(d2) % Derivatives(m, 0:Order, j), &
                                           W(i, sy(j):sy(j)+Order, k) );
                                  
          elseif (coordinate == 3) then 
                
           Wxi(i,j,k) = dot_product( Grids(d3) % Derivatives(m, 0:Order, k), &
                                           W(i, j, sz(k):sz(k)+Order) );                          
          else
          
           write(*,*) " Error Derivative3D"
           stop    
            
          endif   
       
         enddo 
      enddo 
    enddo 
    deallocate( sx, sy, sz ) 
   
   else
    
        write(*,*) " Error Derivative3D"
        stop 
   
   end if    
       
       


end subroutine



!****************************************************************************************
!*  Derivative 2D by means of Lagrange interpolation
!
!   I(x,y)             = sum_{ij} ( W_{ij} Lagrange_i(x) Lagrange_j(y) 
!   (dI/dx)(x_i, y_j)  = sum_{ij} ( W_{ij} dLagrange_i(x_i) ) 
!   (dI/dy)(x_i, y_j)  = sum_{ij} ( W_{ij} dLagrange_j(x_j) ) 
!****************************************************************************************
subroutine FD_Derivative2D( direction, coordinate, derivative_order, W, Wxi ) 

   character(len=*), intent(in) :: direction(1:2)
   integer, intent(in) :: coordinate, derivative_order
   real, intent(in) ::   W(0:, 0:)
   real, intent(out)::   Wxi(0:, 0:) 
   
    integer :: i, j, d1, d2, Nx, Ny
    integer, allocatable :: sx(:), sy(:)
    integer :: k 

    d1 = 0  ;   d1 = findloc( Grids(:) % name, direction(1), dim=1 )
    d2 = 0  ;   d2 = findloc( Grids(:) % name, direction(2), dim=1 )
  
    k = derivative_order
    
 if (d1 > 0 .and. d2 > 0) then
 Nx = Grids(d1) % N
 Ny = Grids(d2) % N
 allocate( sx(0:Nx), sy(0:Ny) ) 
 sx = Stencilv( Order, Nx ) 
 sy = Stencilv( Order, Ny ) 
   
   do i=0, Nx
      do j=0, Ny 

        if     (coordinate == 1) then 
            
         Wxi(i,j) = dot_product( Grids(d1) % Derivatives(k, 0:Order, i), & 
                                 W(sx(i):sx(i)+Order, j) );

        elseif (coordinate == 2) then 
            
         Wxi(i,j) = dot_product( Grids(d2) % Derivatives(k, 0:Order, j), &
                                 W(i, sy(j):sy(j)+Order) );
        else
                write(*,*) " Error Derivative"
                stop   
        endif   
       
      enddo 
   enddo 
   deallocate( sx, sy ) 
   
   else
        write(*,*) " Error Derivative2D"
        write(*,*) "Grids =", Grids(:)% name, "direction =", direction 
        write(*,*) "d1 =", d1, "d2 =", d2
        stop 
   end if  
end subroutine


!****************************************************************************************
!* Derivative 1D
!****************************************************************************************
subroutine FD_Derivative1D( direction, derivative_order, W, Wxi, j ) 

   character(len=*), intent(in) :: direction
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   W(0:)
   real, intent(out)::   Wxi(0:)
   integer, optional, intent(in) :: j 
   
    integer :: i, d, N, i1, i2  
    integer, allocatable :: sx(:)
    integer :: k 

    d = 0    
    d = findloc( Grids(:) % name, direction, dim=1   )
    k = derivative_order
    
    if (d > 0) then
        
        N = Grids(d) % N
        allocate ( sx(0:N) )
        sx = Stencilv( Order, N )   
        
        if (present(j)) then 
                      i1 = j; i2 = j 
        else 
                      i1 = 0; i2 = N 
        end if 
        do i= i1,  i2 
            
            Wxi(i)=dot_product( Grids(d) % Derivatives(k, 0:Order, i), & 
                                W(sx(i):sx(i)+Order) )
        enddo 
        
       deallocate( sx ) 
       
    else
        write(*,*) " Error Derivative1D"; stop 
    endif 

end subroutine


end module 
