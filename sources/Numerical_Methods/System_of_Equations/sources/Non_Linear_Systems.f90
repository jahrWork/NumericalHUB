module Non_Linear_Systems 

 use Jacobian_module 
 use Linear_systems
  
implicit none 

private 
public ::  & 
  Newton,  & ! It solves a vectorial system F(x) = 0
  Newtonc, & ! It solves a vectorial system G(x) = 0
             ! with M implicit equations < N unknowns
             ! e.g.  G1 = x1 - x2 (implicit) with x1 = 1 (explicit)
  max_min, & ! It finds a relative max or min of E(x)
  Gradient_descent ! It finds a relative min 


contains 

!***************************************************************************
!*  Newton solver 
!*                x0   : initial guess and output value 
!*                F(X) : vector function 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
subroutine Newton(F, x0) 
 
 procedure (FunctionRN_RN) :: F 
 real, intent(inout) :: x0(:) 
 
   real ::   Dx( size(x0) ), b(size(x0)), eps
   real :: J( size(x0), size(x0) )
   integer :: iteration, itmax = 1000 
   
   integer :: N 
   integer :: i

   N = size(x0) 

   Dx = 2 * x0 
   iteration = 0 
   eps = 1 
 
   do while ( eps > 1d-8 .and. iteration <= itmax )
    
      iteration = iteration + 1 
      J = Jacobian( F, x0 )       
     
      b = F(x0);
      Dx = Gauss( J, b) 
      x0 = x0 - Dx;  

      eps = norm2( DX )  
         
   end do 
   
   if (iteration == itmax) then 
      write(*,*) " morm2(J) =", maxval(J),  minval(J) 
      write(*,*) " Norm2(Dx) = ", eps, iteration  
   endif 
 
   
end subroutine

!***************************************************************************
!*  This version allows to solve: (see examples test_Newtonc)
!*           x(1) = x(2) 
!*           F(1) = 0
!*           F(2) = x(1) + x(2) - 1 
!*       
!*
!*  Newton solver 
!*                x0   : initial guess and output value 
!*                F(X) : vector function 
!*

!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
subroutine Newtonc(F, x0) 
 procedure (FunctionRN_RNE) :: F 
 real, intent(inout) :: x0(:) 

   real, allocatable  ::   x(:), Dx(:), b(:), J( :, : )
   real ::  eps, F0( size(x0) ), xr( size(x0) ) 
   integer :: iteration, itmax = 1000 
   logical :: equations( size(x0) )
  
  
   integer :: i, N
   
   
   call random_number(xr) 
   equations  = F(xr) /= 0  
   N = count( equations ) ! count valid equations 
       
   allocate ( x(N), Dx(N), b(N), J(N,N) ) 

   iteration = 0 
   eps = 1 
   x = pack( x0, equations )
  
   do while ( eps > 1d-8 .and. iteration <= itmax )
    
      iteration = iteration + 1 
      J = Jacobianc( G, x )
     
      b = -G(x) 
      call LU_factorization( J ) 
      
      Dx = Solve_LU( J,  b ) 
      x = x + Dx
    
      eps = norm2( Dx ) 
    
   end do 
   
   if (iteration==itmax) then 
      write(*,*) " morm2(J) =", maxval(J),  minval(J) 
      write(*,*) " Norm2(Dx) = ", eps, iteration  
   endif
   
   deallocate ( x, Dx, b, J ) 
  
   
contains 
!_____________________________________________
! Packs equations which are different from zero 
!_____________________________________________
function G( y ) 
      real, target :: y(:) 
      real :: G(size(y)) 
   
      integer :: i, m
    
      
       m = 1 
       do i=1, size(x0) 
           if (equations(i)) then 
                    x0(i) = y(m) 
                    m = m + 1 
           end if 
       end do 
    
       F0 = F( x0 ) 
       G = pack( F0, equations ) 
               
end function 

end subroutine


!***************************************************************************
!*  It calculates relative max and min of a function E 
!*    E: RN -> R 
!*                x    : initial guess and output value 
!*                E(X) : vector function 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
 subroutine max_min(E, x, Solver)  
     procedure (FunctionRN_R) :: E 
     real, intent(inout) :: x(:) 
     procedure (NonLinearSolver), optional:: Solver
     
     if (present(Solver)) then 
         
          call Solver(G, x) 
     else 
          call Newton( G, x ) 
          
     end if 
     
   
 contains
!_____________________________________________
! It calculates the set of equations: 
! derivative E with  respect to x_i: G_i = dE/dx_i  
! to find the local extrema by Newton or Solver 
!_____________________________________________ 
 function G(x) 
  real, intent(in) :: x(:) 
  real :: G(size(x)) 
  
   integer :: j, N 
   real :: dx(size(x)), h = 1d-6
  
  N = size(x) 
 
  do j=1, N       
            dx = 0 
            dx(j) = h 
            G(j) = ( E( x + dx ) - E( x - dx ) )/ ( 2 * h ) 
  end do 
  
  
end function    
 
 end subroutine 
 
!***************************************************************************
!*  It calculates relative min of a function E 
!* by the the gradient descent method
!*    E: RN -> R 
!*                x    : initial guess and output value 
!*                E(X) : vector function 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
 subroutine Gradient_descent(E, x)  
     procedure (FunctionRN_R) :: E 
     real, intent(inout) :: x(:) 
    
      integer :: iteration, itmax = 1000  
      real :: gamma, eps 
      real, allocatable :: grad(:), grad0(:), dx(:)   
     
     
   iteration = 0 
   
   eps = 1 
   dx = x 
   dx = 1d-2 
   grad0 = grad_E(x-dx) 
  
 
   do while ( eps > 1d-12 .and. iteration < itmax )
    
      iteration = iteration + 1 
      grad = grad_E(x)
      
      gamma = dot_product( dx, grad-grad0 ) / norm2( grad-grad0 )**2 
 
      dx = -gamma * grad
      x = x + dx 
      eps = norm2( dx )  
      grad0 = grad 

    !write(*,*) " x =", x
    !write(*,*) " gamma =", gamma 
    !write(*,*) " grad E =", grad 
    !read(*,*)
    
   end do 
   
   !if (iteration == itmax) then 
   !
   !   write(*,*) " morm2(grad_E) =", norm2(grad)
   !   write(*,*) " Norm2(Dx) = ", norm2(dx), "iteration =", iteration  
   !endif 
 
   
   
 contains
!_____________________________________________
! It calculates the set of equations: 
! derivative E with  respect to x_i: G_i = dE/dx_i  
! to find the local extrema by Newton or Solver 
!_____________________________________________ 
 function grad_E(x) 
  real, intent(in) :: x(:) 
  real :: grad_E(size(x)) 
  
   integer :: j, N 
   real :: dx(size(x)), h = 1d-2  
  
  N = size(x) 
 
  do j=1, N       
            dx = 0 
            dx(j) = h 
            grad_E(j) = ( E( x + dx) - E( x - dx ) )/ ( 2 * h ) 
  end do 
  
  
end function    
 
 end subroutine  


end module 
