module API_Example_Lagrange_Interpolation

    use Interpolation
    use Lagrange_interpolation
    use plots 
    
    implicit none
    real :: PI = 4 * atan(1d0)
    
    
    contains
    
  
subroutine Lagrange_Interpolation_examples

     call Interpolated_value_example 
     call Interpolant_example 
     call Integral_example 
     
     call Lagrange_polynomial_example 
     call Ill_posed_interpolation_example
     
     call Lebesgue_and_PI_functions
     
     call Chebyshev_polynomials
     call Interpolant_versus_Chebyshev 
   
end subroutine  



!*************************************************************************************
! It interpolates a function at a given abscisa x= xp for a given set of values:
!       (xj, fj) j=0, ..., N  
!******************************************************************************
subroutine Interpolated_value_example 

    integer, parameter :: N = 6 
    real :: x(N) = [ 0.0, 0.1, 0.2, 0.5, 0.6, 0.7 ] 
    real :: f(N) = [ 0.3, 0.5, 0.8, 0.2, 0.3, 0.6 ] 
    
    real :: xp = 0.15 
    real :: yp(N-1)
    integer :: i 
      
    do i=2, N-1 
      yp(i) = Interpolated_value( x , f , xp, i )
      write (*,'(A10, i4, A40, f10.3)') 'Order = ', i, 'The interpolated value at xp is = ', yp(i)
    end do 
    
    
    write (*,'(A20, 10f8.3)') 'xp = ', xp
    write (*,'(A20, 10f8.3)') 'nodes x_j = ', x
    write (*,'(A20, 10f8.3)') 'function f_j = = ', f
    write(*,*) "press enter " 
    read(*,*)
   
end subroutine
 
  

!********************************************************************
!*  It builds an interpolant for a given set of values 
!************************************************************
subroutine Interpolant_example

 integer, parameter :: N=3, M=400 
 real ::  xp(0:M) 
 real :: x(0:N) = [ 0.0, 0.1, 0.2, 0.5 ] 
 real :: f(0:N) = [ 0.3, 0.5, 0.8, 0.2 ]
 real :: I_N(0:N, 0:M)  ! first index: derivative
                        ! second index: point where the interpolant is evaluated 
 real :: a, b 
 integer :: i 
 
  a = x(0); b = x(N) 
  xp = [ (a + (b-a)*i/M, i=0, M) ] 

  I_N = Interpolant(x, f, N, xp) 
 
  write(*,*) "It plots an interpolant and its derivative"
  write(*,*) "press enter ";  read(*,*)
  call plot_parametrics(xp, transpose(I_N(0:1,:)),["I", "dI/dx"],"x","y")
   
 
end subroutine    
    
 



!***************************************************************************************
! Integral_example
!***************************************************************************************
subroutine Integral_example
 integer, parameter :: N=6
 real :: x(0:N), f(0:N), a = 0, b = 1, I0  
 integer :: i 
   x = [ (a + (b-a)*i/N, i=0, N) ] 
   f = sin ( PI * x )
   
   I0 = Integral( x, f, 4 )
   write(*, *) "The integral [0,1] of sin( PI x ) is: ", I0
   write(*, *) "Error = ",  ( 1 -cos(PI) )/PI - I0
   write(*, *) "press enter ";  read(*,*)
end subroutine


!**************************************************************************************************
! Lagrange_polynomial_example
!*****************************************************************************************
subroutine Lagrange_polynomial_example
 integer, parameter :: N=4, M=400 
 real :: x(0:N), xp(0:M), a=-1, b=1 
 real :: Lg(-1:N, 0:N, 0:M)   ! Lagrange polynomials 
                              ! -1:N (integral, lagrange, derivatives) 
                              !  0:N ( L_0(x), L_1(x),... L_N(x)     ) 
                              !  0:M ( points where L_j(x) is evaluated  ) 
 real :: Lebesgue_N(-1:N, 0:M) 
 character(len=2) :: legends(0:N) = [ "l0", "l1", "l2", "l3","l4" ]
 integer :: i 
 
 x  = [ (a + (b-a)*i/N, i=0, N) ]
 xp = [ (a + (b-a)*i/M, i=0, M) ]  
    
 do i=0, M;  Lg(:, :, i) = Lagrange_polynomials( x, xp(i) ); end do 
 Lebesgue_N = Lebesgue_functions( x, xp ) 
   
 write (*, *) "It plots Lagrange functions"
 write(*,*) "press enter ";  read(*,*)
 call plot_parametrics(xp, transpose(Lg(0, 0:N, :)), legends,"x","y")
end subroutine

 


!********************************************************************
!* Ill_posed_interpolation_example 
!***********************************************************
subroutine Ill_posed_interpolation_example 
 
 integer, parameter :: N=64, M=300
 real :: x(0:N), f(0:N)
 real :: I_N(0:N, 0:M)          
 real :: Lebesgue_N(-1:N, 0:M) 
 real :: xp(0:M)   
 real :: a=-1, b=1
 integer :: i  
 
 x  = [ (a + (b-a)*i/N, i=0, N) ] 
 xp = [ (a + (b-a)*i/M, i=0, M) ] 
 f = sin ( PI * x ) 
 
 I_N = Interpolant(x, f, N, xp) 
 Lebesgue_N = Lebesgue_functions( x, xp ) 
 
 write(*, *) "It plots an interpolant with errors at boundaries " 
 write(*, *) "maxval Lebesgue =", maxval( Lebesgue_N(0,:) )  
 
 write(*, *) "press enter "; read(*,*)
 call plot_parametrics(xp, transpose(I_N(0:0, :)), ["I"], "x", "y")
 
end subroutine 




!********************************************************************
!* Lebesgue and PI functions
!********************************************************************
subroutine Lebesgue_and_PI_functions 
 
 integer, parameter :: N=10, M=700
 real :: x(0:N), xp(0:M)
 real :: Lebesgue_N(-1:N, 0:M),  PI_N(0:N, 0:M) 
 real :: a=-1, b=1
 integer :: i, k  
 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 
   xp = [ (a + (b-a)*i/M, i=0, M) ] 
 
   Lebesgue_N = Lebesgue_functions( x, xp ) 
   PI_N = PI_error_polynomial( x, xp ) 
  
   write(*, *) "It plots Lebesgue  functions" 
   write(*, *) "function, first and second derivative" 
   write(*, *) "press enter ";  read(*, *)
   call plot_parametrics(xp, transpose(Lebesgue_N(0:2, :)),  & 
                         ["l", "dldx","d2ldx2"], "x","y")
   
   call plot_parametrics(xp, transpose(PI_N(0:2, :)),        & 
                         ["pi", "dpidx","d2pidx2"], "x","y")                       
  
 
end subroutine 



 
!**************************************************************************************************
!
!**************************************************************************************************
subroutine Chebyshev_polynomials
 
    integer, parameter :: N = 100, M = 5 
    real :: x(0:N), theta(0:N), Tk(0:N, 0:M), Uk(0:N,0:M)
    real :: x0=-1, xf= 1
    integer :: i, k
    character(len=2) :: lTK(0:M) = ["T0","T1","T2","T3","T4","T5"]
    character(len=2) :: lUk(0:M) = ["U0","U1","U2","U3","U4","U5"]
    
    x = [ ( x0 + (xf-x0)*i/N, i=0, N ) ] 
    do k=1, M
           theta = acos(x) 
           Tk(:, k) = cos ( k * theta )
           Uk(:, k) = sin ( k * theta ) / sin (theta)
    end do 
    
    write(*, *) "It plots Chebyshev polynomials" 
    write(*, *) "press enter ";  read(*,*)
    call plot_parametrics(x, Tk(:, 0:M), lTk, "x","y")
    call plot_parametrics(x, Uk(:, 0:M), lUk, "x","y")
end subroutine



!**************************************************************************************************
! Comparison between Interpolant and Chebyshev truncated series 
!**************************************************************************************************
subroutine Interpolant_versus_Chebyshev 
   
    integer, parameter :: N = 6  ! # of Chebyshev terms or poly order 
    integer, parameter :: M = 500! # of points to plot 
    real :: x(0:N), f(0:N)
    real :: I_N(0:N, 0:M)          ! Interpolant 
    real :: P_N(0:M)               ! Truncated series 
    real :: xp(0:M), theta(0:M)    ! domain to plot 
    real :: Error(0:M, 2)          ! Error :interpolant and truncated  
    real :: Integrand(0:M)  
    character(len=8) :: legends(2) =  ["Error_I", "Error_P"]
    
    integer :: i, k
    real :: c_k, a=-1, b = 1, gamma
    
!** equispaced points to plot 
    xp = [ (a + (b-a)*i/M, i=0, M) ] 
    theta = acos( xp ) 
  
!** Chebyshev truncated series    
    do k=0, N 
        Integrand = sin( PI * xp) * cos ( k * theta ) 
        
        if (k==0) then;    gamma =  PI; 
                  else;    gamma =  PI / 2; 
        end if 
        c_k = Integral( theta , Integrand ) / gamma
        P_N = P_N - c_k * cos( k * theta ) 
    end do 
   
    x = [ (cos(PI*i/N), i=N, 0, -1) ] 
    f = sin( PI * x )
 
!** Interpolant based on Chebyshev extrema     
    I_N = Interpolant(x, f, N, xp) 
    
    Error(:, 1) = sin( PI * xp ) - I_N(0, :) 
    Error(:, 2) = sin( PI * xp ) - P_N 
    call plot_parametrics(xp, Error(:, 1:2), legends, "x","y")
   
end subroutine



end module 
    
    
    
    