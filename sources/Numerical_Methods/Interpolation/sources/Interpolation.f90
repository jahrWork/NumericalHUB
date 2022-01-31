
!**********************************************************************************************************
!            Piecewise polynomial interpolation  and  to perform weighted averages                                                 
!
! Author:  Juan A. Hernandez   (juanantonio.hernandez@upm.es )             
!**********************************************************************************************************
module Interpolation

    use Lagrange_interpolation
    use Chebyshev_interpolation 
    use Fourier_interpolation 
    implicit none
    
  private 
  public ::           & 
  Interpolated_value, & ! It interpolates at xp from (x_i, y_i)
  Integral,           & ! It integrates from x_0 to x_N    
  Interpolant           ! It interpolates I(xp) from (x_i, y_i)

contains

!*************************************************************************************************************
!* Computes the piecewise polynomial interpolation of I(xp) from the data x(:) and y(:). 
!  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*************************************************************************************************************
real pure function Interpolated_value(x, y, xp, degree)
    real, intent(in) :: x(0:), y(0:), xp
    integer, optional, intent(in) :: degree
    integer :: N, s, j
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk !
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 

    N = size(x) - 1

    
    if(present(degree))then
        Nk = degree
    else
        Nk = 2
    end if

    allocate( Weights(-1:Nk, 0:Nk))
  
    j = max(0, maxloc(x, 1, x < xp ) - 1) 

    if( (j+1) <= N ) then !  the (j+1) cannot be accessed
        if( xp > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   


    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For Nk=2 
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

  
    Weights(-1:Nk, 0:Nk)  = Lagrange_polynomials( x = x(s:s+Nk), xp = xp )
    interpolated_value = sum ( Weights(0, 0:Nk) * y(s:s+Nk) ) 
 
    deallocate(Weights)

end function


!************************************************************************************************************
!* Computes the piecewise polynomial integration of y(x) from the data x(:) and y(:).
!*  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!**************************************************************************************************************
real function Integral(x, y, degree)
    real, intent(in) :: x(0:), y(0:)
    integer, optional, intent(in) :: degree

    integer :: N, j, s 
    real :: summation, Int, xp
    
!   maximum order of derivative and width of the stencil
    integer :: Nk 
    
 !  Lagrange coefficients and their derivatives at xp    
    real, allocatable :: Weights(:,:,:)

    N = size(x) - 1

    if(present(degree))then
                          Nk = degree
    else
                          Nk = 2
    end if


    allocate(Weights(-1:Nk, 0:Nk, 0:N))

    summation = 0   

    do j=0, N
     if (mod(Nk,2)==0) then 
                          s = max( 0, min(j-Nk/2, N-Nk) )   
     else 
                          s = max( 0, min(j-(Nk-1)/2, N-Nk) )
     endif 
     xp = x(j) 
     Weights(-1:Nk, 0:Nk, j)=Lagrange_polynomials(x = x(s:s+Nk), xp = xp ) 
 
     Int = sum ( Weights(-1, 0:Nk, j) * y(s:s+Nk) ) 

     summation   = summation  + Int 

   enddo 

   Integral = summation

   deallocate(Weights)

end function

!******************************************************************************************************************
!* Computes the average of the function f(x) weighted by the function w(x):  <f> = Int[f(x)w(x)dx] / Int[w(x)dx]. 
!*  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*******************************************************************************
real function weighted_average( f, w, x, degree )
    real, intent(in) :: f(0:), w(0:), x(0:)
    integer, optional, intent(in) :: degree


    if(present(degree))then
        weighted_average = Integral(x, f*w , degree) / Integral(x, w, degree)
    else
        weighted_average = Integral(x, f*w) / Integral(x, w)
    end if

end function


!*************************************************************************************************************
!* Computes the piecewise polynomial interpolation of I(xp) from the data x(:) and y(:). 
!  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*************************************************************************************************************
function Interpolant(x, y, degree, xp )
    real, intent(in) :: x(0:), y(0:), xp(0:) 
    integer, intent(in) :: degree
    real :: Interpolant(0:degree, 0:size(xp)-1)
    
    integer :: N, M, s, i, j, k
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk 
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 
    
    N = size(x) - 1
    M = size(xp) - 1 
    Nk = degree
    allocate( Weights(-1:Nk, 0:Nk))
    
do i=0, M 
        
    j = max(0, maxloc(x, 1, x < xp(i) ) - 1) 

    if( (j+1) <= N ) then ! the (j+1) cannot be accessed
        if( xp(i) > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   

    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For Nk=2 
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

  
    Weights(-1:Nk, 0:Nk)=Lagrange_polynomials( x = x(s:s+Nk), xp = xp(i) )
     
    do k=0, Nk 
        Interpolant(k, i) = sum ( Weights(k, 0:Nk) * y(s:s+Nk) ) 
    end do    
     
end do  
 
    deallocate(Weights)

end function




end module
