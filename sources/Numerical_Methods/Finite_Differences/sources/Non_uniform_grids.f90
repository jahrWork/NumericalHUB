module Non_uniform_grids 


 use Lagrange_interpolation
 use Non_Linear_Systems


 implicit none 

    contains 



!********************************************************************************************************
!   This module allows to generate  non uniform grids. 
!   These routines arecalled by the high order finite diference schemes (routines) 
!
! Author : Juan A Hernandez (juanantonio.hernandez@upm.es) 
!*********************************************************************************************************
subroutine Uniform_grid(x) 
      real, intent(out)   :: x(0:)
      
 
  real :: a = -1, b = 1;  
  integer :: N, i
  real :: x0, xf 

    N = size(x)-1 
    x0 = x(0) 
    xf = x(N) 
    
    do i = 0, N
        x(i) = a  + (b-a) *  i / N
    enddo
     x = [ ((x0+xf)*0.5 + (xf-x0)*0.5*x(i), i=0, N) ]

end subroutine 
    


!************************************************************************************
! Unkonwn nodes y(1:N) are allowed to form a polynomial of degree N-1
! Its polynomial error is:  
!                           PI(x) = (x-y1)*(x-y2)*(x-y3)... (x_yN) 
!
! Extrema of PI(x) ( dpi(x)/dx = 0) are obtained to give: x1, x2,.... x_n-1
!
! These nodes are completed with x_0 = a  and x_n = b to form the non uniform grid 
!
! Important: to assure convergece MUST be a = -1 and b = 1 
!************************************************************************************
subroutine Non_uniform_grid( x, Order) 
      real, intent(out)   :: x(0:)
      integer, intent(in) :: Order
             
  
  real :: y(1:size(x)-1), Errors(0:size(x)-1, 0:2)  
  real :: a = -1, b = 1;  
  integer :: N, i, Degree
   real :: x0, xf 


    N = size(x)-1 
    x0 = x(0) 
    xf = x(N)
    write(*,*) " Determining non uniform grid.. " 

!   initial guess for the unknowns 
    do i = 1, N
        y(i) = a  + (b-a) *  i / N
    enddo

!   Solves increasing polynomial degrees to ensure convergence
    do Degree = 2, Order-1
       call Newton( Extrema_equations, y ) 
       write(*,*) " Order = ", Degree 
    enddo

!   Computes the dual nodes of the Gauss nodes, which are the interesting ones
    call Extrema_points(Order-1, y, x, Errors )
    
    x = [ ((x0+xf)*0.5 + (xf-x0)*0.5*x(i), i=0, N) ]
    

contains
!************************************************************************************
!* Pi(x_0) = Pi(x_1) = Pi(x_2) ..... = Pi(x_n) 
!************************************************************************************
 function  Extrema_equations(y)
        real, intent(in)  :: y(:)
        real :: Extrema_equations( size(y) )

     integer :: N, i
     real    :: x_extrema(0:size(y)), Extrema_Errors(0:size(y), 0:2), E( size(y) )            

   !  write(*,*) "Extrema "
     N = size(y)

     call Extrema_points(Degree, y, x_extrema, Extrema_Errors)

     do i = 1, N
        E(i) = abs(Extrema_Errors(i,0)) - abs(Extrema_Errors(i-1,0))
     enddo
     
     Extrema_equations = E 

end function

!************************************************************************************
! Extrema of Pi(x) <=> dPi/dx = 0 
! Newton =>  x_i = x_i - dPi/dx / d2Pi/dx2, i=1... N-1  
!************************************************************************************
subroutine Extrema_points(Degree, y, x_extrema, Extrema_Errors)
  integer, intent(in) :: Degree
  real, intent(in)    :: y(:)
  real, intent(out)   :: x_extrema(0:), Extrema_Errors(0:,0:)


  integer :: N, i, s, q, Seeds( size(y) )  
  real    :: dE, Errors(0:size(y))


    N = size(y) 
    call Generate_seeds( Degree, Seeds ) 
    q = Degree;

    x_extrema(0) = a
    x_extrema(N) = b
    do i=1, N-1
        x_extrema(i) = ( y(i+1) + y(i) )/2 
    end do 
    


!** extrema points from unkonwn points y1,y2, yN
    do i = 1, N-1 
                  s = Seeds(i)
                  dE = 1
                  do while (abs(dE) > 1d-8)
                           Errors(0:q+1) =  Lagrange_Error_Polynomial( y(s:s+q), x_extrema(i) )
                           dE = Errors(1) / Errors(2)
                          x_extrema(i) = x_extrema(i) - dE 
                  enddo   
                  Extrema_Errors(i,:) = Errors(0:2)
    enddo

    Extrema_Errors(0,:)   =  Lagrange_Error_Polynomial( y(1:1+q), x_extrema(0) )
    Extrema_Errors(N,:)   =  Lagrange_Error_Polynomial( y(N-q:N), x_extrema(N) )


end subroutine 


end subroutine 


!******************************************************************************
! Different seedss to descenter the stencil close to de boundaries
! Centered schemes at inner points. 
!******************************************************************************
subroutine  Generate_seeds( Degree, Seeds ) 
     integer, intent(in) :: Degree
     integer, intent(out) ::  Seeds(:)


  integer :: N, i, N1, N2 



    N = size(Seeds) 
   

    if (mod(Degree,2)==0) then 
                              N1 = Degree/2;    N2 = Degree/2; 
                  
    else 
                              N1 = (Degree-1)/2; N2 = (Degree+1)/2; 
    endif   
               
    Seeds(1:N1)      = 1;   
    Seeds(N1+1:N-N2) = [ ( i+1, i=0, N-N1-N2 ) ];  
    Seeds(N-N2+1:N)  = N - Degree ;

end subroutine 
 

end module 


