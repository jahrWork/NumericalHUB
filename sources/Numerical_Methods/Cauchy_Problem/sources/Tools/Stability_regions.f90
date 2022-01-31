module Stability_regions 
   
    use Temporal_scheme_interface 
    implicit none
    
private 
public :: Absolute_Stability_Region ! For a generic temporal scheme

contains

!***************************************************************************
!  Stability region                                                    
!  It determines the Regions absolute stability
!  of  numerical schemes 
!   
!  INPUTS: 
!           Scheme
!           x, y : partition of the complex domain z
!  OUTPUTS: 
!           Region : region of absolute stability 
!    
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 2019
!         Francisco Javier Escoto Lopez
!*************************************************************************
subroutine Absolute_Stability_Region(Scheme, x, y, Region)

procedure(Temporal_Scheme) :: Scheme
real, intent(in):: x(0:), y(0:)
real, intent(out):: Region(0:,0:)

integer :: i, j, k, l, ierr
integer :: Nx, Ny, steps
complex :: img_unity=(0.d0, 1.d0), w
real    :: t1 = 0.d0, t2=1d0, pi,escala
complex, allocatable :: roots(:), A(:)
real :: U1(2), U2(2)



Nx = size(x)-1
Ny = size(y)-1

! number of steps of the temporal scheme 
call Steps_determination(steps, Scheme) 
write(*,*) "# of steps = ", steps

! coefficients of the characteristic stability polynomial 
! roots of the polynomial for some point of the complex plane w = lambda dt = (x_i, y_j )  
allocate ( A(0:steps), roots(steps) )

pi = dacos(-1.d0)
A = -1.d0
do i = 0, Nx
    do j = 0, Ny
        
        w = cmplx ( x(i), y(j) )  
        do k = 1, steps
            t1 = 0
            do l = 1, steps

                if (k == l) then                        
                    U1(1) = 1
                    U1(2) = 0.0d0
                else
                    U1 = 0.d0
                endif
                
                call Scheme(ODE_RA, t1, t1+1, U1, U2, ierr)
                t1 = t1 + 1
            enddo
            A(k-1) = cmplx ( U2(1), U2(2) ) 
        enddo
    
        !*** Roots  initialization    
        do k=1, steps     
            roots(k) = exp( 2*pi*(k-1)*img_unity / steps )
        end do     
                    
        call Polynomial_roots(A, roots)
            
        Region(i, j) = maxval(abs(roots))

        
    end do   
end do


contains

   function ODE_RA(U,t) result(F)
           real :: U(:), t
           real :: F(size(U))
   
      complex ::  Fc
   
      call ODE_RAc( cmplx(U(1), U(2)), Fc ) 
   
      F(1:2) = [ real(Fc), imag(Fc) ]
   
   end function
   
   subroutine ODE_RAc(U, F)
         complex, intent(in)  ::  U
         complex, intent(out) :: F 
   
       F = w * U 
   
   end subroutine
   
end subroutine







subroutine Steps_determination(steps, Scheme)
   procedure(Temporal_Scheme) :: Scheme
   integer, intent(out)    :: steps

   real :: U1(1), U2(1), pi = acos(-1d0)
   real    :: t0 = 0.d0, t1, h = 1d-1
   integer :: i, ierr
   

    i = 0
    
    U1 = 1.d0


do      
    t1 = t0 + h
    call Scheme(Linear_ODE, t0, t1, U1, U2, ierr)
    
    !write(*,*) "U2=", U2
    if ( abs(U2(1)) <= epsilon(1d0) ) then

        steps = i
        exit

    endif

    U1 = 0
    t0 = t1
    i = i+1
    
end do

contains
    function Linear_ODE(U,t) result(F)  
        real ::  U(:), t
        real :: F(size(U))

        
            F(1) = U(1)

    end function

end subroutine




 subroutine Polynomial_reduction(A, b)

     complex, intent(inout)  :: A(:), b
     
     complex :: C, D
     integer :: n, i

     n = size(A)

     C = A(n-1) + A(n)*b
     A(n-1) = A(n)

     do i = n-2,1,-1

         D = A(i) + C*b
         A(i) =  C
         C = D

     enddo

 end subroutine

 subroutine Polynomial_computation(A,x,y)

     complex, intent(in) :: A(:), x
     complex, intent(out):: y

     integer :: n,i

     n = size(A)

     y = 0.d0

     do i = 1,n

         y = y + A(i) * (x**(i-1.d0))

     enddo

 end subroutine

 subroutine Polynomial_derivative(A,x,dydx)

     complex, intent(in) :: A(:), x
     complex, intent(out):: dydx

     integer :: n, i

     n = size(A)

     dydx = 0.d0

     do i = 1,n

         dydx = dydx + (i-1.d0) * A(i) * (x**(i-2.d0))

     enddo

 endsubroutine

 subroutine Polynomial_roots(A, roots)

     complex, intent(inout) :: A(:)
     complex, intent(inout):: roots(:)

     integer :: n, i
     
     n = size(A)

     do i = 1, n-1

         call Root_solution(A(1:n+1-i), roots(i))
         call Polynomial_reduction(A(1:n+1-i), roots(i))

     enddo

 end subroutine

 subroutine Root_solution(A,root)

     complex, intent(in) :: A(:)
     complex, intent(inout):: root


     complex :: y, dydx
     integer :: n, i, itmax = 2000
     
     n = size(A)

     call Polynomial_computation(A,root,y)

     do i=1, itmax

         call Polynomial_derivative(A,root,dydx)

         root = root - y/dydx

         call Polynomial_computation(A,root,y)

         if (abs(y) < 1.d-9) exit

         if (i==itmax) then

             write(*,*) "Maximum iteration number achieved:", i
             

         endif

     enddo

 endsubroutine




end module
