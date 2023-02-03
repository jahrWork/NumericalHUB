module API_Example_Systems_of_Equations

    use Linear_systems
    use Non_Linear_Systems 
    use Numerical_Recipes
    use Utilities

    implicit none


    
    contains    
    
subroutine Systems_of_Equations_examples
    
    call LU_Solution
    call Newton_Solution
    call Implicit_explicit_equations
    call Test_Power_Method 
    call Test_eigenvalues_PM  
    call Vandermonde_condition_number
    call Test_Linear_regression
    
end subroutine   
  
!****************************************************************
! 
!****************************************************************
subroutine LU_Solution
    real :: A(4,4), b(4), x(4)
    integer :: i 
   
    A(1,:) = [ 4, 3, 6, 9]
    A(2,:) = [ 2, 5,  4, 2]
    A(3,:) = [ 1, 3, 2, 7]
    A(4,:) = [ 2, 4, 3, 8]
    b = [ 3, 1, 5, 2]
    
    write (*,*) 'Linear system of equations  '
    write (*,*) 'Matrix of the system A= '
    do i=1, 4;   write(*,'(100f8.3)') A(i, :);  end do
    write (*,'(A20, 100f8.3)') 'Independent term b=  ', b 
    write(*,*)  
    
    call LU_factorization( A )    
    x = Solve_LU( A , b )
    
    write (*,'(A20, 100f8.3)') 'The solution is = ', x
    write(*,*) "press enter " ; read(*,*)
end subroutine 


!****************************************************************
! 
!****************************************************************
subroutine Newton_Solution

    real :: x0(3) = [1., 1., 1.  ] 

    call Newton( F, x0 ) 
    
    write(*,*)  'Zeros of F(x) by Newton method '
    write(*,*)  'F(1) = x**2 - y**3 - 2'
    write(*,*)  'F(2) = 3 * x * y - z'
    write(*,*)  'F(3) = z**2 - x'
    write(*,'(A20, 100f8.3)')  'Zeroes of F(x) are x = ', x0
    write(*,*) "press enter " 
    read(*,*)
   
end subroutine 


function F(xv)
    real, intent(in) :: xv(:)
    real:: F(size(xv))
 
    real :: x, y, z
 
    x = xv(1); y = xv(2); z = xv(3)
        
    F(1) = x**2 - y**3 - 2
    F(2) = 3 * x * y - z
    F(3) = z**2 - x
   
end function

!******************************************************
!* Test of NewtonC : Implicit and explicit equations 
!     F1 : three implicit equations 
!     F2 : two implicit equations + one explicit 
!*******************************************************
subroutine Implicit_explicit_equations

    real :: x(3)
    write(*,*)  'Zeros of F(x) by Newton method '
    write(*,*)  'F(1) = x**2 - y**3 - 2'
    write(*,*)  'F(2) = 3 * x * y - z'
    write(*,*)  'F(3) = z**2 - x'
    
       x = 1 
       call Newtonc( F = F1, x0 = x  )
       write(*,'(A35, 100f8.3)')  'Three implicit equations, x = ', x
       
       x = 1 
       call Newtonc( F = F2, x0 = x  ) 
       write(*,'(A35, 100f8.3)')  'Two implicit + one explicit, x = ', x
       write(*,*) "press enter "; read(*,*) 
contains  
    
function F1(xv) result(F)
    real, target :: xv(:)
    real :: F(size(xv))
 
    real, pointer :: x, y, z
    
    x => xv(1); y => xv(2); z => xv(3)
        
    F(1) = x**2 - y**3 - 2
    F(2) = 3 * x * y - z
    F(3) = z**2 - x
   
end function 

function F2(xv) result(F)
    real, target :: xv(:)
    real :: F(size(xv))
 
    real, pointer :: x, y, z
    
    x => xv(1); y => xv(2); z => xv(3)
    
    x =  z**2 
    
    F(1) = 0 ! forall xv 
    F(2) = x**2 - y**3 - 2
    F(3) = 3 * x * y - z
   
end function 

end subroutine    
    

!****************************************************************
! Power method  
!****************************************************************
subroutine Test_Power_Method 

      integer, parameter :: N = 3 
      real :: A(N, N), lambda, U(N)=1   
      
      A(1,:) = [ 7, 4, 1 ] 
      A(2,:) = [ 4, 4, 4 ] 
      A(3,:) = [ 1, 4, 7 ] 
      
      write (*,*) 'Power method  '
      call power_method(A, lambda, U) 
      
      write(*,'(A10, f8.3, A10, 3f8.5)')  "lambda= ", lambda, "U= ", U
      write(*,*) "press enter "; read(*,*) 
     
end subroutine 

!****************************************************************
! Test eigenvalues 
!****************************************************************
subroutine Test_eigenvalues_PM 

 integer, parameter :: N = 3 
 real :: A(N, N), lambda(N), U(N, N) 
 integer :: i 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 4 ] 
 A(3,:) = [ 1, 4, 7 ] 
  
 call Eigenvalues_PM(A, lambda, U) 
 
 do i=1, N 
    write(*,'(A8, f8.3, A15, 3f8.3)')  & 
         "lambda = ", lambda(i),  "eigenvector = ", U(:,i)
 end do   
  
end subroutine 


!****************************************************************
! Test eigenvalues 
!****************************************************************
subroutine   test_SVD2 

 integer :: i, j, k  
 integer, parameter :: N = 3 
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N)  
 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 4 ] 
 A(3,:) = [ 1, 4, 7 ] 
 
 do i=1, N 
    write(*,'(A8, 3f8.3)') "A = ", A(i,:) 
 end do  
 
 call SVD(A, sigma, U, V) 

 do i=1, N 
    write(*,'(A8, f8.3, A15, 3f8.3, A15, 3f8.3)') "sigma = ", sigma(i), "U = ", U(i, :), " V = ",  V(i, :)
 end do   
  
end subroutine
!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_Vandermonde

 integer :: i, j, k  
 integer, parameter :: N = 20
 real :: A(N, N), sigma(N), U(N, N), V(N,N)      
 
 
 
 do i=1, N; do j=1, N
    A(i,j) = (i/real(N))**j 
 end do; end do 

 call SVD(A, sigma, U, V) 
  
 do i=1, N 
    write(*,'(A8, e15.7)') "sigma = ", sigma(i)
 end do 
 
 
end subroutine

!****************************************************************
! Vandermonde_condition_number
!****************************************************************
subroutine Vandermonde_condition_number
   
    integer, parameter :: N = 10
    real :: A(N, N), kappa
    integer :: i, j  
   
    do i=1, N; do j=1, N; 
        A(i,j) = (i/real(N))**(j-1) 
    end do; end do 
    
    kappa = Condition_number(A) 
    
    write(*,*)  'Condition number of Vandermonde matrix '
    write(*,'(A40, e10.3)') " Condition number (power method) =", kappa
    write(*,*) "press enter "; read(*,*) 
end subroutine

!*********************************************************************************
!*
!*********************************************************************************
subroutine print_matrix( legend, A) 
character(len=*), intent(in) :: legend 
real, intent(in) :: A(:,:) 
  
integer :: N, i 

N = size( A, dim=1) 

write(*,*) legend 
do i=1, N 
     write(*,'(100f8.3)')  A(i,:) 
end do



end subroutine 


!*********************************************************************************
!*
!*********************************************************************************
subroutine print_eigenvectors( legend, U) 
character(len=*), intent(in) :: legend 
real, intent(in) :: U(:,:) 
  
integer :: N, i 

N = size( U, dim=1) 

write(*,*) legend 
do i=1, N 
     write(*,'(100f8.3)')  U(i,:) 
end do

do i=1, N-1 
     write(*,'(A, i5, a, 100ES29.6)') " i=", i, " uiT ui+1 =", dot_product(U(:, i), U(:, i+1) )  
end do



end subroutine 


!*********************************************************************************
!*
!*********************************************************************************
subroutine print_SVD( legend, U, S, V) 
character(len=*), intent(in) :: legend 
real, intent(in) :: U(:,:), S(:,:), V(:,:)  
  
integer :: N, M, i, j  

M = size( U, dim=1) 
N = size( V, dim=2) 

write(*,*) " M = ", M, " N = ", N 
write(*,*) legend 
 

do i=1, N 
    if (i>M) then 
        
         write(*, '(A1, A32, A1, A40, A1, 5f8.2, A1 )' )  "|", "                ",  "|" , "                ", "|", V(:,i), "|" 
    
    else
 
       write(*, '(A1, 4f8.2, A1, 5f8.2, A1, 5f8.2, A1 )' )  "|", U(i,1:M),  "|" , S(i,1:N), "|", V(:,i), "|" 
   end if 
end do
write(*,*)


end subroutine 



!****************************************************************
! Vandermonde_condition_number
!****************************************************************
subroutine Vandermonde_SVD_condition_number
   
    integer :: i, j, k  
    integer, parameter :: N = 10
    real :: A(N, N), U(N, N), S(N, N), V(N, N)
    real :: A_SVD(N, N), sigma(N)
    real :: kappa 

    do i=1, N; do j=1, N; 
        A(i,j) = (i/real(N))**(j-1) 
    end do; end do 
    
    call SVD(A, sigma, U, V) 
    S = 0 
    do i=1, N;  S(i,i) = sigma(i); end do 
    A_SVD = matmul(U , matmul( S ,  transpose(V) ) )
    
    kappa = Condition_number(A) 
    
    write(*,*)  'Condition number of Vandermonde matrix '
    write(*,'(A40, e10.3)') " Condition number (SVD) =", maxval(sigma)/minval(sigma) 
    write(*,'(A40, e10.3)') " Condition number (power method) =", kappa
    write(*,*) "press enter "; read(*,*) 
end subroutine

    
!****************************************************************
! Test SVD
!****************************************************************
subroutine Test_SVD 

 integer, parameter :: M = 4, N = 5
 real :: A(M, N), B(M,N), sigma(N), U(M, M), V(N, N), S(M,N) 
 real :: sigma_min = 1e-5
 integer :: i
 A(1,:) = [ 1, 0, 0, 0, 2 ] 
 A(2,:) = [ 0, 0, 3, 0, 0 ] 
 A(3,:) = [ 0, 0, 0, 0, 0 ]
 A(4,:) = [ 0, 2, 0, 0, 0 ]
  
 write(*,*) "-----Test SVD decomposition----(press enter) "; read(*,*) 
 call SVD(A, sigma, U, V, sigma_min) 
 call print_matrix(" matrix A =", A )
 call print_eigenvectors(" check U eigenvectors =", U )
 call print_eigenvectors(" check V eigenvectors =", V )
 
 B = 0 
 do i=1, N-1 
     B = B + sigma(i) * Tensor_product( U(:,i), V(:,i) ) 
 end do 
 call print_matrix(" approximated matrix A =", B )
 
 S = 0; do i=1, M;  S(i,i) = sigma(i);  end do
 call print_SVD(" SVD  A = U S VT = ", U, S, V  )
 
end subroutine    
    
    
!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_SVD1 

 integer :: i, j, k  
 integer, parameter :: N = 3
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N), A_SVD(N,N), D(N,N), C(N,N), B(N,N)       
 
 
 A(1,:) = [ 7, 4, 1 ] 
 A(2,:) = [ 4, 4, 5 ] 
 A(3,:) = [ 1, 4, 7 ]
  
 
 call SVD(A, sigma, U, V) 
 
 D = 0 
 do i=1, N 
     D(i,i) = sigma(i) 
 end do 
 
 do i=1, N 
    write(*,'(A8, 100f8.3)') "A = ", A(i,:) 
 end do  
 
 do i=1, N 
    write(*,'(A8, f8.3)') "sigma = ", sigma(i)
 end do 
 
 !write(*,*) 
 ! do i=1, N 
 !    do j=1, N 
 !        write(*,'(A5, 2i3, 100f8.3)') "V dot V  = ", i, j, dot_product(V(:,i), V(:,j) )  
 !    end do 
 !end do
 
 C = matmul( V , matmul( matmul(D,D) , transpose(V) ) )
 do i=1, N 
    write(*,'(A20, 100f8.3)') "transpose(A) * A  = ", C(i,:) 
 end do 
 
 write(*,*) 
 B = matmul( transpose(A), A )
  do i=1, N 
    write(*,'(A20, 100f8.3)') "transpose(A) * A  = ", B(i,:) 
  end do 
 
 !write(*,*) 
 ! do i=1, N 
 !    do j=1, N 
 !        write(*,'(A5, 2i3, 100f8.3)') "U dot U  = ", i, j, dot_product(U(:,i), U(:,j) )  
 !    end do 
 ! end do
  
  !write(*,*) " B ="
  !B = matmul( transpose(U), U ) 
  !do i=1, N 
  !       write(*,'(A15, 100f8.3)') "transpose(U)  U  = ",  B(i,:)  
  !end do
  !
  !write(*,*) " B ="
  !B = matmul( U, transpose(U) ) 
  !do i=1, N 
  !       write(*,'(A15, 100f8.3)') " U transpose(U)   = ",  B(i,:)  
  !end do
  
  
 write(*,*) 
 B = matmul( A, transpose(A) )
  do i=1, N 
    write(*,'(A20, 100f8.3)') "A * transpose(A)  = ", B(i,:) 
  end do 
  
 write(*,*)  
 C = matmul( U , matmul( matmul(D,D) , transpose(U) ) )
 do i=1, N 
    write(*,'(A20, 100f8.3)') "A * transpose(A)  = ", C(i,:) 
 end do 
 
 
 
 !A_SVD = matmul( U , matmul( D ,  transpose(V) ) )
 
 A_SVD = matmul( transpose(U) , matmul( A ,  V ) )

 write(*,*) 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "Sigma = ", A_SVD(i, :)
 end do 
 

 A_SVD = matmul(U , matmul( D ,  transpose(V) ) )

 write(*,*) 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "A_SVD = ", A_SVD(i, :)
 end do 
 
 
 
end subroutine

!****************************************************************
! Test SVD
!****************************************************************
subroutine   test_SVD3 

 integer :: i, j, k  
 integer, parameter :: N = 2
 real :: A(N, N)
 real :: sigma(N), U(N, N), V(N,N), A_SVD(N,N), D(N,N), ATA(N,N), AAT(N,N)      
 
 
 A(1,:) = [ 0,  1 ] 
 A(2,:) = [ -1, 0 ] 
 
 
 do i=1, N 
    write(*,'(A8, 2f8.3)') "A = ", A(i,:) 
 end do  
 
 call SVD(A, sigma, U, V) 
 
 do i=1, N 
    write(*,'(A8, f8.3, A15, 2f8.3, A15, 2f8.3, A5, 2f8.3 )') "sigma = ", sigma(i), "U = ", U(:,i), " V = ",  V(:,i), " norm =", norm2(U(:,i)), norm2(V(:,i)) 
 end do   
 
 D = 0 
 do i=1, N 
     D(i,i) = sigma(i)
 end do 
 
  
 A_SVD = matmul( U , matmul( D , transpose(V) ) )
 ATA = matmul( V , matmul( D , transpose(V) ) )
 AAT = matmul( U , matmul( D , transpose(U) ) )
 
 do i=1, N 
    write(*,'(A5, 100f8.3)') "U = ", U(i,:) 
 end do
 do i=1, N 
    write(*,'(A5, 100f8.3)') "V = ", V(i,:) 
 end do
 
 
 do i=1, N 
    write(*,'(A20, 100f10.5)') "transpose(A) * A  = ", ATA(i,:) 
 end do 
 
 
  
 do i=1, N 
    write(*,'(A20, 100f10.5)') "A * transpose(A)   = ", AAT(i,:) 
 end do 
 
 do i=1, N 
    write(*,'(A8, 100f10.5)') "A_SVD = ", A_SVD(i, :)
 end do 
 
 
end subroutine


!****************************************************************
! Test SVD
!****************************************************************
subroutine Test_Linear_regression  

 integer, parameter :: N = 6 
 real :: x(N)  = [5, 15, 25, 35, 45, 55]
 real :: y(N) =  [5, 20, 14, 32, 22, 38]
 real :: coef(2) 
! y = 0.54 x +  5.6333
 
 write(*,*) "--------- Linear regression----------- " 
 write(*,*) "given a set of points (xi, yi)" 
 write(*,*) "obtain the best fit y = c1 x + c0 " 
 write(*,*) " press enter"
 read(*,*) 
 coef = linear_regression(x, y) 
 write(*,*) " coef =", coef 
 write(*,*)
 
end subroutine   
 




!*********************************************************
!* Example to test relative max and min
!********************************************************* 
subroutine Test_max_min 

      real :: x(2) 
   
     x = 0 
     call max_min(Energy, x) 
     write(*,*) " Solution =", x 

contains
real function Energy(x) 
      real, intent(in) :: x(:) 
 
      real :: x1, x2 
      
      x1 = x(1) 
      x2 = x(2) 
 
      Energy = x1**2 - x2**2 - 2*x1*x2 + 2*x1 + 2*x2 
 
 
end function    

end subroutine

    
!*********************************************************
!* Example to test relative max and min
!********************************************************* 
subroutine Test_gradient_descent 

      real :: x(2) 
   
     x = 0  
     call Gradient_descent(Energy, x) 
     write(*,*) " Solution =", x 
     
contains  



end subroutine

real function Energy(x) 
      real, intent(in) :: x(:) 
 
      real :: x1, x2 
      
      x1 = x(1) 
      x2 = x(2) 
 
      Energy = x1**2 - x2**2 - 2*x1*x2 + 2*x1 + 2*x2 
 
 
end function    

!*********************************************************
!* Example to test the tensor product 
!********************************************************* 
subroutine Test_tensor_product  

      real :: u(2) = [1,2], v(2)= [3,4]  
      real :: A(2,2) 
      integer :: i 
   
    
     A = u .x. v 
     do i=1, 2 
      write(*,*) " A =", A(i,:)  
     end do 
     
 


end subroutine

 
end module 
