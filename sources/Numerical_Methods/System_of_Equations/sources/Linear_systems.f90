module Linear_systems 

implicit none 
private 
public ::               & 
  LU_factorization,     & ! A  = L U (lower, upper triangle matrices)
  Solve_LU,             & ! It solves L U x = b 
  Inverse,              & ! Inverse of A 
  Gauss,                & ! It solves A x = b by Guass elimination
  Condition_number,     & ! Kappa(A)  = norm2(A)  * norm2( inverse(A) )
  Tensor_product,       & ! A_ij = u_i v_j 
  Power_method,         & ! It determines to largest eigenvalue of A 
  Inverse_Power_method, & ! It determines the smallest eigenvalue of A
  Eigenvalues_PM,       & ! All eigenvalue of A by the power method
  SVD,                  & ! A = U S transpose(V)  
  linear_regression       ! y = c1 x + c0 
public :: operator  (.x.)

interface operator(.x.) 
       module procedure Tensor_product 
end interface 


contains 
!***************************************************************************
!*   Factorization of A by means of upper and lower triangle matrices 
!*          input  : A matrix 
!*          output : A  holds L U 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************
subroutine LU_factorization( A ) 
      real, intent(inout) :: A(:, :) 

  integer :: N
  integer :: k, i, j 
  N =size(A, dim =1) 
 
  A(1, :) = A(1,:) 
  A(2:N,1) = A(2:N,1)/A(1,1) 

  do k=2, N 
   do j=k, N
       A(k,j) =  A(k,j) - dot_product( A(k, 1:k-1), A(1:k-1, j) )
   end do   
   
   do i=k+1, N
       A(i,k) = (A(i,k) - dot_product( A(1:k-1, k), A(i, 1:k-1) ) )/A(k,k) 
   end do 
  end do 
end subroutine

!***************************************************************************
!* LU x = b 
!***************************************************************************
function Solve_LU( A, b ) 
  real, intent(in) :: A(:, :), b(:) 
  real :: Solve_LU( size(b) )  

   real :: y (size(b)), x(size(b))
   integer :: i, N
  
   N = size(b) 

   y(1) = b(1)
   do i=2,N
          y(i) = b(i) - dot_product( A(i, 1:i-1), y(1:i-1) )
   enddo

    x(N) = y(N) / A(N,N)
    do i=N-1, 1, -1
     x(i) = (y(i) - dot_product( A(i, i+1:N), x(i+1:N) ) )/ A(i,i)
    end do

    Solve_LU = x 

end function 

!***************************************************************************
!* It solves A B  = I with A = L U  and B the inverse of A 
!***************************************************************************
function Inverse(A) result(B)
  real, intent(inout) :: A(:, :) 
  real :: B( size(A,dim=1), size(A,dim=1) )  

   real :: one( size(A,dim=1) ) 
   integer :: i, N
  
   N = size(one) 
   
   call LU_factorization(A)
   
   do i=1, N 
       one = 0
       one(i) = 1 
       
       B(:, i) =  Solve_LU(A, one ) 
   end do 

end function 

!**********************************************************
! Solutions to a system of linear equations A*x=b
! Method: Gauss elimination (with scaling and pivoting)
! Alex G. (November 2009)
! Juan A. Hernandez, juanantonio.hernandez@upm.es (Feb 2021)   
!***********************************************************
function  Gauss(A, b) result(x) 
   real, intent(inout) :: A(:,:), b(:) 
   real :: x(size(b)) 

   integer :: N, i, j, k, l  
   real :: S(size(b)), c, pivot 
   N = size(b) 
 

! begin forward elimination
  do k = 1, N-1
      
!    s(i) will have the largest element from row i
     do i = k, N ! loop over rows
       S(i) = maxval( abs(A(i,k:N)) ) 
     end do
      
!   find a row with the largest pivoting element
     pivot = abs( A(k,k) / s(k) )
     l = k
     do j = k+1, N 
       if( abs(A(j,k) / s(j) ) > pivot) then
          pivot = abs( A(j,k) / s(j) )
          l = j
       end if
     end do

!    Check if the system has a sigular matrix
     if(pivot == 0.0) then
       write(*,*) " The matrix is sigular "
       stop 
     end if

!    pivoting: swap rows k and l (if needed)
     if (l /= k) then
       call Swap( A(k,k:N), A(l,k:N) ) 
       call Swap( b(k:k), b(l:l) )  
     end if

!    the elimination (after scaling and pivoting)
     do i = k + 1, N  ! for all rows below pivot:
       c = A(i, k) / A(k,k)  
       A(i, 1:k) = 0 
       do j = k+1, N
           A(i, j) = A(i, j) - c * A(k, j) 
       end do     
       b(i) = b(i) - c * b(k) 
     enddo 
    
  end do
  
  
! back substiturion
  x(N) = b(N) / A(N,N) 
  do i = N-1, 1, -1  
     x(i) = ( b(i) - dot_product( A(i,i+1:N), x(i+1:N) ) )/ A(i,i) 
  enddo 
 
end function



!***************************************************************************
!* Swap two vectors 
!***************************************************************************
subroutine Swap( A, B ) 
  real, intent(inout) :: A(:), B(:)  

  real :: temp( size(A) ) 

   temp = A
   A = B
   B = temp 

end subroutine






!****************************************************************
!* Power method  to obtain the largest eigenvalue of a square matrix A 
!*          input  : A matrix 
!*
!*          output : lambda eigenvalue 
!*                   U eigenvector 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!****************************************************************
subroutine Power_method(A, lambda, U) 
      real, intent(in) :: A(:,:)
      real, intent(out) :: lambda, U(:) 
  
  integer :: N, k, k_max = 10000 
  real, allocatable :: U0(:), V(:)  
  
    N = size( A, dim=1) 
    allocate( U0(N), V(N) ) 
    U = [ (k, k=1, N) ] 
    
    k = 1 
    do while( norm2(U-U0) > 1d-12 .and. k < k_max )  
         U0 = U
         V = matmul( A, U ) 
         U = V / norm2(V) 
         k = k + 1         
    end do 
    
    lambda = dot_product( U, matmul(A, U) )    
    
    
end subroutine 

!****************************************************************
! Inverse Power method to obtain the smallest eigenvalue of a square matrix A 
!*          input  : A matrix 
!*
!*          output : lambda eigenvalue 
!*                   U eigenvector 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!****************************************************************
subroutine Inverse_power_method(A, lambda, U) 
      real, intent(inout) :: A(:,:)
      real, intent(out) :: lambda, U(:) 
  
  integer :: N, k, k_max = 10000 
  real, allocatable :: U0(:), V(:), Ac(:, :)   
    N = size(U) 
    allocate ( Ac(N,N), U0(N), V(N) ) 
    Ac = A 
    call LU_factorization(Ac) 
    U = [ (k, k=1, N) ] 

    k = 1 
    do while( norm2(U-U0) > 1d-12 .and. k < k_max )  
         U0 = U 
         V = solve_LU(Ac, U) 
         U = V / norm2(V) 
         k = k + 1         
    end do     
    lambda = norm2(matmul(A, U))    
    
end subroutine 

!****************************************************************
!* Eigenvalues of a square symmetric matrix A.
!* It calculates by the power method the largest eigenvalue. 
!* Then, it builds a new matrix by removing the last eigenvalue. 
!*
!*          input  : A(:,:) 
!*                   optional lambda_min 
!*                   minimum value of calculated eigenvalues
!*                   if not present, lamba_min is the machine epsilon
!*
!*          output : lambda(k), k  eigenvalue
!*                   U(:,k),  eigenvector associated to lambda(k)
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!****************************************************************
subroutine Eigenvalues_PM(A, lambda, U, lambda_min) 
      real, intent(inout) :: A(:,:)
      real, intent(out) :: lambda(:), U(:,:)
      real, optional,  intent(in) :: lambda_min
  
      integer :: k, N
      real :: l_min 
      logical :: next 
      
      next = .true. 
      N = size(A, dim=1) 
      if ( present(lambda_min) ) then 
                                      l_min = lambda_min  
      else 
                                      l_min = epsilon(1.) 
      end if 
      
      k= 1;   lambda = 0 
      do while (next)  
          
          call Power_method(A, lambda(k), U(:, k) ) 
          A = A - lambda(k) * Tensor_product( U(:, k), U(:, k) ) 
          
          next = lambda(k) > l_min .and. k < N
          k = k+ 1   
      end do 
      
end subroutine 
          
!****************************************************************
! SVD decomposition of a real matrix A 
!*                   A = U S transpose(V)  
!*
!*    input  : 
!*             A(:,:)     input square matrix M x N 
!*             sigma_min  optional value of minimum singular value
!*                        to be obtained.  if not present sigma_min is epsilon  
!*
!*    output : 
!*            sigma(k) singular values
!*                     sigma(k)**2 are  eigenvalues of transpose(A) * A 
!*            V(:,k)   eigenvector of transpose(A) * A (NxN)
!*            U(:,k)   eigenvector of  A * transpose(A) (MxM) 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!****************************************************************
subroutine SVD(A, sigma, U, V, sigma_min) 
      real, intent(in) :: A(:,:)
      real, intent(out) :: sigma(:), U(:,:), V(:,:) 
      real, optional, intent(in) :: sigma_min   
  
      integer :: i, r, N, M ! r< N, rank of AT A 
      real, allocatable :: B(:,:) 
      real :: lambda_min
      
      M = size(A, dim=1); N = size(A, dim=2) 
      if (present(sigma_min)) then 
                           lambda_min = sigma_min**2 
      else 
                           lambda_min = epsilon(1.)  
      end if 
      
      B = matmul( transpose(A), A ) 
    
      call Eigenvalues_PM( B, sigma, V, lambda_min ) 
    
      sigma = sqrt( abs(sigma) ) ! to avoid negative round-off 
      
      r = count( sigma > sqrt(lambda_min) )
      do i=1, r 
           U(:,i) = matmul( A, V(:, i) ) / sigma(i)
      end do  
      
      if (r<M) call Gram_Schmith( r, U )
      if (r<N) call Gram_Schmith( r, V ) 
      
    
end subroutine 

subroutine Gram_Schmith(r, U) 
  integer, intent(in) :: r 
  real, intent(inout) :: U(:,:) 
  
  integer :: k, i, j, N 
  N = size(U, dim=2)
     
  do k=r+1, N 
       
     call random_number( U(:,k) ) 
   
   ! Uk must be orthogonal to U1, U2, ... Uk-1 
     do j=1, k-1
         U(:,k) = U(:,k) - dot_product( U(:,k), U(:,j) ) / norm2(U(:,j)) * U(:,j) 
      end do 
      U(:,k) = U(:,k) / norm2( U(:,k) )
      
  end do 

end subroutine   
    


!****************************************************************
!* Condition number of a real  matrix A 
!*                   Kappa(A)  = norm2(A)  * norm2( inverse(A) ) 
!*
!* The norm2 of a real matrix A is its largest singular value 
!* From the SVD decomposition of a matrix
!*      A = U D transpose(V), 
!*      inverse(A) = V inverse(D) transpose(U) 
!* Hence, norm2( inverse(A) ) is the smallest singular value. 
!* The smallest value is calculated by means of the Inverse power method. 
!*
!*    input  : A(:,:) 
!*
!*    output : 
!*            Condition_number 

!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es)  
!****************************************************************
real function Condition_number(A) 
      real, intent(in) :: A(:,:)
    
  
      integer :: i, j, k, N 
      real, allocatable :: B(:,:), U(:) 
      real :: sigma_max, sigma_min, lambda  
      
      N = size(A, dim=1)
      allocate( U(N), B(N,N) ) 
      B = matmul( transpose(A), A ) 
      
      
     call Power_method( B,  lambda, U ) 
     sigma_max = sqrt(lambda) 
     
     call Inverse_power_method( B,  lambda, U )
     sigma_min = sqrt(lambda) 
     
     Condition_number = sigma_max / sigma_min  
   
    
end function  

!****************************************************************
!* Linear regression 
!*                   y = c1 x + c0  
!*
!* Given a set of points ( y_i, x_i) i=1, .. M 
!* determine the best fit regression by means of SVD decomposiotn. 
!* A x  = b with martrix A of MxN and M>>N (overdetermined)
!* Best fitting: least squeares 
!*      x = V S^-1 U^T b 
!*
!*    input  : x(:), y(:)  
!*
!*    output : 
!*            c(:) 
!*
!*  Author: 2021 Juan A Hernandez (juanantonio.hernandez@upm.es)  
!****************************************************************
function linear_regression(x, y) result(c) 
       real, intent(in) :: x(:), y(:)
       real :: c(2) 
       
   real, allocatable :: A(:, :), U(:,:), V(:,:), sigma(:), b(:) 
   integer :: M, N  
       
     N = 2;  M = size(x)
     allocate( A(M,N), b(M), U(M, M), V(N, N), sigma(N) ) 
     
     A(:,1) = 1 
     A(:,2) = x 
       
     call SVD( A, sigma, U, V ) 
     
     b = matmul( transpose(U), y )/ sigma 
     c = matmul( V, b )
     
end function 








!************************************************************************
!* Band width = 2 * rj  + 1 
!************************************************************************
subroutine LU_band_factorization( A , rj) 
  real, intent(inout) :: A(:, :) 
  integer, intent (in) :: rj 
  
  integer :: N, r
  integer :: k, i, j
  
      
      N =size(A, dim =1) 
      r= rj+1
    
  A(1, :) = A(1,:) 
  A(2:N,1) = A(2:N,1)/A(1,1) 

 
  do k=2, r
      
     do i=k, N
         A(k,i) = A(k,i) - dot_product( A(k, 1:k-1), A(1:k-1, i)  )
     end do 
     
     do j=k+1, N
         A(j,k) = ( A(j,k) - dot_product( A(1:k-1, k), A(j, 1:k-1) ) )/A(k,k)
     end do     
    
  end do
     
  do k=r+1, N-r

   do i=k, k+r     
      A(k,i) = A(k,i) - dot_product(A(k, k-r:k-1), A(k-r:k-1, i) )
   end do 
   
   do j=k+1, k+r-1
      A(j,k) = (A(j,k)-dot_product(A(k-r+1:k-1, k), A(j, k-r+1:k-1) ))/A(k,k)
   end do     
   
  end do 
 
   do k=N-r+1, N
      
    do i=k, N     
        A(k,i) =  A(k,i) - dot_product(A(k, 1:k-1), A(1:k-1, i)  ) 
    end do 
    
    do j=k+1, N
       A(j,k) = (A(j,k) - dot_product(A(1:k-1, k), A(j, 1:k-1) ) )/A(k,k)                     
    end do 
                      
   end do

end subroutine


!************************************************************************
!* Band width = 2 * rj  + 1 
!************************************************************************
function Solve_LU_band( A, b, rj) 
  real, intent(in) :: A(:, :), b(:) 
  real :: Solve_LU_band( size(b) )
  integer, intent (in) :: rj

   real :: y (size(b)), x(size(b))
   integer :: i, N, r
  
   N = size(b) 
   r = rj +1

   y(1) = b(1)
   
   do i=2,r
          y(i) = b(i) - dot_product( A(i, 1:i-1), y(1:i-1) )   
   enddo
   
    do i=r+1,N
          y(i) = b(i) - dot_product( A(i, i-r:i-1), y(i-r:i-1) )
   enddo

   x(N) = y(N) / A(N,N)

   do i=N-1, N-r, -1
     x(i) = (y(i) - dot_product( A(i, i+1:N), x(i+1:N) ) )/ A(i,i)
   end do
    
   do i=N-r-1, 1, -1
     x(i) = (y(i) - dot_product( A(i, i+1:i+r), x(i+1:i+r) ) )/ A(i,i)
   end do

   Solve_LU_band = x 

end function 


!************************************************************************
!* Tensor product A_ij = U_i V_j 
!************************************************************************
function Tensor_product(U, V) result(A) 
     real, intent(in) :: U(:), V(:) 
     real A( size(U), size(V) ) 

     integer ::i, j, N, M 
     
     N = size(U) 
     M = size(V) 
     
     
     do i=1, N; do j=1, M 
         
         A(i,j) = U(i) * V(j) 
         
     end do;  end do 
     
     

end function 


end module 
