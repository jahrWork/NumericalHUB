

module Utilities 

 implicit none 
 private
 public ::                  &
           Tensor_product,  &     ! A_ij = u_i v_j 
           my_reshape,      &     ! Creates a pointer
           split,           &     ! split a char string based on delimiter 
           Binary_search,   &     ! Binary search in an ordered list 
           Linear_operatorF,& 
           Tensor_product_R3
  
 public :: operator  (.x.)

 interface operator(.x.) 
       module procedure Tensor_product 
 end interface

 
contains



!***********************************************************************
! Given a vector function L: RN-> RN . 
! If the function is linear, it gives the system matrix A
! If L = A U + b, it gives the matrix A 
!***********************************************************************
function Linear_operatorF( L_operator, N) result(A)
      interface
         function L_operator(U)
            real, intent(in) :: U(:) 
            real :: L_operator( size(U) ) 
         end function 
      end interface 
      integer, intent(in) :: N  
      real :: A(N, N) 
      
     integer :: i   
     real, target ::  U(N), b(N), L(N)  
   
     U=0 
     b = L_operator(U) 
     
     do i=1, N
              U = 0 
              U(i) = 1 
              L = L_operator(U) - b 
              A(:, i) = L 
     end do 
             
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

!************************************************************************
!* Tensor product A_ij = U_i V_j 
!************************************************************************
function Tensor_product_R3(U, V, W) result(A) 
     real, intent(in) :: U(:), V(:), W(:) 
     real A( size(U), size(V), size(W) ) 

     integer :: i, j, k, N, M, L 
     
     N = size(U) 
     M = size(V) 
     L = size(W)
     
     
     do i=1, N; do j=1, M; do k=1, L  
         
         A(i,j,k) = U(i) * V(j) * W(k)  
         
     end do;  end do; end do  
     
end function     

!************************************************************************
!* It creates a pointer pU 
!************************************************************************
subroutine my_reshape( U, N1, N2,  pU )
     integer, intent(in) :: N1, N2
     real, target, intent(in) :: U(1:N1, 1:N2)  
     real, pointer, intent(out) :: pU(:,:) 
          
     pU => U  
     
end subroutine 

 
!*****************************************************************************
!* It splits a string into N different substrings based on the delimiter
!*****************************************************************************   
subroutine split( string, delimiter, substrings, N  ) 
   character(len=*), intent(in) :: string 
   character(len=1), intent(in) :: delimiter
   character(len=*), intent(out) :: substrings(:) 
   integer, intent(out) :: N 


      integer :: i, index, k 
   
    i = 1   
    k = 1 
    do while (i<=len(string) ) 
        
        index = scan( string(i:), delimiter) 
        
        if (index>0) then 
                          substrings(k) = string(i:i+index-2)
                          i = i + index 
                          k = k + 1 
        else 
             substrings(k) = string(i:len(string))
             N = k 
             exit 
        endif 
        
    end do 
    
   
end subroutine 


    
!*****************************************************************************
!* It converts a integer  to an integer  
!*****************************************************************************   
character(len=2)  function  integer_to_char( i ) result(c) 
                      integer, intent(in) :: i
 
        if (i<10) then 
                      write(c, '(i1)') i
        else 
                      write(c, '(i2)') i
        end if 
        
      
end function    
 


!*****************************************************************************
!* It searches the element x in an ordered list V ( x = V(ix) ) 
!* If the element x is not in the list it returns -1 
!*****************************************************************************
subroutine Binary_search(x, V, ix)
    real, intent (in) :: x, V(0:)
    integer, intent (out) :: ix

  integer:: m, upper_limit, lower_limit, N
  N = size(V) - 1 
  
  lower_limit = 0; upper_limit = N;
 
  do while( upper_limit > lower_limit ) 
     
     m = ( upper_limit + lower_limit ) / 2
     
     if ( x > V(m) ) then 
                          lower_limit = m + 1 
     else 
                          upper_limit = m 
     endif 
     
  end do 
 
  !if (x == V(m) ) then 
  !      ix = m 
  !else 
  !      ix = -1 
  !endif 
  ix = lower_limit 
 
 
end subroutine
   
    
   

!*****************************************************************************
!* It calculates trace of matrix A 
!*****************************************************************************   
real function trace(A) 
  real, intent(in) :: A(:,:) 

  integer :: i
  real :: S 
  
  S = 0 
  do i=1, size(A, dim=1) 
      
      S = S + A(i,i) 
      
  end do 
  
  trace = S 
  

end function   

!*****************************************************************************
!* It allocates a Vandermonde matrix of dimension M 
!*****************************************************************************
function Vandermonde(M) result(A) 
  integer, intent(in) :: M 
  real, allocatable :: A(:,:) 

  integer :: i, j 
 
  allocate ( A(M, M) ) 
  
  do i=1, M; do j=1, M 
      A(i,j) = ( i / real(M) ) **(j-1)
  end do; end do     
  

end function


 
!****************************************************************************************************
!* It looks for the given name in an database (array of names) and if returns its position 
!****************************************************************************************************
logical function  equal(a, b)  
   class(*), intent(in) :: a, b 

   
   select type (a) 
       
       type is (integer) 
            select type (b) 
            type is (integer)
            equal = (a == b)  
            end select
            
        
       type is (character(len=*)) 
            select type (b) 
            type is (character(len=*))
            equal = (trim(a) == trim(b))  
            end select
     
       type is (real) 
            select type (b) 
            type is (real)
            equal = (a == b)  
            end select
            
       class default
             write(*,*) "ERROR " 
             stop 
       
   end select
   
   

   
end function 

  
 


!***********************************************************
! it reads a filename and save its data in matrix A 
!***********************************************************
    function load_matrix( filename ) result(A) 
        character(len=*), intent(in) :: filename 
        real, allocatable :: A(:,:) 
        
        integer :: i_file = 15       ! file unit to read 
        character(len=400) :: header ! text header row 
        integer ::  N, M             ! rows and columns 
        integer :: i 
               
        open( unit = i_file, file = filename) 
        
        read(unit = i_file, fmt = '(A)') header 
       ! write(*,'(2A)' ) "header =", trim(header)  
        
        M = columns( header, "," )
        N = 0
        do while(1) 
            read(unit = i_file, fmt ='(A)', end = 10) header
            N = N + 1
         !   write(*,*) " N = ", N 
         !   write(*,*) trim(header) 
        end do 
10      write(*,*) " # of rows =", N, "# of columns = ", M  
                
        allocate( A(N,M) )
        
        rewind(i_file)
        read(i_file, *) ! header 
        do i=1, N 
            read( i_file, *) A(i, :) 
        !    write(*,'(100f10.3)') A(i,:) 
        end do 
        
        close(i_file)
        
    end function 
    
!************************************************************    
! it couns the number of colummns by counting delimiters c
!*************************************************************    
    integer function columns( string, c) result(M) 
       character(len=*), intent(in) :: string 
       character, intent(in) :: c 
       
       integer :: i, N 
       
       N = len_trim(string)
       
       M = 1
       do i=1, N 
           if (string(i:i) == c ) then 
               M = M + 1 
           end if 
       end do 
       
    end function 


!***********************************************************
! it reads a filename and save its data in matrix A 
!***********************************************************
subroutine save_matrix( filename, header, cformat, A ) 
        character(len=*), intent(in) :: filename, header, cformat   
        real, intent(in) :: A(:,:) 
        
        integer :: i_file = 15       ! file unit to write 
        integer ::  N, M             ! rows and columns 
        integer :: i 
               
        open( unit = i_file, file = filename) 
        N = size( A, dim=1) 
        M = size( A, dim=2) 
        
        write(unit = i_file, fmt = '(A)') header 
      
        do i=1, N 
            if (cformat == "integer") then 
                  write( i_file, fmt = '(100i12)' ) int(A(i, :))
            else 
                  write( i_file, fmt = '(100e15.7)' ) A(i, :)
            end if 
            
        end do 
       
        close(i_file)
        
end subroutine   
    
    

end module 
    
    
    

    