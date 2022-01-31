module Dependencies
   
implicit none
	
	 interface  
     
       function F_RNXD_RD(X) result(F) 
            real, intent(in) :: X(:, :)
            real :: F( size(X, dim=1) ) 
       end function 
       
     end interface
	
   interface  
  
       real function F_RN_R(u)
            real, intent(in) :: u(:)
       end function
       
       real function F_IBVP1D(x, t, u, ux, uxx) 
            real, intent(in) :: x, t, u, ux, uxx                 
       end function 
    
    end interface
    
contains
   
    
 
function is_Jacobian_zero( Nv, Nd, F )  result(D) 
    integer, intent(in) :: Nv, Nd 
    procedure (F_RNXD_RD) :: F 
    logical :: D(Nv, Nd) 
        
    integer :: i, j
    real ::  U0(Nv, Nd), e(Nv, Nd), dF(Nv), h=0.04  
     
    do j=1, Nd 
        do i=1, Nv 
         call random_number(U0)
          e = 0  
          e(i, j) = 1 
          dF = F(U0 + e * h) - F(U0) 
          D(i, j) = any( dF /= 0 ) 
        end do 
    end do 

end function    
    

    
    
function Function_Dependencies( F, N, P_ux ) result(D)

    procedure (F_RN_R) :: F
    integer, intent(in) :: N, P_ux 
    logical :: D(1:N-P_ux+1) 
    
    integer :: i, j
    real :: h( N ), dF, u( N ), dx, r(2)
    logical :: Daux( 1:N, 1:2)
    
    
    call random_number(r)
    dx = 0.04
    
    do i = P_ux, n
        
        do j = 1, 2
            u = r(j)
            h = 0
            h(i) = 1
            dF=F(u + h * dx) - F(u)
                if(dF /= 0) then
                    Daux(i,j) = .true.
                else
                    Daux(i,j) = .false.
                end if
        end do

        if(Daux(i,1) .eqv. Daux(i,2)) then
            D(i-P_ux+1) = Daux(i,1)
        else 
            D(i-P_ux+1) = .true.
        end if
    end do

end function



  
function IBVP1D_Dependencies( F ) result(d)
    procedure (F_IBVP1D) :: F
    logical :: D( 2 )
    
    integer :: N = 5
    integer :: P_ux = 4
    
    D = Function_Dependencies( F_v, N, P_ux)
    write(*,*) " IBVP1D_Dependencies :", D 
        
contains

    real function F_v( u )
    
        real, intent(in) :: u(:)
        
        F_v = F(u(1), u(2), u(3), u(4), u(5))
        
    end function
end function


end module