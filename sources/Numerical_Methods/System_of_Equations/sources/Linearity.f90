module Linearity
   
implicit none
    
    abstract interface  

       real function F1Ds(x, u, ux, uxx) 
            real, intent(in) :: x, u, ux, uxx                 
       end function 
       
        function F1D_system(x, u, ux, uxx) 
            real, intent(in) :: x, u(:), ux(:), uxx(:)   
            real :: F1D_system(size(u)) 
       end function 
       
       
       real function F2Ds(x, y, u, ux, uy, uxx, uyy, uxy)
            real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
       end function
       
        function F2D_system(x, y, u, ux, uy, uxx, uyy, uxy)
            real, intent(in) :: x, y, u(:), ux(:), uy(:), &
                                uxx(:), uyy(:), uxy(:)
            real :: F2D_system(size(u)) 
        end function
                   
    end interface
  

    contains

logical function Function_Linearity( F, Nv, N ) result(L)
       integer, intent(in) :: Nv, N      
       interface 
            function F(U, Nv) 
              real, intent(in) :: U(:)
              integer, intent(in) :: Nv 
              real :: F(Nv) 
            end function
       end interface 
       
      
    integer :: i
    real :: U(N), V(N), e1(Nv), e2(Nv)
    real, parameter :: eps = 1e-6
    logical :: F_linear(Nv) 
    
    L = .true. 
    do i=1, N 
            U = 0 
            V = 0
            call random_number(U(i))
            call random_number(V(i))
        !** First requirement e1 = 0
            e1 = F(U + V, Nv) -  F(U, Nv) - F(V, Nv) 
    
       !** Second requirement e2 = 0
           e2 = F(3 * U, Nv) - 3 * F(U, Nv)
           

           F_linear  = abs(e1) < eps .AND. abs(e2) < eps
           L = L .and. all(F_linear == .true. ) 
    end do 
    

end function  
    
    
logical function Linearity_BVP1D( F ) result(L)
    procedure (F1Ds) :: F
 
       
    L = Function_Linearity( Fv, 1, 3)
    write(*,*) "The actual scalar BVP 1D is linear  = ", L 
   
contains

    function Fv( u, Nv )
        real, intent(in) :: u(:)
        integer, intent(in) :: Nv 
        real :: Fv(Nv) 
        
        Fv = F( 0.0, u(1), u(2), u(3) )
        
    end function
    
end function


logical function Linearity_BVP1D_system( F, Nv ) result(L)
    procedure (F1D_system) :: F
    integer, intent(in) :: Nv 
       
    L = Function_Linearity( Fv, Nv, 3*Nv)
    write(*,*) "The actual vector BVP 1D is Linear = ", L 
   
contains

    function Fv( u, Nv )
        real, intent(in) :: u(:)
        integer, intent(in) :: Nv 
        real :: Fv(Nv) 
        
        Fv = F( 0.0, u(1:Nv), u(Nv+1:2*Nv), u(2*Nv+1:3*Nv) )
              
    end function
    
end function



logical function Linearity_BVP2D( F ) result(L)
    procedure (F2Ds) :: F
 
       
    L = Function_Linearity( Fv, 1, 6)
    write(*,*) "The actual scalar BVP 2D is linear  = ", L 
   
contains

    function Fv( u, Nv )
        real, intent(in) :: u(:)
        integer, intent(in) :: Nv 
        real :: Fv(Nv) 
        
        Fv = F( 0.0, 0.0, u(1), u(2), u(3), u(4), u(5), u(6) )
        
    end function
    
end function



logical function Linearity_BVP2D_system( F, Nv ) result(L)
    procedure (F2D_system) :: F
    integer, intent(in) :: Nv 
       
    L = Function_Linearity( Fv, Nv, 6*Nv)
    write(*,*) "The actual vector BVP 2D is Linear = ", L 
   
contains

    function Fv( u, Nv )
        real, intent(in) :: u(:)
        integer, intent(in) :: Nv 
        real :: Fv(Nv) 
        
        Fv = F( 0.0, 0.0, u(1:Nv), u(Nv+1:2*Nv), u(2*Nv+1:3*Nv),   & 
                u(3*Nv+1:4*Nv), u(4*Nv+1:5*Nv), u(5*Nv+1:6*Nv)  ) 
        
    end function
    
end function




end module
    