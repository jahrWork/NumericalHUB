    
module API_examples_dependencies
   
      use Dependencies
      use Dependencies_IBVP2D 
      use Dependencies_BC 
     implicit none
	


contains 
 
subroutine Test_Dependencies2D_system

 integer, parameter :: Nv = 2 
 logical :: D(Nv, 5)
 logical :: D_BC(Nv, 2), D_BC1(2)
 
 
   D = IBVP2D_Dependencies_system( Nv, Wave_equation2D )
   write(*,*) " L(v,w)  = [w,  vxx +vyy+ vx ] "
   write(*,*) "Variable =", 1, " vx, vy, vxx, vyy, vxy ", D(1, :) 
   write(*,*) "Variable =", 2, " wx, wy, wxx, wyy, wxy ", D(2, :)
   
   
   D_BC = BC_IBVP2D_Dependencies_system( Nv, Wave_BC2D, 0., 0., 0., 0. )
   write(*,*) " BC(v,w)  = [v,  w ] "
   write(*,*) "Variable =", 1, " vx, vy ", D_BC(1, :) 
   write(*,*) "Variable =", 2, " wx, wy ", D_BC(2, :) 
   
   
   D_BC = BC_IBVP2D_Dependencies_system( Nv, Wave_BC2D2, 0., 0., 0., 0. )
   write(*,*) " BC(v,w)  = [vx+vy, wy] "
   write(*,*) "Variable =", 1, " vx, vy ", D_BC(1, :) 
   write(*,*) "Variable =", 2, " wx, wy ", D_BC(2, :) 
   !
   D_BC1 = BC_IBVP2D_Dependencies( Heat_BC2D, 0., 0., 0., 0. )
   write(*,*) " BC(v,w)  = ux "
   write(*,*) "Variable =", 1, " ux, uy ", D_BC1(:) 
   
   

end subroutine 

function Wave_equation2D( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
      real, intent(in) ::  x,y,t,u(:),ux(:),uy(:),uxx(:),uyy(:),uxy(:)
      real :: L(size(u))
            
            real :: v, vx, vxx, vyy, w, wx  
            
            v = u(1);  vx = ux(1); vxx = uxx(1); vyy = uyy(1)  
            w = u(2);  wx = ux(2);
        
            L(1)  = w 
            L(2)  = vxx +vyy + vx 
                 
end function 

function Wave_BC2D(x, y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
  
        real :: v, w 
        v = u(1) 
        w = u(2)
        
        if (x==0 .or. x==0 .or. y==0.or. y==0 ) then
                            BC = [v, w] 
        else
             write(*,*)  "Error in BC2D_waves"; stop 
        endif
      
end function

function Wave_BC2D2(x, y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
  
        real :: v, vx, vy, w, wx ,wy  
        v = u(1); vx = ux(1) ; vy = uy(1)  
        w = u(2); wx = ux(2) ; wy = uy(2)  
        
        if (x==0 .or. x==0 .or. y==0.or. y==0 ) then
                            BC = [vx+vy, wy] 
        else
             write(*,*)  "Error in BC2D_waves"; stop 
        endif
      
end function

function Heat_BC2D(x, y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u, ux, uy 
  real :: BC
  
        
        
        if (x==0 .or. x==0 .or. y==0.or. y==0 ) then
                            BC =ux
        else
             write(*,*)  "Error in BC2D_waves"; stop 
        endif
      
end function 
  


end module 
    
    
    
      