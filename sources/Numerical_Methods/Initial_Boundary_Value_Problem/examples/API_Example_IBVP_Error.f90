module API_Example_Error_IVBP

  
    use Initial_Boundary_Value_Problem1D
    use Temporal_Schemes
    use Collocation_methods
    use Temporal_error
    use plots
    
implicit none

    contains
    

subroutine Error_Heat_equation_1D

       integer, parameter :: Nx = 20, Nt = 400 
       real :: x(0:Nx), Time(0:Nt), U(0:Nt,0:Nx), U0(0:NX)
       real :: Error_x(0:Nt, 0:Nx), Error_t(0:Nt, 0:Nx), R(0:Nx)
       integer, parameter :: Ne = 5
       real :: log_Et(Ne), log_Nt(Ne), log_Ex(Ne/2), log_Nx(Ne/2)
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 1, dt  
       integer :: Order = 2 
       integer :: i 
       real :: ordert
    
       
     write (*, '(A50)') 'Error 1D Solution (Heat equation)'
     
     dt = (tf-t0)/Nt
     Time = [ (t0 + dt*i, i=0, Nt ) ] 
     x(0) = x0; x(Nx) = xf
       
     call Temporal_Error_IBVP1D(       & 
         
          Heat_equation, Heat_BC, Heat_IC, Euler, time, x, Order, U, Error_t)
     
     call Spatial_Error_IBVP1D(  & 
         
          Heat_equation, Heat_BC, Heat_IC, time, x, Order, U, Error_x )   
     
     call Truncation_Spatial_Error_IBVP1D(  & 
         
          Heat_equation, Heat_BC, Heat_IC, x, Order, R )
     
    call Grid_Initialization( "nonuniform", "x",  x, Order)
    U0 = Heat_IC(x) 
    call Temporal_convergence_rate( time, F_Cauchy, U0, Euler, ordert, log_Et, log_Nt ) 
    
    call plot(x, U, "Temperature along space for different tn")   
    call plot(time, Error_t, "Error along time with Euler scheme for different positions xj")
    call plot(log_Nt, log_Et, "Euler:Convergence rate for different time steps") 
    
    call plot(x, R,"Truncation space error")
    call plot(x, Error_x, "Error along space (temporal scheme:RK4) for different  tn")  
    
    
    
    
contains 
 
function F_Cauchy( U, t ) result(F) 
                      real ::  U(0:), t, F(0:size(U)-1) 
           
  F = Spatial_discretization1D( & 
      Heat_equation, Heat_BC, x, U, t )
         
      
end function

real function Heat_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            Heat_equation =   uxx
                 
end function 

function Heat_BC(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u, ux 
  real :: BC 
  
        if (x==x0) then
                            BC = u - exp( -50*(t-0.5)**2 ) 

        else if (x==xf) then
                            BC = u

        else
             write(*,*)  "Error in Heat_BC"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function

function Heat_IC(x) result(U) 
  real, intent(in) :: x(:) 
  real :: U(size(x)) 
    
   U = exp( -10*x**2 ) 

end function

end subroutine 


end module
    