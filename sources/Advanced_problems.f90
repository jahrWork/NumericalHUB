module Advanced_problems
    
     
       use API_Example_Chebyshev_interpolation
     
       use MUSE_orbits
       use API_Example_Fourier_series
       use API_Example_IBVP_Fourier
       use API_Example_IBVP_Chebyshev
      
       use API_examples_dependencies
     
       
       use Navier_Stokes_cavities
       use Boundary_layer
       use Special_IBVP
       use Burgers 
          
    
       implicit none 
       private
       public :: Spectral_and_advanced_problems
    
contains 
        
    
        
subroutine Spectral_and_advanced_problems 

integer :: option 
       
option = 1
do while (option>0) 
    
     write(*,*) "Welcome to Advanced methods" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Fourier problems  "
     write(*,*) " 2. Chebyshev problems  "
     write(*,*) " 3. Heat equation in solids with different conductivities"
     write(*,*) " 4. Navier Stokes problems"
     write(*,*) " 5. Shock capturing in Burgers equation with upwinding methods"
     
     read(*,*) option 
     
     select case(option)
         
     case(1)  
           call Fourier_problems 
         
     case(2)
          call Chebyshev_problems
          
     case(3)   
           call N_solids_with_different_conductivity
      
     case(4)       
         call Cavities_examples
         call Compressible_boundary_layer_examples
        
     case(5)       
         call Burgers_upwinding
         call Burgers_Cauchy
         
     case default
              
     end select 
     
end do
    
end subroutine 



    
end module 
    