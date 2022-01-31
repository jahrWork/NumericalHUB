!**********************************************************
!* Book:  How to learn Applied maths (amazom.es)
!*        researchgate: Juan A Hernandez Ramos (UPM, Spain)
!**********************************************************    
program main_NumericalHUB

       
       use API_Example_Systems_of_Equations
       use API_Example_Lagrange_Interpolation
       use API_Example_Cauchy_Problem
       use API_Example_Finite_Differences
       use API_Example_Boundary_Value_Problem
       use API_Example_Initial_Boundary_Value_Problem
       use API_Example_IBVP_and_BVP
       
       use my_examples
       use MUSE_orbits
       use API_Example_Fourier_series
       use API_Example_IBVP_Fourier
      
       use API_examples_dependencies
       
       use API_Example_Error_IVBP
          
       implicit none 
       integer :: option = 1  

        
       
do while (option>0) 
    
     write(*,*) "Welcome to NumericalHUB" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Systems of equations  "
     write(*,*) " 2. Lagrange interpolation  " 
     write(*,*) " 3. Finite difference   "
     write(*,*) " 4. ODE Cauchy problems   "
     write(*,*) " 5. Boundary value problems  "
     write(*,*) " 6. Initial-boundary value problems  "
     write(*,*) " 7. Mixed problems: IBVP+BVP  "
     write(*,*) " 8. High order  ODE schemes "
     write(*,*) " 9. Advanced methods and problems "
     
     read(*,*) option 
     
     select case(option)
     case(1) 
         call Systems_of_Equations_examples
         
     case(2) 
         call Lagrange_Interpolation_examples 
         
     case(3) 
         call Finite_difference_examples
       
     case(4)
         call Cauchy_problem_examples
      
     case(5) 
         call BVP_examples
      
     case(6) 
         call IBVP_examples 
      
     case(7) 
         call Nonlinear_Plate_Vibrations
       
     case(8) 
         call Advanced_Cauchy_problem_examples
         
     case(9) 
         call Advanced_problems 
         
         case default
              
     end select 
     
end do
  

contains 
    
subroutine Advanced_problems 

integer :: option = 1  
           
do while (option>0) 
    
     write(*,*) "Welcome to Advanced methods" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Orbits and numerical methods (Master MUSE)  "
     write(*,*) " 2. SVD applications  "
     write(*,*) " 3. Fourier problems  "
     write(*,*) " 4. Chebyshev problems (not yet implemented)  "
     write(*,*) " 5. Heat equation in solids with different conductivities"
     
     read(*,*) option 
     
     select case(option)
     case(1) 
           call Orbits_and_Numerical_Methods
         
     case(2) 
          call Test_SVD
          call Vandermonde_SVD_condition_number
          call Test_Linear_regression
         
     case(3)  
           call Fourier_problems 
         
     case(4)
     case(5)   
           call N_solids_with_different_conductivity
      
        
         case default
              
     end select 
     
end do
    
end subroutine 
    
 
end program  

