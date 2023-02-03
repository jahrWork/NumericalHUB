
!***********************************************************************
! It integrates in time the Initial value boundary problem.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Boundary_Value_Problems
 use Initial_Boundary_Value_Problem1D
 use Initial_Boundary_Value_Problem2D
 use Utilities
 use Temporal_Schemes

 implicit none 
 private
 public :: Initial_Boundary_Value_Problem, & 
           Spatial_discretization,         &
           Spatial_truncation_error,       & 
           Linear_operator

 interface  Initial_Boundary_Value_Problem
    module procedure IBVP1DS, IBVP2DS 
 end interface
 
 interface  Spatial_discretization
    module procedure Spatial_discretization1DS, Spatial_discretization2DS
 end interface
 
 interface  Spatial_truncation_error
    module procedure Spatial_truncation_error1DS, Spatial_truncation_error2DS
 end interface
 
 interface  Linear_operator
    module procedure Linear_operator1DS, Linear_operator2DS  
 end interface
 
end module 