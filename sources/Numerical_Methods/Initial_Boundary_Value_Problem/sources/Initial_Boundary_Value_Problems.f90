
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
           Spatial_discretization

 interface  Initial_Boundary_Value_Problem
    module procedure IBVP1D, IBVP1D_system, IBVP2D, IBVP2D_system
 end interface
 
 interface  Spatial_discretization
   module procedure Spatial_discretization1D, Spatial_discretization1D_system, & 
                    Spatial_discretization2D, Spatial_discretization2D_system
 end interface

end module 