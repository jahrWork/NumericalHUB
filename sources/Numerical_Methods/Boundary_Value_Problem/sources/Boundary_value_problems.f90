module Boundary_value_problems 
 use Boundary_value_problems1D
 use Boundary_value_problems2D
 use Boundary_value_problems3D

 implicit none 
 private
 public :: Boundary_Value_Problem ! It solves a boundary value problem

 interface Boundary_Value_Problem
      module procedure Boundary_Value_Problem1D,        &
                       Boundary_Value_Problem1D_system, &
                       Boundary_Value_Problem2D,        & 
                       Boundary_Value_Problem2D_system, & 
                       Boundary_Value_Problem3D_system
 end interface
end module
