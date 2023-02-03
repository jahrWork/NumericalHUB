module Boundary_value_problems
    
 use Boundary_value_problems1D
 use Boundary_value_problems2D

 implicit none 
 private
 public :: Boundary_Value_Problem ! It solves a boundary value problem

 interface Boundary_Value_Problem
     module procedure  Boundary_Value_Problem1DS, & 
                       Boundary_Value_Problem2DS 
      
 end interface
end module
