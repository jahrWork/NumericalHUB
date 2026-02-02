     call set_solver(family_name="wABM")
     do j=1, Np  
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,         &
                            "$x$", "$y$", "(a)", path(1) )
    
     call set_solver(family_name="ABM")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,        &
                               "$x$", "$y$", "(b)", path(2) )
     
end subroutine
