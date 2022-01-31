module show_orbits 
    
    use dislin 
    implicit none 
   
    
contains 
subroutine plot_orbits(N, Nv, U1, U2) 
    integer, intent(in) :: N, Nv  
    real, intent(in) :: U1(0:N, 1:Nv), U2(0:N, 1:Nv) 
    
    real :: xf, yf 
    real :: xmax, ymax 
    
    xmax = maxval( abs(U2(:,1) ) );   ymax = maxval( abs(U2(:,2) ) )
    call metafl("xwin")
    CALL PAGE (4000, 4000)
    call scrmod("reverse")
    call disini  
    call graf(-xmax, xmax, -xmax, xmax/10 , -ymax, ymax, -ymax, ymax/10)    
  
    call color("blue");  
    call curve( U1(:,1), U1(:,2), N+1) 
    
    call color("red");  
    call curve( U2(:,1), U2(:,2), N+1)
    
    call color("blue");  
    call marker(21); call incmrk(-1)
    xf = U1(N,1); yf = U1(N,2) 
    call curve([xf, xf], [yf, yf], 2 )
   
    
    call color("red");  
    call marker(21);  call incmrk(-1) 
    xf = U2(N,1); yf = U2(N,2) 
    call curve([xf, xf], [yf, yf], 2 )
   
    
    call disfin
    
end subroutine  


subroutine plotm(U, mu, N) 
    real, intent(in) :: U(:,:), mu 
    integer, intent(in) :: N 

    
    call metafl("xwin")
    CALL PAGE (4000, 4000)
    call scrmod("reverse")
    call disini  
    call graf(-1.5, 1.5, -1.5, 0.5, -1.5, 1.5, -1.5, 0.5) 
    call curve( U(:,1), U(:,2), N+1) 
  
    call color("blue"); call marker(21); call incmrk(-1)
    call curve([-mu, -mu], [0., 0.], 2 )
    call color("red");  call marker(21);  call incmrk(-1) 
    call curve([1-mu, 1-mu], [0., 0.], 2 )
    
    call disfin
    
end subroutine   
    
end module     