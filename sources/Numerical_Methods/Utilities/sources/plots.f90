
module plots

use Utilities
implicit none

! Overloading. First argument can be a vector or a matrix 
interface plot_parametrics 
      module procedure plot_parametrics1, plot_parametrics2 
end interface

interface plot
      module procedure plot1, plotn 
end interface

private 
public :: plot, plot_parametrics, plot_contour
public :: open_unit_file, ini_dislin, float_to_str
public :: plot_legends, plot_title, dislin_ini 


contains
    
!************************************************************************************** 
!  It plots y-x graph for different parametrics
!    
! INPUT: 
!         x(:,:), 
!         y(:,:) : first index represents points
!                  second index represents the index of the parametric  
!         legends: array of character strings with the legends of the parametrics 
!         x_label : x-axis label 
!         y_label : y-axis label 
!         title   : title of the plot      
!         path    : optional file path to write a data or plt file of different parametrics   
!                   when this argument is present, this subroutine creates: 
!                    1) ascii data file with the information of parametrics 
!                    2) latex file to paint those parametrics from latex  
!         graph_type: TODO   
!
! Author: Javier Escoto, Juan A. Hernandez (2019)  
!**************************************************************************************
subroutine plot_parametrics2(x, y, legends, x_label, y_label, title, path, graph_type) 
   real, intent(in) :: x(:, :), y(:,:) 
   character(len=*), intent(in) :: legends(:), x_label, y_label
   character(len=*), optional, intent(in) :: title, path, graph_type  
   
   !  character(len=3000) :: cbuf 
     real :: xmin, xmax, ymin, ymax 
     integer :: M, Np, j   
     character(len=20) :: gtype
   
     
    Np = size(legends)  
    M = size(x, dim=1)  
    xmin = minval( x ); xmax = maxval( x ); 
    ymin = minval( y ); ymax = maxval( y ); 
    
    if (present(graph_type)) then
        gtype = graph_type
    else 
        gtype = "blackwhite"    
    end if 
    
    
    if (present(path)) then 
       call Graph_Tex_and_plt( x, y, x_label, y_label, legends, title, path, gtype )
    end if 
    
    
    call scrmod("reverse")
    call page(3000, 3000) 
    call metafl("xwin")
    call disini  
    call texmod("on") 
    
    call labels("FEXP", "y")
    call name(x_label, 'X-axis' )
    call name(y_label, 'Y-axis' )
    call graf(xmin, xmax, xmin, (xmax-xmin)/5, ymin, ymax, ymin, (ymax-ymin)/5 )
    
    do j=1, Np  
           call chncrv('color')
           call curve( x(:, j), y(:,j), M ) 
    end do 
    
    call plot_legends(legends)
    
    if (present(title)) call plot_title([title,""]) 
    
      
    call disfin 


end subroutine     
    
    

!********************************************************************************** 
!  Same than before subroutine but first argument is a vector instead of a matrix 
!**********************************************************************************
subroutine plot_parametrics1(x, y, legends, x_label, y_label, title, path, graph_type) 
   real, intent(in) :: x(:), y(:,:) 
   character(len=*), intent(in) :: legends(:), x_label, y_label
   character(len=*), optional, intent(in) :: title, path, graph_type 
   
     integer :: N, M, j   
     real, allocatable :: xp(:,:) 
     
     N = size(y, dim=1)
     M = size(y, dim=2) 
     allocate ( xp(N, M) ) 
     
     do j=1, M 
         xp(:,j) = x 
     end do 
   
     call plot_parametrics2( xp, y, legends, x_label, y_label, title, path, graph_type) 

end subroutine 

!***********************************************************************************
! It creates a plt file 
!***********************************************************************************
subroutine Graph_Tex_and_plt ( x, y, x_label, y_label, legends, title, path, graph_type )
      real, intent(in) :: x(:,:), y(:,:) 
      character(len=*), intent(in) :: x_label, y_label, legends(:), title, path 
      character(len=*), optional, intent(in) :: graph_type
      
      integer :: i, j, N, M  
      character(len=50) :: substrings(10), file_name 
      character(len=200) :: folder
      real, allocatable :: matrix(:,:)
     
      
      call split( trim(path), "/", substrings, N  ) 
      
      folder = substrings(1) 
      do i=2, N-1
        folder = trim(folder) //"/"// trim( substrings(i) )  
      end do
      file_name = substrings(N) 
      write(*,*) "folder =", trim(folder)
      write(*,*) "file_name =", trim(file_name) 
    
      N = size(x, dim=1) 
      M = size(x, dim=2) 
      allocate( matrix(N, 2*M) ) 
      
      do j=0, M-1 
          matrix(:, 2*j+1) = x(:, j+1) 
          matrix(:, 2*j+2) = y(:, j+1) 
      end do     
     
      ! *** file_name.plt creation
      call write_graph_plt( matrix, trim(folder), trim(file_name), x_label, legends )
      
      ! *** file_name.tex creation
      call tex_graph(N, folder, file_name, x_label, y_label, legends, title, graph_type ) 
      

end subroutine
    

!*********************************************************************
! It writes a matrix of data. Rows are points and columns variables
!*********************************************************************
subroutine write_graph_plt( data_array, folder, plt_name, x_label, legends ) 
    real, intent(in) :: data_array(0:,:)
    character(len=*), intent(in) :: folder, plt_name, x_label, legends(:) 
    
    integer:: N, M, i, j, funit
    character(len = 300 ) :: Header
    character(len=30) :: cx, temp
    
    N = size(data_array, dim=1) -1 
    M = size(data_array, dim=2)
    
    ! *** Header of the plt file writing
    Header = ' VARIABLES = '
    
    do j = 1, size(legends)
        write(temp, '(i1)') j
        if (j==1) then 
            cx = trim(x_label)
        else 
            cx = trim(x_label)//trim(temp)
        end if 
        Header = trim(Header) //"  " // trim(cx) //"  "// trim(legends(j))  
    end do
    
    ! *** Writing of the plt in the given folder
    funit = open_unit_file( trim(folder)//"/", trim(plt_name)//".plt" )
  !  open(20, file = trim(folder) // "\"  // trim(plt_name) // ".plt" )
         write(funit, '(500A)') trim(Header)
         write(funit, *) 'ZONE I=', N+1,',DATAPACKING=POINT'
         
         do  i = 0 , N
             write (funit, '(500e23.15)') data_array(i,1:M)
         end do
    
     close(funit)
     
end subroutine
!*********************************************************************
! It writes file_name.tex to process with latex 
!********************************************************************* 
 subroutine tex_graph( N, folder, file_name, x_label, y_label, legends, title, graph_type  )
     integer, intent(in) :: N 
     character(len=*), intent(in) :: folder, file_name, x_label, y_label, legends(:), title, graph_type 
         
     character(len=200) :: line(20) 
     character(len=200) :: clegends ! " $ \sin x $, $ \cos x $ " 
     character(len=200) :: ccolumns ! "0/1,2/3"
     character(len=200) :: tex_file,  data_file 
     character(len=10)  :: cmarks 
     integer ::i, M, funit 
     
     data_file = trim(folder) // "/" // trim(file_name) // ".plt" 
       
     M = size(legends) 
     clegends = trim(legends(1)) 
     ccolumns = "0/1"
     do i=2, M 
         clegends = trim(clegends) // ", " // trim(legends(i))  
         ccolumns = trim(ccolumns) // "," // number_to_str(2*i-2) // "/" // number_to_str(2*i-1)
     end do 
     cmarks = integer_to_str( min(N/25+1, 25) ) 
     
     
     line(1) = "  \begin{tikzpicture}[]                                                                               " 
     line(2) = "     \begin{axis}[ width = \textwidth,                                                                "
     line(3) = "        title style={at={(0.5,-0.45)},anchor=south}, title={" // trim(title) // "},                   "
     line(4) = "        xlabel={ "// trim(x_label) // " },                                                            "
     line(5) = "        ylabel={ "// trim(y_label) // " },                                                            "   
     if (graph_type=="color") then 
     line(6) = "        {}, legend style={font=\small},  {no marks, legend columns = 1 } ]                            "
     else if(graph_type=="blackwhite")   then 
     line(6) = "        {}, legend style={font=\small},  {mark repeat ="//trim(cmarks)//", cycle list name = black white, legend columns = 1 }, legend pos = north east  ]"
     else 
           write(*,*) " Error tex_graph"; stop 
     end if  
     
     
     line(7) = "        \foreach \x/\y in { " // trim(ccolumns) // "}                                                 "
     line(8) = "        \addplot table [skip first n=2, x index=\x, y index=\y]{" // trim(data_file) // "   };        "
     line(9) = "        \legend{ " // trim(clegends) // " }                                                           "
     line(10)= "     \end{axis}                                                                                       "
     line(11)= "  \end{tikzpicture}                                                                                   "
 
     tex_file = trim(folder) // "\" // trim(file_name) // ".tex"
     
     funit = open_unit_file( trim(folder)//"/", trim(file_name)//".tex" )
   !  open(20, file = tex_file)
     do i=1, 11 
         write(funit,'(A)') trim(line(i))  
     end do 
     
     close(funit) 

 
 end subroutine 
 
    
 
 character(len=1) function number_to_str(n) result(c) 
      integer, intent(in) :: n 

    write(c,'(i1)')  n
    
 end function 
 
 
    
!************************************************************************************** 
!  It plots contour graphs (isolines and color maps)  
!    
! INPUT: 
!         x(:)      : x cordinates  
!         y(:)      : y coordinates 
!         z(:,:)    : matrix to plot z= z(x,y) 
!         levels(:) : iso line levels to show 
!         x_label   : x-axis label 
!         y_label   : y-axis label 
!         legend    : title of the plot  
!         path      : optional file path to write a plt file. Format three columns: x, y, z   
!                     when this argument is present, this subroutine creates: 
!                    1) ascii data file with the information of contour map 
!                    2) latex file to paint the counter map from latex  
!         graph_type: "isolines" or "color"   
!
! Author: Javier Escoto, Juan A. Hernandez (2019, 2021)  
!**************************************************************************************
subroutine plot_contour(x, y, z, x_label, y_label, levels, legend, path, graph_type) 
   real, intent(in) :: x(:), y(:), z(:,:)
   real, optional, intent(in) :: levels(:)
   character(len=*), intent(in) ::  x_label, y_label
   character(len=*), optional, intent(in) :: legend, path, graph_type  
   
    
     real :: xmin, xmax, ymin, ymax, zmin, zmax 
     integer ::  i, Nl 
     character(len=20) :: g_type
     real, allocatable :: iso_levels(:) 

    xmin = minval( x ); xmax = maxval( x ); 
    ymin = minval( y ); ymax = maxval( y ); 
    zmin =0; zmax = 0; 
    
    if (present(levels)) then 
        Nl = size(levels) 
        allocate(iso_levels(Nl)) 
        iso_levels = levels 
        zmin = minval(levels); zmax = maxval(levels); 
    else 
        Nl = 10
        allocate(iso_levels(Nl)) 
    endif 
    
     if(zmax==zmin) then 
        zmin = minval(z); zmax = maxval(z); 
        iso_levels = [ (zmin + (zmax-zmin)*(i-1)/real(Nl-1), i=1, Nl) ] 
     end if
     
    call scrmod("reverse")
    call page(3000, 3000) 
    call metafl("xwin")
    call disini  
    call texmod("on") 
    
    call name(x_label, 'X-axis' )
    call name(y_label, 'Y-axis' )
    call labdig(-2, "y");
    call labdig(-2, "x");
    
    call color("white")
   
    if (present(legend)) call titlin(legend, 1)
        
    if (.not.present(graph_type)) then 
              g_type = "color" 
    else 
              g_type = graph_type
    end if 
    
    if (g_type == "color") then 
        
      call graf3(xmin, xmax, xmin, (xmax-xmin)/5,  & 
                 ymin, ymax, ymin, (ymax-ymin)/5,  &
                 zmin, zmax, zmin, (zmax-zmin)/5 )
       call shdmod('cell','CONTUR')
       call conshd(x, size(x), y, size(y), z, iso_levels,  Nl )
       
    else if (g_type == "isolines") then 
       
     call graf(xmin, xmax, xmin, (xmax-xmin)/5,  & 
               ymin, ymax, ymin, (ymax-ymin)/5 )  
      do i=1, Nl 
       if (mod(i, Nl/10) == 0 ) then 
           call labels("FLOAT", "CONTUR")
       else                   
           call labels("NONE", "CONTUR")
       end if    
       call contur( x, size(x), y, size(y), z, iso_levels(i) ) 
      end do 
      
    end if   
    
    if (present(legend)) call title 
    call disfin 

    if (present(path)) then 
      call contour_Tex_and_plt( x, y, z, iso_levels, x_label, y_label, legend, path, g_type )
    end if
    
    deallocate(iso_levels) 

end subroutine   



!***********************************************************************************
! It creates a plt file 
!***********************************************************************************
subroutine contour_Tex_and_plt ( x, y, z, levels, x_label, y_label, legend, path, graph_type )
      real, intent(in) :: x(:), y(:), z(:,:), levels(:)   
      character(len=*), intent(in) :: x_label, y_label, legend, path, graph_type 
      
      integer :: i, j, k, N, M  
      character(len=100) :: substrings(10), file_name 
      character(len=200) :: folder
      real, allocatable :: matrix(:,:)
      
      call split( trim(path), "/", substrings, N  ) 
      
      folder = substrings(1) 
      do i=2, N-1
        folder = trim(folder) //"/"// trim( substrings(i) )  
      end do
      file_name = substrings(N) 
      write(*,*) "folder =", trim(folder)
      write(*,*) "file_name =", trim(file_name) 
    
      
      N = size(x);  M = size(y) 
      allocate( matrix(N*M, 3) ) 
      
      k = 0 
      do i=1, N 
          do  j=1, M 
              k = k + 1 
              matrix(k, :) = [ x(i), y(j), z(i,j) ] 
          end do 
      end do 
     
      ! *** file_name.plt creation
      call write_graph_plt( matrix, trim(folder), trim(file_name), x_label, [legend] )
      
      ! *** file_name.tex creation
      call tex_contour(folder, file_name, x, x_label, y, y_label, levels, legend, graph_type ) 

end subroutine


!*********************************************************************
! It writes file_name.tex to process with latex
!********************************************************************* 
 subroutine tex_contour( folder, file_name, x, x_label, y, y_label, levels, legend, graph_type  )
     real, intent(in) :: x(:), y(:), levels(:) 
     character(len=*), intent(in) :: folder, file_name, x_label, y_label, legend, graph_type 
     
     character(len=800) :: line(20) 
     character(len=400) :: clevels ! " 0.2, 0.7, 0.8  " 
     character(len=200) :: tex_file,  data_file 
     character(len=10) :: cN, cM 
     character(len=7) :: cxmax, cxmin, cymax, cymin
     integer :: i, Nl, funit 
    
     
     data_file = trim(folder) // "/" // trim(file_name) // ".plt" 
       
     
     Nl = size(levels)
     clevels = float_to_str( levels(1) ) 
     do i=2, Nl
         clevels = trim(clevels) // "," // float_to_str( levels(i) ) 
     end do 
     cN = integer_to_str( size(x) ) 
     cM = integer_to_str( size(y) ) 
     cxmax = float_to_str( maxval(x) ); cxmin = float_to_str( minval(x) );
     cymax = float_to_str( maxval(y) ); cymin = float_to_str( minval(y) );
     
     line(1) = "  \begin{tikzpicture}[]  " 
     if (graph_type=="color") then 
     line(2) = "     \begin{axis}[ height = 0.8\textwidth, width = 0.8\textwidth,                                     " 
     line(3) = "        title style={at={(0.5,-0.5)}, yshift=-0.1}, title={" // trim(legend) // "},                   "
     else 
     line(2) = "     \begin{axis}[ height = \textwidth, width = \textwidth,                                           "
     line(3) = "        title style={at={(0.5,-0.4)}, yshift=-0.05}, title={" // trim(legend) // "},                   "
     end if 
    
     line(4) = "        xlabel={ "// trim(x_label) // " },    xmin= "//trim(cxmin)//", xmax = "//trim(cxmax)//",         "
     line(5) = "        xlabel={ "// trim(y_label) // " },    ymin= "//trim(cymin)//", ymax = "//trim(cymax)//",         "
     
     if (graph_type=="color") then 
         line(6) = "        view={0}{90}, colormap/bluered, colorbar ]                                  "
     else 
         line(6) = "         view={0}{90}, colormap/blackwhite ]                                    "
     end if 
     
     if (graph_type=="color") then  
     line(7) = "  \addplot3[smooth, surf,                                                                             "  
     line(8) = "            shader=interp, z buffer=sort,                                                              "
     
     else 
     line(7) = "  \addplot3[contour/labels=false, contour gnuplot = { levels={ " // trim(clevels) // " }},             "  
     line(8) = "            mesh/ordering = y varies,                                                                "
    
     end if 
     line(9) = "            mesh/rows =  "// trim(cN) // ",                                                           " 
     line(10) ="            mesh/cols =  "// trim(cM) // "                                                            " 
     line(11) ="           ]                                                                                          "  
     line(12) ="           table [skip first n=2] {" // trim(data_file) // "};                                        "
     line(13) ="     \end{axis}                                                                                       "
     line(14) =" \end{tikzpicture}                                                                                    "  
        
 
     tex_file = trim(folder) // "\" // trim(file_name) // ".tex"
     
     funit = open_unit_file( trim(folder)//"/", trim(file_name)//".tex" )
     !open(funit, file = tex_file)
     do i=1, 14 
         write(funit,'(A)') trim(line(i))  
     end do 
     
     close(funit) 

 
 end subroutine 
 
 
 

character(len=7) function float_to_str(x) result(c) 
      real, intent(in) :: x

    write(c,'(f7.2)')  x
  
    
 end function



character(len=3) function integer_to_str(N) result(c) 
      integer, intent(in) :: N 

    write(c,'(i3)')  N
    
end function

 
integer function open_unit_file( folder, file ) result(unit) 
 character(len=*), intent(in) :: folder, file 
  
      integer :: funit = 20 
      character(len=200) :: bfolder 

 do while (1) 
     
    open(unit = funit, file = trim(folder)//trim(file), err = 30) 
    unit = funit 
    exit  
    
30  call replace_char( "/", "\", folder, bfolder) 
    call execute_command_line("mkdir" // adjustl(trim(bfolder) ) )  
 end do    
 
end function 

subroutine replace_char( f, r, string, new_string ) 
        character(len=1), intent(in) :: f, r 
        character(len=*), intent(in) :: string 
        character(len=*), intent(out) :: new_string 
  
      integer :: i,N
      character(len=1) c
      
      N = len(trim(string))
      do i=1, N 
          c = string(i:i) 
          if (c==f) then 
                   new_string(i:i) = r
          else 
                   new_string(i:i) = c
          end if 
      end do 
  
end subroutine  
 



subroutine ini_dislin( xmin, xmax,  ymin, ymax ) 
    real, intent(in) ::  xmin, xmax, ymin, ymax 

    
    call metafl("xwin")
    CALL PAGE (4000, 4000)
    call scrmod("reverse")
    call disini  
    call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
   
end subroutine

!****************************************************************
!  plot y versus x 
!  If legend contains "space" plots y(i,:) parmetrics are time 
!  If legend contains "time"  plots y(:,i) parmetrics are space 
!****************************************************************
subroutine  plotn(x, y, legend) 
     real, intent(in) :: x(0:), y(:,0:) 
     character(len=*), optional, intent(in) :: legend 

    real :: xmin, xmax, ymin, ymax
    integer :: Nx, Np, i, j  
    Nx = size(y,dim=2)-1; Np = size(y,dim=1);
     
    xmin = minval(x); xmax = maxval(x) 
    ymax = maxval(y); ymin = minval(y) 
    call scrmod("reverse") 
    call metafl("xwin")
   
    call page(4000, 4000)
    call disini 
    call winfnt("Courier New Bold")
    call texmod("on") 
    call labels("fexp", "x")
    call labels("fexp", "y")
    !call labdig(-2, "y");
    !call labdig(-2, "x");
    call height(46)
    call hname(60)
    call graf(xmin, xmax, xmin, (xmax-xmin)/5,  & 
              ymin, ymax, ymin, (ymax-ymin)/5 )
   
    if (present(legend)) then 
    !  call height(60)
      call height(46)
      call titlin(legend, 1)
      call title
    end if 
    
    
    do i=1, Np 
       call incmrk(-1);  call marker(21);call color("red"); 
       call curve( x, y(i,:), Nx+1)
       call incmrk(0);  ;call color("white"); 
       call curve( x, y(i,:), Nx+1) 
    end do 
    
   ! call height(60)
    call height(46)
    call title
    call disfin 

end subroutine 

!****************************************************************
!  plot y versus x 
!****************************************************************
subroutine  plot1(x, y, legend) 
     real, intent(in) :: x(0:), y(0:)
     character(len=*), intent(in) :: legend 

    real :: xmin, xmax, ymin, ymax
    integer :: i, N 
         
    N = size(x)-1
    xmin = minval(x); xmax = maxval(x) 
    ymax = maxval(y); ymin = minval(y) 
    call scrmod("reverse") 
    call metafl("xwin")
   
    call page(4000, 4000)
    call disini 
    call winfnt("Courier New Bold")
    call texmod("on") 
    call labels("fexp", "x")
    call labels("fexp", "y")
    !call labdig(-2, "y");
    !call labdig(-2, "x");
    call height(46)
    call hname(60)
    if ( index(legend,"time")>0 )   call name("$t$", "x" )
    if ( index(legend,"space")>0 )  call name("$x$", "x" )
    call graf(xmin, xmax, xmin, (xmax-xmin)/5,  & 
              ymin, ymax, ymin, (ymax-ymin)/5 )
   
    call titlin(legend, 1)
    call title
    call curve(x, y, size(x))
    call color("red");  call marker(21);  call incmrk(-1) 
    call curve(x, y, size(x))
       
    call disfin 

end subroutine 


subroutine dislin_ini(xmax, xmin, ymax, ymin) 
 real, intent(in) :: xmax, xmin, ymax, ymin 
 
    call scrmod("reverse") 
    call metafl("xwin")
   
    call page(4000, 4000)
    call disini 
    call graf(xmin, xmax, xmin, (xmax-xmin)/5, ymin, ymax, ymin, (ymax-ymin)/5 )
    call height(40);
    call winfnt("Courier New Bold")
    call color("white") 

end subroutine 
    
    
subroutine plot_legends(legends) 
 character(len=*), intent(in) :: legends(:)

    character(len=3000) :: cbuf 
    integer :: j, Np 
 
   
    call winfnt("Courier New Bold")
    Np = size(legends)
    call height(40);
    call color("white") 
   
    call legini(cbuf, Np, 30)
    call legclr 
    call legtit("")
 
    do j=1, Np
          call leglin(cbuf, legends(j), j)
    end do 
    
    call legpos(10, 10)
    call legend(cbuf, Np)
    
end subroutine     


subroutine plot_title(tit) 
 character(len=*), intent(in) :: tit(:)

    integer :: i, N 
    
    call color("white") 
    N  = size(tit) 
    call height(30);
    call linesp(3.) 
    
    do i=1, N 
      call titlin( tit(i), i)
    end do 
    
    call title
    
   

end subroutine


end module