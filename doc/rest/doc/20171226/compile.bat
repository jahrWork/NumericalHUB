
echo off

del %1.exe 


set PATH=..\..\fortran\mingw64\bin;%PATH%

set COPTIONS= -I./lib -J./lib -Wall -freal-4-real-8 -ffree-line-length-500 -fstack-arrays -fbounds-check
set LOPTIONS= -I./release -I./lib -J./lib 
set LIBRARIES= -luser32 -lgdi32 -lopengl32 --static
 

rem echo %COPTIONS%

gfortran -c %COPTIONS% ./sources/dislin_d.f90 -o ./release/dislin_d.o

if %3.==. (
    
    if %2.==. (
    
        echo ---------------one file to compile----------------------
        gfortran -c  %COPTIONS%  ./sources/%1.f90 -o ./release/%1.o
        gfortran     %LOPTIONS%  ./release/%1.o  ./lib/dismg64_d.a ./lib/Numerical_methods.a  -o  ./release/%1.exe %LIBRARIES% 
        
    ) else  ( 
    
        echo ---------------two files to compile-----------------
        gfortran -c  %COPTIONS%   ./sources/%2.f90 -o ./release/%2.o 
        gfortran -c  %COPTIONS%   ./sources/%1.f90 -o ./release/%1.o
        gfortran     %LOPTIONS%   ./release/%1.o  ./release/%2.o ./lib/dismg64_d.a ./lib/Numerical_methods.a -o  ./release/%1.exe %LIBRARIES%
    )

) else (

    echo three files to compile 
    gfortran -c %COPTIONS%   ./sources/%3.f90 -o ./release/%3.o
    gfortran -c %COPTIONS%   ./sources/%2.f90 -o ./release/%2.o
    gfortran -c %COPTIONS%   ./sources/%1.f90 -o ./release/%1.o
    gfortran    %LOPTIONS%   ./release/%1.o  ./release/%2.o ./release/%3.o ./lib/dismg64_d.a  ./lib/Numerical_methods.a -o  ./release/%1.exe %LIBRARIES%
  
)

pause 



.\release\%1.exe


