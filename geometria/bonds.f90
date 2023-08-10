!-----------------------------------------------------------------------------------
!bonds: A program for calculation of bond lenght and numbers of bonds.
!-----------------------------------------------------------------------------------
!>@file   Main program
!>@author Rogério Ribeiro Macêdo
!>        Laboratory of Computational Chemistry - LaQC
!>        Federal University of Itajuba - UNIFEI
!>        <https://en.unifei.edu.br/research-centers-labs/computational-chemistry-laboratory/>
!>@email  rogerioribeiromacedo@gmail.com
!-----------------------------------------------------------------------------------
!   Copyright 2023-2023 Rogério Ribeiro Macêdo
!
!   This program is free software: you can redistribute it and/or modify 
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------
!Written by Rogério Ribeiro Macêdo
!Laboratory of Computational Chemistry - LaQC
!e-mail: rogerioribeiromacedo@gmail.com
!-----------------------------------------------------------------------------------
program bonds
    use mod_functions
    implicit none
    character(len=5), parameter :: PROG_NAME = 'bonds'
    character(:), allocatable :: file_name
    type(atom_structure), allocatable, dimension(:) :: atoms
    
    ! Capture the file name
    file_name = get_filename()
    if (len(file_name) > 0 ) then

        atoms = get_structure(file_name)
    end if

contains
    subroutine write_file()
        
    end subroutine write_file



end program bonds