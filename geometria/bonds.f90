!-----------------------------------------------------------------------------------
!bonds: A program for calculation of bond lenght and numbers of bonds in
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
    implicit none
    character(len=5), parameter :: PROG_NAME = 'bonds'
    character(len=256) :: file_name = ""
    
    ! Capture the file name
    file_name = get_filename()
    if (len(file_name) > 0 ) then
        call write_file()
    end if

contains
    subroutine write_file()
        character(len=256) :: linha
        integer :: fim_arquivo
        ! Opening a file
        open(10, file=file_name, status="old")
        
        write(*, "(a)") "", "Structure: "
        do 
            read(10, fmt="(a)", iostat=fim_arquivo) linha
            if (fim_arquivo < 0) then
                exit
            else
                write(*, "(a)") trim(linha)
            end if
        end do
        
        ! Closing a file
        close(10)
    end subroutine write_file


    character(len=256) function get_filename() result(file)
        character(len=:), allocatable :: argument
        integer :: arglen, argc

        file = ""
        argc = command_argument_count()
        if (argc > 0) then
            ! Get the arguments                        
            call get_command_argument(1,length=arglen)
            allocate(character(arglen) :: argument)
            call get_command_argument(1,value=argument)
            argument = trim(argument)

            if (argument(len(argument)-3:len(argument)) == ".xyz") then
                file = adjustr(trim(argument))
            else
                write(*, fmt="(a)") "Rememnber, the file must be in a xyz format, and with a xyz extension."
            end if
        else
            write(*, fmt="(a)") "You must inform the file name (.xyz)"
        end if
            
    end function get_filename
end program bonds