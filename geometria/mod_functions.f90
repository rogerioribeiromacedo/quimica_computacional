module mod_functions
    implicit none

    ! Position of atoms
    type :: atom_structure
        character(len=2) :: symbol
        real :: x, y, z            
    end type atom_structure

contains
    ! Get the structure of atom
    type(atom_structure) function get_structure(file_name) result(list_atoms)
        character(:), allocatable, intent(in) :: file_name
        character(:), allocatable :: linha        
        integer :: i = 0, fim_arquivo, num_linha = 1, fim = 0, inicio = 1
        integer :: value
        type(atom_structure), dimension(:), allocatable :: atoms
        
        ! Opening a file
        open(10, file=file_name, status="old")
        
        do 
            read(10, fmt="(a)", iostat=fim_arquivo) linha
            if (fim_arquivo < 0) then
                exit
            else
                if (num_linha >= 3) then
                    i = i + 1
                    write(*, fmt="(a)") linha

                    ! x
                    fim = scan(linha, " ")
                    write(*, fmt="(a, I4)", advance="No") "Fim: ", fim
                    inicio = fim
                    linha = adjustr(adjustl(trim(linha(inicio: len(linha)))))
                    print *, ""
                    write(*, fmt="(a)", advance="No") "Valor: ", linha
                    !atoms%symbol = linha(1, fim)


                    !atoms%x      = 
                    
                else
                    if (num_linha == 1) then
                        num_linha = num_linha + 1
                        write (linha,'(I3)') value
                        allocate(atoms(value))
                    else
                        num_linha = num_linha + 1
                    end if
                end if
            end if
        end do
        
        ! Closing a file
        close(10)
        
        list_atoms = atom_structure("", 1, 2, 3)

    end function get_structure

   ! Get the filename from a command line
    function get_filename()
        character(len=:), allocatable :: argument
        character(:), allocatable :: get_filename
        integer :: arglen, argc

        get_filename = ""
        argc = command_argument_count()
        if (argc > 0) then
            ! Get the arguments                        
            call get_command_argument(1,length=arglen)
            allocate(character(arglen) :: argument)
            call get_command_argument(1,value=argument)
            argument = trim(argument)

            if (argument(len(argument)-3:len(argument)) == ".xyz") then
                get_filename = adjustr(trim(argument))
            else
                write(*, fmt="(a)") "Rememnber, the file must be in a xyz format, and with a xyz extension."
            end if
        else
            write(*, fmt="(a)") "You must inform the file name (.xyz)"
        end if
    end function get_filename
      
end module mod_functions
