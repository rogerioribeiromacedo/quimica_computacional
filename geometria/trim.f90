program trim
    implicit none
    character(:), allocatable :: teste

    teste = "ABC"
    write(*, "(a)") teste
    write(*, "(i3)") len(teste)
end program trim