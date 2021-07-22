program main
    implicit none
    real :: x

    x = helo(10.34)
    write(*, *) "x is ", x
end program main


real function helo(n) 
    implicit none
    real :: n
    helo = n

end function helo