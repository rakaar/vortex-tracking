program main
    implicit none
    real :: x, helo

    x = helo(10.34)/2
    write(*, *) "x is ", x
end program main


real function helo(n) 
    implicit none
    real :: n
    helo = n

end function helo