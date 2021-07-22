program main
    implicit none
    ! Todo - Declare variables at the end
    real, dimension(2) :: x_guess, x_value
    real :: mod_psi_infinite, psi_x_guess_real, psi_x_guess_img 
    logical :: NR_stop
    ! Calculation of mod_psi_infinite as per root(rho_infinite/m)
    ! rho_infinite = |psi|^2 integration
    ! m is mass of boson = 4
    mod_psi_infinite = (1/4)**0.5

    
    ! initial guess
    x_guess(0) = 0.1
    x_guess(1) = 0.1


    NR_stop = .FALSE.
    do while(NR_stop == .FALSE.)
        jacboian_inverse = jacobian(x_guess)
        psi_x_guess_real = psi(x_guess(0), x_guess(1))(0)
        psi_x_guess_img = psi(x_guess(0), x_guess(1))(1)

        x_value(0) = x_guess(0) - (jacboian_inverse(0,0)*psi_x_guess_real + jacboian_inverse(0,1)*psi_x_guess_img)
        x_value(1) = x_guess(1) - (jacboian_inverse(1,0)*psi_x_guess_real + jacboian_inverse(1,1)*psi_x_guess_img)
        
        NR_stop = should_stop(x_value, mod_psi_infinite)
    enddo

end program main

logical function should_stop(x_guess, mod_psi_infinite)
    ! x_guess is a (x,y) array 
    logical :: found_vortex
    found_vortex = .FALSE.

    delta = 0.000001
    x = x_guess(0)
    y = x_guess(1)

    mod_psi_x_guess = (psi(x,y)(0)**2 + psi(x,y)(1)**2 )**0.5
    if (mod_psi_x_guess - delta*mod_psi_infinite <= 0) then
        found_vortex = .TRUE.
    endif 

    should_stop = found_vortex
    
end function should_stop


! TODO - psi will be obtained from the data, currently assuming a wavfunc at 2.4.2 - https://gfd.whoi.edu/wp-content/uploads/sites/18/2018/03/Grisouard_report_136567.pdf
function psi(x,y) result(psi_values)
    implicit none
    real, dimension(2) :: psi_values
    real :: x_coord, y_coord

    ! x,y co-ordinate = (i-1)/2.0
    ! wave func = x + (0.1x - 0.2y)i
    x_coord =  (i-1)/2.0
    y_coord =  (i-1)/2.0
    
    psi_values(0) = x_coord
    psi_values(1) = y_coord

    ! function will return (psi_real, psi_imaginary) array at a given x,y position
end function psi



function jacobian(x_guess) result(j_inverse)
    implicit none
    
    real :: j_inverse(2,2), x, y
    h = 0.000001 ! a very small number
    x = x_guess(0)
    y = x_guess(1)

    dau_x_psi_r = (psi(x+h,y)(0) - psi(x-h,y)(0))/2h
    dau_x_psi_i = (psi(x+h,y)(1) - psi(x-h,y)(1))/2h

    dau_y_psi_r =  (psi(x,y+h)(0) - psi(x,y-h)(0))/2h   
    dau_y_psi_i =  (psi(x,y+h)(1) - psi(x,y-h)(1))/2h

    det = dau_x_psi_r*dau_y_psi_i - dau_x_psi_i*dau_y_psi_r
    j_inverse(0,0) = dau_y_psi_i/det
    j_inverse(0,1) = -dau_y_psi_r/det
    j_inverse(1,0) = -dau_x_psi_i/det
    j_inverse(1,1) = dau_x_psi_r/det

    jacobian = j_inverse

end function jacobian