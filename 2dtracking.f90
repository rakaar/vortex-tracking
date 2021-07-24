program main
    implicit none
    ! Todo - Declare variables at the end
    real, dimension(2) :: x_guess, x_value
    real :: mod_psi_infinite, psi_x_guess_real, psi_x_guess_img, psi_real, psi_img, jacobian(2,2), jacboian_inverse(2,2)
    logical :: NR_stop, should_stop
    ! Calculation of mod_psi_infinite as per root(rho_infinite/m)
    ! rho_infinite = |psi|^2 integration
    ! m is mass of boson = 4
    mod_psi_infinite = (1/4)**0.5

    
    ! initial guess
    x_guess(0) = 0.1
    x_guess(1) = 0.1


    NR_stop = .FALSE.
    do while(NR_stop .eqv. .FALSE.)
        jacboian_inverse = jacobian(x_guess(0), x_guess(1))
        psi_x_guess_real = psi_real(x_guess(0), x_guess(1))
        psi_x_guess_img = psi_img(x_guess(0), x_guess(1))

        x_value(0) = x_guess(0) - (jacboian_inverse(0,0)*psi_x_guess_real + jacboian_inverse(0,1)*psi_x_guess_img)
        x_value(1) = x_guess(1) - (jacboian_inverse(1,0)*psi_x_guess_real + jacboian_inverse(1,1)*psi_x_guess_img)
        
        NR_stop = should_stop(x_value, mod_psi_infinite)
        write(*,*) "loop done"
    enddo

end program main

logical function should_stop(x_guess, mod_psi_infinite)
    ! x_guess is a (x,y) array 
    logical :: found_vortex
    real, dimension(2) :: x_guess
    real :: mod_psi_infinite, delta,x,y, mod_psi_x_guess, psi_real, psi_img

    found_vortex = .FALSE.

    delta = 0.000001
    x = x_guess(0)
    y = x_guess(1)

    mod_psi_x_guess = (psi_real(x,y)**2 + psi_img(x,y)**2 )**0.5
    if (mod_psi_x_guess - delta*mod_psi_infinite <= 0) then
        found_vortex = .TRUE.
    endif 

    should_stop = found_vortex
    
end function should_stop


! TODO - psi will be obtained from the data, 
! currently assuming a wavfunc at 2.4.2 - https://gfd.whoi.edu/wp-content/uploads/sites/18/2018/03/Grisouard_report_136567.pdf
real function psi_real(x,y)
    implicit none
    real :: x,y
    write(*,*) "x is",x
    write(*,*) "psi_real is",(x-1)/2.0

    psi_real = (x-1)/2.0 !psi_real = x

end function psi_real

real function psi_img(x,y)
    implicit none
    real :: x,y, a1, a2

    a1 = 0.000000221
    a2 = -0.0000000433

    write(*,*) "x is",x
    write(*,*) "psi_real is",a1*((x-1)/2.0) + a2*((y-1)/2.0)

    psi_img = a1*((x-1)/2.0) + a2*((y-1)/2.0)

end function psi_img




function jacobian(x,y) result(j_inverse)
    implicit none
    
    real :: j_inverse(2,2), x, y, h, psi_real, psi_img, dau_x_psi_r, dau_x_psi_i, dau_y_psi_r, dau_y_psi_i, det
    real, dimension(2) :: x_guess

    h = 0.000001 ! a very small number
    

    dau_x_psi_r = (psi_real(x+h,y) - psi_real(x-h,y))/(2*h)
    dau_x_psi_i = (psi_img(x+h,y) - psi_img(x-h,y))/(2*h)

    dau_y_psi_r =  (psi_real(x,y+h) - psi_real(x,y-h))/(2*h)   
    dau_y_psi_i =  (psi_img(x,y+h) - psi_img(x,y-h))/(2*h)

    det = dau_x_psi_r*dau_y_psi_i - dau_x_psi_i*dau_y_psi_r
    j_inverse(0,0) = dau_y_psi_i/det
    j_inverse(0,1) = -dau_y_psi_r/det
    j_inverse(1,0) = -dau_x_psi_i/det
    j_inverse(1,1) = dau_x_psi_r/det

    ! jacobian = j_inverse

end function jacobian