PROGRAM rho_vzu

IMPLICIT NONE

integer:: Nx,Ny
integer:: ios
integer:: i1,i2,i3
integer:: ifile

! 2d tracking 
integer:: x_guess, y_guess, h, x, y
logical :: NR_stop, should_stop
real :: dau_x_psi_r, dau_x_psi_i, dau_y_psi_r, dau_y_psi_i, psi_part,det, psi_x_guess_real, psi_x_guess_img
real :: condensate_density, mass_of_boson, mod_psi_inf, delta, w_ps
real, dimension(2,2) :: j_inverse

! Some constants
integer, parameter:: GP = KIND(0.0D0)

! Reading psi
complex(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: psi
integer:: r_iopoint
real(KIND=GP):: r_time,r_dt
character(100)::fnn,prcjj
integer:: filen,jj

Nx = 256
Ny = 256

filen = 1

allocate(psi(1:Nx,1:Ny))


do ifile=1,filen
    write(fnn,'(i8)') ifile
    write(*,*) "fnn is ",fnn

    do i3 = 1,4
        jj = (i3-1)*(Nx/4)+1
        write(prcjj,'(i8)') i3

        write(*,*) 'wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat'
        open(UNIT=11,FILE='wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
        FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

        READ(11) r_iopoint,r_time,r_dt
        do i2 = jj,i3*(Nx/4)
            READ(11) (psi(i1,i2),i1=1,Nx)
        enddo

        close(11)
    enddo
enddo

!write(*,*) psi

write(*,*) "#################################################################################"
write(*,*) "testing --------------------------------------------------------->",psi(1,2)
write(*,*) (ABS(psi(1,2)))**2
write(*,*) real(psi(1,2))
write(*,*) aimag(psi(1,2))
write(*,*) "testing --------------------------------------------------------->",psi(256,256)
write(*,*) "testing --------------------------------------------------------->",psi(100,120)

write(*,*) "#################################################################################"


! choose a random point
x_guess = 200
y_guess =  200

h = 1

NR_stop = .false.

! for stopping condition, TODO: Need to know proper values
condensate_density = 1
mass_of_boson = 4
mod_psi_inf = (condensate_density/mass_of_boson)**0.5
delta = 0.1

x = x_guess
y = y_guess
do while(NR_stop .eqv. .FALSE.)
    write(*,*) "***start of loop-x,y", x,y
    ! 1 for real, 2 for imaginary
    dau_x_psi_r = (psi_part(x+h,y,psi,1) - psi_part(x-h,y,psi,1))/(2*h)
    ! write(*,*) "dau_x_psi_r", dau_x_psi_r
    dau_x_psi_i = (psi_part(x+h,y,psi,2) - psi_part(x-h,y,psi,2))/(2*h)
    ! write(*,*) "dau_x_psi_i", dau_x_psi_i
    
    dau_y_psi_r =  (psi_part(x,y+h,psi,1) - psi_part(x,y-h,psi,1))/(2*h)  
    ! write(*,*) "dau_y_psi_r", dau_y_psi_r
    
    dau_y_psi_i =  (psi_part(x,y+h,psi,2) - psi_part(x,y-h,psi,2))/(2*h)
    ! write(*,*) "dau_y_psi_i", dau_y_psi_i
    
    ! write(*,*) "dausssss",dau_x_psi_r, dau_x_psi_i, dau_y_psi_r, dau_y_psi_i

    det = dau_x_psi_r*dau_y_psi_i - dau_x_psi_i*dau_y_psi_r
    ! write(*,*) "det1",dau_x_psi_r*dau_y_psi_i
    ! write(*,*) "det2",dau_x_psi_i*dau_y_psi_r
    ! write(*,*) "det",det

    j_inverse(1,1) = dau_y_psi_i/det
    j_inverse(1,2) = -dau_y_psi_r/det
    j_inverse(2,1) = -dau_x_psi_i/det
    j_inverse(2,2) = dau_x_psi_r/det

    ! write(*,*) "Jacobian Inverse in the loop is this ", j_inverse


    psi_x_guess_real = psi_part(x,y,psi,1)
    psi_x_guess_img = psi_part(x,y,psi,2)
    x = x - (j_inverse(1,1)*psi_x_guess_real + j_inverse(1,2)*psi_x_guess_img)
    y = y - (j_inverse(2,1)*psi_x_guess_real + j_inverse(2,2)*psi_x_guess_img)
    write(*,*) "***end of loop-x,y", x,y

    ! suppose at end of loop it jumps away
    if(x > 256) then
        x = 256
    else if (x < 1) then
        x = 1
    endif

    if(y > 256) then
        y = 256
    else if (y < 1) then
        y = 1
    endif

   
    write(*,*) ">>>>mod_psi",ABS(psi(x,y)**2)
    call sleep(1)
    ! write(*,*) ">>>>delta*mod_psi_inf",delta*mod_psi_inf
    
    ! write(*,*) ">>>>>>>>>>>>>>>>>>>>>>>mod_psi - delta x mod_psi_inf",ABS(psi(x,y)) - delta*mod_psi_inf
    
    
    ! write(*,*) "w_ps", w_ps
    ! w_ps = h_bar * grad_psi_real * grad_psi_img (take h_bar = 1)
    w_ps = dau_x_psi_r*dau_y_psi_i - dau_y_psi_r*dau_x_psi_i 
    ! condition to stop: mod_psi < delta x mod_psi_inf
    ! if((ABS(psi(x,y)) - delta*mod_psi_inf <= 0) .and. w_ps /= 0) then
    !     NR_stop = .true.   
    ! endif

    ! a simpler condition for now, till proper values of condensate, mass are known
     if(ABS(psi(x,y))**2 <= 0.1) then
        NR_stop = .true.   
    endif
enddo



deallocate(psi)

END PROGRAM

real function psi_part(x_cord, y_cord, complex_num, real_or_img)
    integer, parameter:: GP = KIND(0.0D0)
    COMPLEX(KIND=GP), DIMENSION(256,256):: complex_num
    integer :: x_cord, y_cord, final_x, final_y, real_or_img
    
    final_x = x_cord
    final_y = y_cord
    
    if(x_cord > 256) then
        final_x = 256
    else if (x_cord < 1) then
        final_x = 1
    endif

    if(y_cord > 256) then
        final_y = 256
    else if (y_cord < 1) then
        final_y = 1
    endif

    ! write(*,*) "final_x, final_y",final_x, final_y
    ! write(*,*) real_or_img
    if(real_or_img .eq. 1) then
        psi_part = real(complex_num(final_x, final_y))
    else
        psi_part = aimag(complex_num(final_x, final_y))
    endif
    
end function psi_part
