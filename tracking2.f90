PROGRAM rho_vzu

IMPLICIT NONE

INTEGER(KIND=8) :: rec_len
INTEGER:: Nx,Ny,Nz
INTEGER:: Nxh,Nyh,Nzh
INTEGER:: Nxhp,Nyhp,Nzhp
INTEGER:: Nxpp,Nypp,Nzpp
INTEGER:: nshell
INTEGER:: ios
INTEGER:: i1,i2,i3,ic1,ic2,ic3,cstep
INTEGER:: ix,iy,iz,ir
INTEGER:: factor1,factor2
INTEGER:: ifile
INTEGER:: cloop,cntfile
INTEGER:: cxmin,cxmax,cymin,cymax,czmin,czmax
INTEGER:: ip, order

integer :: index_i, index_j
! Some constants
integer, parameter:: GP = KIND(0.0D0)
real(kind=GP), parameter:: zero  = 0.0_GP
real(kind=GP), parameter:: one   = 1.0_GP
real(kind=GP), parameter:: mone  = -1.0_GP
real(kind=GP), parameter:: pi    = 4.0_GP*atan(1.0_GP)
real(kind=GP), parameter:: half  = 1.0_GP/2.0_GP
real(kind=GP), parameter:: oney4 = 1.0_GP/4.0_GP
real(kind=GP), parameter:: two   = 2.0_GP
real(kind=GP), parameter:: oney3 = 1.0_GP/3.0_GP 
real(kind=GP), parameter:: oney6 = 1.0_GP/6.0_GP
complex(kind=GP), parameter:: czero = CMPLX(0.0,0.0,KIND=GP)!complex(0.0_GP,0.0_GP)
complex(kind=GP), parameter:: zi    = CMPLX(0.0,1.0,KIND=GP)!complex(0.0_GP,1.0_GP)

REAL(KIND=GP):: lengthx,lengthy,lengthz

REAL(KIND=GP):: dx,dy,dz
REAL(KIND=GP):: du1,du2,du3
REAL(KIND=GP):: su
REAL(KIND=GP):: rs
REAL(KIND=GP):: k0,a,b,c
COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: psi
REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: rho
!REAL*4, ALLOCATABLE, DIMENSION(:,:,:):: rho

INTEGER:: r_iopoint
REAL(KIND=GP):: r_time,r_dt

CHARACTER(100)::fnn,prcjj
INTEGER:: filen,nregrpfile,ii,jj,kk,nplnperproc_orignal
INTEGER:: nsplit_orignal,nsplit_current
INTEGER:: i2_loc,i3_loc,myrank
integer:: x,y,z,h, length_of_below_threshold, phase_grad_floored
real :: w_ps, dau_x_psi_r, dau_x_psi_i, dau_y_psi_r, dau_y_psi_i,psi_part, max_psi_square, threshold,find_phase_grad, phase_grad 
real :: integer_check_diff
real:: find_angle_numerator, find_mod_prod
integer, dimension(256,2) :: below_threshold_pts


Nx = 256
Ny = 256

Nxh = Nx/2
Nyh = Ny/2

Nxhp = Nx/2+1
Nyhp = Ny/2+1
Nxpp = Nx+2
Nypp = Ny+2

lengthx = two*pi
lengthy = two*pi

dx = lengthx/Nx
dy = lengthy/Ny


nsplit_orignal = 4
nsplit_current = 1 
myrank = 0
filen = 7

ALLOCATE(psi(1:Nx,1:Ny))
ALLOCATE(rho(1:Nx,1:Ny))



DO ifile=filen,filen

WRITE(fnn,'(i8)') ifile
write(*,*) "fnn is ",fnn

DO i3 = 1,4
jj = (i3-1)*(Nx/4)+1

WRITE(prcjj,'(i8)') i3

write(*,*) 'wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat'
OPEN(UNIT=11,FILE='wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',&
FORM='UNFORMATTED', status='OLD', IOSTAT=ios)

READ(11) r_iopoint,r_time,r_dt
DO i2 = jj,i3*(Nx/4)
!print*,i2
READ(11) (psi(i1,i2),i1=1,Nx)
ENDDO

CLOSE(11)

ENDDO



ENDDO

! find max psi^2
max_psi_square = ABS(psi(1,1))**2
do index_i = 1,Nx
    do index_j = 1,Ny
        if(ABS(psi(index_i, index_j))**2 > max_psi_square) then
            max_psi_square = ABS(psi(index_i, index_j))**2
        endif
    enddo
enddo
write(*,*) "max psi square", max_psi_square

! filter all those within threshold
threshold = 0.15! as per paper it is 0.15
do index_i = 1,Nx
    do index_j = 1,Ny
        if(ABS(psi(index_i, index_j))**2 - max_psi_square*threshold < 0) then
            
            ! check if phase grad is integer or not
            phase_grad = find_phase_grad(index_i, index_j, psi)/2*pi ! divide by pi to check whether its integral multiples of pi or not
            phase_grad_floored = floor(phase_grad)
            integer_check_diff = phase_grad - phase_grad_floored
            if(integer_check_diff < 0.01 .OR. integer_check_diff > 0.99) then
                write(*,*) "[", index_i,",", index_j,"]",","
                ! write(*,*) "psi^2",ABS(psi(index_i, index_j))**2
            
                ! write(*, *) "this might be vortex"
                ! write(*,*) "grad", phase_grad
                ! write(*,*) "------------------------------"
            endif
            
        endif
    enddo
enddo


DEALLOCATE(psi)
DEALLOCATE(rho)

END PROGRAM

real function find_angle_numerator(x1, y1, x2, y2, complex_num)
    integer, parameter:: GP = KIND(0.0D0)
    COMPLEX(KIND=GP), DIMENSION(256,256):: complex_num
    integer :: x1, y1, x2, y2
    real :: a,b

    a = psi_part(x1, y1, complex_num,1) * psi_part(x2, y2, complex_num, 1)  
    
    b =  psi_part(x1,y1,complex_num,2) * psi_part(x2,y2,complex_num,2)

    find_angle_numerator = a + b
    ! find_angle_numerator = 2+3
end function find_angle_numerator

real function find_mod_prod(x1, y1, x2, y2, complex_num)
    integer, parameter:: GP = KIND(0.0D0)
    COMPLEX(KIND=GP), DIMENSION(256,256):: complex_num
    integer :: x1,y1,x2,y2
    real :: v1, v2

    v1 = (psi_part(x1,y1,complex_num,1)**2 + psi_part(x1,y1,complex_num,2)**2)**0.5
    v2 = (psi_part(x2,y2,complex_num,1)**2 + psi_part(x2,y2,complex_num,2)**2)**0.5

    find_mod_prod = v1*v2
end function find_mod_prod

real function find_phase_grad(x_cord, y_cord, complex_num)
    ! TODO: Need to figure out how to calculate phase gradient
    integer, parameter:: GP = KIND(0.0D0)
    COMPLEX(KIND=GP), DIMENSION(256,256):: complex_num
    integer :: x_cord, y_cord, final_x, final_y
    real :: deno1, deno2, deno3, deno4, deno5, deno6, deno7, deno8
    

    if(x_cord >= 256 .or. x_cord <=1  .or. y_cord >= 256 .or. y_cord <= 1) then
       find_phase_grad = 10000000 ! we can't have vortex at the edges, so returning this absurd value
    endif

  
    ! find_phase_grad =  atan(psi_part(x_cord+1,y_cord,complex_num,2)/psi_part(x_cord+1,y_cord,complex_num,1)) + &
    !               atan(psi_part(x_cord+1,y_cord+1,complex_num,2)/psi_part(x_cord+1,y_cord+1,complex_num,1)) + &
    !               atan(psi_part(x_cord,y_cord+1,complex_num,2)/psi_part(x_cord,y_cord+1,complex_num,1)) + &
    !               atan(psi_part(x_cord-1,y_cord+1,complex_num,2)/psi_part(x_cord-1,y_cord+1,complex_num,1)) + &
    !               atan(psi_part(x_cord-1,y_cord,complex_num,2)/psi_part(x_cord-1,y_cord,complex_num,1)) + &
    !               atan(psi_part(x_cord-1,y_cord-1,complex_num,2)/psi_part(x_cord-1,y_cord-1,complex_num,1)) + &
    !               atan(psi_part(x_cord,y_cord-1,complex_num,2)/psi_part(x_cord,y_cord-1,complex_num,1)) + &
    !               atan(psi_part(x_cord+1,y_cord-1,complex_num,2)/psi_part(x_cord+1,y_cord-1,complex_num,1))

    ! num =  psi_part(x_cord+1, y_cord-1, complex_num, 1)*psi_part(x_cord+1, y_cord, complex_num, 1)  + psi_part(x_cord+1, y_cord-1, complex_num, 2)*psi_part(x_cord+1, y_cord, complex_num, 2)

    ! deno = ((psi_part(x_cord+1, y_cord-1, complex_num, 1)**2 + psi_part(x_cord+1, y_cord-1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord+1, y_cord, complex_num, 1)**2 + psi_part(x_cord+1, y_cord, complex_num, 2)**2)**0.5)

    ! find_phase_grad =   acos((psi_part(x_cord+1, y_cord, complex_num, 1)*psi_part(x_cord+1, y_cord+1, complex_num, 1)  + psi_part(x_cord+1, y_cord, complex_num, 2)*psi_part(x_cord+1, y_cord+1, complex_num, 2))/(((psi_part(x_cord+1, y_cord, complex_num, 1)**2 + psi_part(x_cord+1, y_cord, complex_num, 2)**2)**0.5) * ((psi_part(x_cord+1, y_cord+1, complex_num, 1)**2 + psi_part(x_cord+1, y_cord+1, complex_num, 2)**2)**0.5)))) 
                        ! acos((psi_part(x_cord+1, y_cord+1, complex_num, 1)*psi_part(x_cord, y_cord+1, complex_num, 1)  + psi_part(x_cord+1, y_cord+1, complex_num, 2)*psi_part(x_cord, y_cord+1, complex_num, 2))/(((psi_part(x_cord+1, y_cord+1, complex_num, 1)**2 + psi_part(x_cord+1, y_cord+1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord, y_cord+1, complex_num, 1)**2 + psi_part(x_cord, y_cord+1, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord, y_cord+1, complex_num, 1)*psi_part(x_cord-1, y_cord+1, complex_num, 1)  + psi_part(x_cord, y_cord+1, complex_num, 2)*psi_part(x_cord-1, y_cord+1, complex_num, 2))/(((psi_part(x_cord, y_cord+1, complex_num, 1)**2 + psi_part(x_cord, y_cord+1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord-1, y_cord+1, complex_num, 1)**2 + psi_part(x_cord-1, y_cord+1, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord-1, y_cord+1, complex_num, 1)*psi_part(x_cord-1, y_cord, complex_num, 1)  + psi_part(x_cord-1, y_cord+1, complex_num, 2)*psi_part(x_cord-1, y_cord, complex_num, 2))/(((psi_part(x_cord-1, y_cord+1, complex_num, 1)**2 + psi_part(x_cord-1, y_cord+1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord-1, y_cord, complex_num, 1)**2 + psi_part(x_cord-1, y_cord, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord-1, y_cord, complex_num, 1)*psi_part(x_cord-1, y_cord-1, complex_num, 1)  + psi_part(x_cord-1, y_cord, complex_num, 2)*psi_part(x_cord-1, y_cord-1, complex_num, 2))/(((psi_part(x_cord-1, y_cord, complex_num, 1)**2 + psi_part(x_cord-1, y_cord, complex_num, 2)**2)**0.5) * ((psi_part(x_cord-1, y_cord-1, complex_num, 1)**2 + psi_part(x_cord-1, y_cord-1, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord-1, y_cord-1, complex_num, 1)*psi_part(x_cord, y_cord-1, complex_num, 1)  + psi_part(x_cord-1, y_cord-1, complex_num, 2)*psi_part(x_cord, y_cord-1, complex_num, 2))/(((psi_part(x_cord-1, y_cord-1, complex_num, 1)**2 + psi_part(x_cord-1, y_cord-1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord, y_cord-1, complex_num, 1)**2 + psi_part(x_cord, y_cord-1, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord, y_cord-1, complex_num, 1)*psi_part(x_cord+1, y_cord-1, complex_num, 1)  + psi_part(x_cord, y_cord-1, complex_num, 2)*psi_part(x_cord+1, y_cord-1, complex_num, 2))/(((psi_part(x_cord, y_cord-1, complex_num, 1)**2 + psi_part(x_cord, y_cord-1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord+1, y_cord-1, complex_num, 1)**2 + psi_part(x_cord+1, y_cord-1, complex_num, 2)**2)**0.5))) + &
                        ! acos((psi_part(x_cord+1, y_cord-1, complex_num, 1)*psi_part(x_cord+1, y_cord, complex_num, 1)  + psi_part(x_cord+1, y_cord-1, complex_num, 2)*psi_part(x_cord+1, y_cord, complex_num, 2))/(((psi_part(x_cord+1, y_cord-1, complex_num, 1)**2 + psi_part(x_cord+1, y_cord-1, complex_num, 2)**2)**0.5) * ((psi_part(x_cord+1, y_cord, complex_num, 1)**2 + psi_part(x_cord+1, y_cord, complex_num, 2)**2)**0.5)))

    deno1 = find_mod_prod(x_cord+1, y_cord, x_cord+1, y_cord+1, complex_num)
    deno2 = find_mod_prod(x_cord+1, y_cord+1, x_cord, y_cord+1, complex_num)
    deno3 = find_mod_prod(x_cord, y_cord+1, x_cord-1, y_cord+1, complex_num)
    deno4 = find_mod_prod(x_cord-1, y_cord, x_cord-1, y_cord+1, complex_num)
    deno5 = find_mod_prod(x_cord-1, y_cord, x_cord-1, y_cord-1, complex_num)
    deno6 = find_mod_prod(x_cord-1, y_cord-1, x_cord, y_cord+1, complex_num)
    deno7 = find_mod_prod(x_cord, y_cord-1, x_cord+1, y_cord+1, complex_num)
    deno8 = find_mod_prod(x_cord+1, y_cord, x_cord+1, y_cord-1, complex_num)


    find_phase_grad = acos(find_angle_numerator(x_cord+1, y_cord, x_cord+1, y_cord+1, complex_num)/deno1)+&
                      acos(find_angle_numerator(x_cord+1, y_cord+1, x_cord, y_cord+1, complex_num)/deno2)+&
                      acos(find_angle_numerator(x_cord, y_cord+1, x_cord-1, y_cord+1, complex_num)/deno3)+&
                      acos(find_angle_numerator(x_cord-1, y_cord, x_cord-1, y_cord+1, complex_num)/deno4)+&
                      acos(find_angle_numerator(x_cord-1, y_cord, x_cord-1, y_cord-1, complex_num)/deno5)+&
                      acos(find_angle_numerator(x_cord-1, y_cord-1, x_cord, y_cord+1, complex_num)/deno6)+&
                      acos(find_angle_numerator(x_cord, y_cord-1, x_cord+1, y_cord+1, complex_num)/deno7)+&
                      acos(find_angle_numerator(x_cord+1, y_cord, x_cord+1, y_cord-1, complex_num)/deno8)





        ! find_phase_grad = acos()  

end function find_phase_grad

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
