program main
    implicit none
    
    integer, parameter:: GP = KIND(0.0D0)
    INTEGER:: r_iopoint, ios,Nx, Ny, i3, jj, i1, i2
    REAL(KIND=GP):: r_time,r_dt
    CHARACTER(100)::fnn,prcjj

    COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:):: psi

    ALLOCATE(psi(1:Nx,1:Ny))

    Nx = 256
    Ny = 256
    ! fnn = 1
    ! prcjj = 1 ! temp, later loop

    ! WRITE(prcjj,'(i8)') 1
    ! WRITE(fnn,'(i8)') 1

    ! OPEN(UNIT=11,FILE='wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',FORM='UNFORMATTED', status='OLD', IOSTAT=ios)
    
    ! READ(11) r_iopoint,r_time,r_dt
    ! write(*,*) r_iopoint,r_time,r_dt

    ! actual
    WRITE(fnn,'(i8)') 1
    DO i3 = 1,4
        WRITE(prcjj,'(i8)') i3
    
    !     write(*,*) 'wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat'
    !    OPEN(UNIT=11,FILE='wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',FORM='UNFORMATTED', status='OLD', IOSTAT=ios)
    !     READ(11) r_iopoint,r_time,r_dt

        jj = (i3-1)*(Nx/4)+1
        
        ! DO i2 = jj,i3*(Nx/4)
        !     !print*,i2
        !     write(*,*) 'i1,i2', i1,i2
        !     ! READ(11) (psi(i1,i2),i1=1,Nx)
        !     write(*,*) "done"
        ! ENDDO

        DO i2 = jj,i3*(Nx/4)
            !print*,i2
            ! READ(11) (psi(i1,i2),i1=1,Nx)
              OPEN(UNIT=11,FILE='wf'//TRIM(ADJUSTL(fnn))//'p'//TRIM(ADJUSTL(prcjj))//'.dat',FORM='UNFORMATTED', status='OLD', IOSTAT=ios)
              READ(11) r_iopoint,r_time,r_dt

            do i1 = 1,Nx
                 write(*,*) 'i1,i2', i1,i2
                 
                read(11) psi(i1, i2)
            enddo
            write(*,*) "done"
        ENDDO

        CLOSE(11)

    ENDDO

    write(*,*) psi
    DEALLOCATE(psi)
    ! write(*,*) "working"
    
end program main