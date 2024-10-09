program main 
use machina_basic
use cg_calc
    implicit none

    real(kind=f8):: j1, j2 
    integer :: nm1, nm2, nJ, nM
    real(kind=f8),allocatable,dimension(:,:,:,:) :: cg 
    real(kind=f8),allocatable,dimension(:) :: m1,m2,J,M
    integer :: i,jj,k,l
    real(kind=f8) :: last_m1, last_m2, last_J, last_M, last_cg
    last_m1 = -999.0
    last_m2 = -999.0
    last_J = -999.0
    last_M = -999.0
    last_cg = -999.0

    write(*,'(a)') 'Please input j1 and j2: '
    read(*,*) j1,j2 
    call dim(j1,j2,nm1,nm2,nJ,nM)
    call get_value(j1,j2,nm1,nm2,nJ,nM,m1,m2,J,M)
    call init_cg(j1,j2,nm1,nm2,nJ,nM,cg) 

    write(*,'(/,a)') 'The Clebsch-Gordan coefficients are calculated using the Racah formula.'
    write(*,'(/,a)') '< m1 m2 | J M >'
    write(*,'(a,/)') '======================================='
    do i = 1, nm1 
        do jj = 1, nm2 
            do k = 1, nJ 
                do l = 1, nM 
                    if (cg(i,jj,k,l)/=0.0) then
                        if (m1(i) /= last_m1 .or. m2(jj) /= last_m2 .or. J(k) /= last_J .or. M(l) &
                        /= last_M .or. cg(i,jj,k,l) /= last_cg) then
                            write(*,'(a,f5.1,f5.1,a,f5.1,f5.1,a,f10.5)') '<', m1(i), m2(jj), &
                            '  |', J(k), M(l), '  > = ', cg(i,jj,k,l)
                            last_m1 = m1(i)
                            last_m2 = m2(jj)
                            last_J = J(k)
                            last_M = M(l)
                            last_cg = cg(i,jj,k,l)
                        end if
                    end if
                end do 
            end do 
        end do 
    end do
    


end program main 