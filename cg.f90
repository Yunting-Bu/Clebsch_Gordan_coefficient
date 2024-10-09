module cg_calc
    use machina_basic, only: i4, f8
    implicit none
    
contains

function factorial(f_upper_r)result(factorial_r)
    implicit none
    real(kind=f8) :: factorial_r
    integer :: loop_index       !循环指数
    integer :: f_upper
    real(kind=f8),intent(in):: f_upper_r 
    f_upper = int(f_upper_r)
    factorial_r = 1
    if (f_upper == 0) then
        factorial_r = 1     !定义：0！=1
    else
        do loop_index=1,f_upper
            factorial_r = factorial_r*loop_index
        end do
    end if 
end function factorial

subroutine dim(j1,j2,nm1,nm2,nJ,nM)
    implicit none
    
    real(kind=f8),intent(in):: j1, j2 
    integer,intent(out) :: nm1, nm2, nJ, nM 

    nm1 = int(2*j1)+1
    nm2 = int(2*j2)+1
    nJ = int((j1+j2)-abs((j1-j2)))+1
    if ((j1+j2) == nint(j1+j2)) then
        nM = nJ**2
    else 
        nM = nJ**2 + nJ 
    end if

end subroutine

subroutine get_value(j1,j2,nm1,nm2,nJ,nM,m1,m2,J,M)
    implicit none
    real(kind=f8),intent(in):: j1, j2 
    integer,intent(in) :: nm1, nm2, nJ, nM 
    real(kind=f8),allocatable,intent(out),dimension(:) :: m1,m2,J,M
    integer :: i,jj,k,l, start_pos
    real(kind=f8) :: m_val

    allocate(m1(nm1))
    allocate(m2(nm2))
    allocate(J(nJ))
    allocate(M(nM))
    do i = 1, nm1
        m1(i) = j1 - i + 1  
    end do
    do jj = 1, nm2
        m2(jj) = j2 - jj + 1  
    end do
    do k = 1, nJ
        J(k) = j1 + j2 - k + 1  
    end do
    start_pos = 1
    do k = 1, nJ
        do l = 0, int(2 * J(k))
            m_val = -J(k) + real(l, kind=f8)
            M(start_pos) = m_val
            start_pos = start_pos + 1
        end do
    end do

end subroutine



subroutine init_cg(j1,j2,nm1,nm2,nJ,nM,cg)
    implicit none
    
    real(kind=f8),intent(in):: j1, j2 
    integer,intent(in) :: nm1, nm2, nJ, nM 
    real(kind=f8),allocatable,dimension(:,:,:,:),intent(out) :: cg 
    real(kind=f8),allocatable,dimension(:) :: m1,m2,J,M
    real(kind=f8) ::cg_temp
    integer :: i,jj,k,l

    allocate(cg(nm1,nm2,nJ,nM))
    call get_value(j1,j2,nm1,nm2,nJ,nM,m1,m2,J,M)

    do i = 1, nm1 
        do jj = 1, nm2 
            do k = 1, nJ 
                do l = 1, nM 
                    if (M(l) == m1(i)+m2(jj)) then 
                        if (abs(M(l)) <= abs(J(k)) ) then
                            call calc_cg(j1,j2,m1(i),m2(jj),J(k),M(l),cg_temp)
                            cg(i,jj,k,l) = cg_temp
                        else 
                            cg(i,jj,k,l) = 0.0_f8
                        end if
                    else 
                        cg(i,jj,k,l) = 0.0_f8
                    end if
                end do 
            end do
        end do
    end do 

end subroutine

subroutine calc_cg(j1,j2,m1,m2,J,M,cg)
    implicit none
    
    real(kind=f8),intent(in) :: j1,j2,m1,m2,J,M
    real(kind=f8),intent(out) :: cg 
    integer :: z, z_min, z_max, z_temp
    real(kind=f8) :: sum

    z_temp = max(0, int(j1 - m1), int(j2 + m2), int(j1 + j2 - J))
    z_max = min(int(j1 + j2 - J), int(j1 - m1), int(j2 + m2))
    z_min = min(z_max,z_temp)
    sum = 0.0
    do z = z_min, z_max
        sum = sum+(-1)**z*(1/(factorial(real(z,kind=f8))*factorial(j1+j2-J-z)*factorial(j1-m1-z)*factorial(j2+m2-z) &
                            * factorial(J-j2+m1+z)*factorial(J-j1-m2+z)))
    end do
    cg = sqrt((2*J+1)*factorial(j1+j2-J)*factorial(J+j2-j1)*factorial(J+j1-j2)/factorial(J+j1+j2+1)) &
         * sqrt(real(factorial(J+M)*factorial(J-M)*factorial(j1+m1)*factorial(j1-m1)*factorial(j2+m2)*factorial(j2-m2))) &
         * sum
end subroutine


                
    
end module cg_calc

