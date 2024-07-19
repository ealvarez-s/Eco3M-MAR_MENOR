      program toto
      INTEGER :: n,errorflag,i,j,k,kmax=30,kmin=1
      REAL*8, DIMENSION(:,:), allocatable :: inputmat,invrsmat,matprod &
       ,mat2dv,transpmat2dv
      REAL*8, DIMENSION(:), allocatable :: profil,depth,coef,matlhs,profil2
      REAL*8 sum1

      n=3
      allocate(    inputmat(n,n)   ) ; inputmat=0.
      allocate(    invrsmat(n,n)   ) ; invrsmat=0.
      allocate(     matprod(n,n)   ) ; matprod=0.
      allocate(      mat2dv(kmax,n)) ; mat2dv=0.
      allocate(transpmat2dv(n,kmax)) ; transpmat2dv=0.
      allocate(        profil(kmax)) ; profil=0.
      allocate(       profil2(kmax)) ; profil2=0.
      allocate(           matlhs(n)) ; matlhs=0.
      allocate(         depth(kmax)) ; depth=0.
      allocate(             coef(n)) ; coef=1.

      do k=kmax,1,-1
       depth(k)=3000.*real(k-kmax)/real(kmax-1)
       profil(k)=exp(depth(k)/5000.)
!      profil(k)=1.+depth(k)+depth(k)**2
!      write(66,*)depth(k),profil(k)
      enddo

      do i=1,n
       do k=kmax,kmin,-1
!       mat2dv(k,i)=depth(k)**(i-1)
        mat2dv(k,i)=exp(depth(k)/(50.+(i-1)*1000.))
       enddo
      enddo

!     matlhs=matmul(mat2dv,coef)
!     do k=kmax,1,-1
!      write(6,*)profil(k)
!     enddo
      

      transpmat2dv=transpose(mat2dv)
        matlhs=matmul(transpmat2dv,profil)
      inputmat=matmul(transpmat2dv,mat2dv)
      call s_matrix_reverser(inputmat,invrsmat,n,errorflag)

      write(6,*)'inputmat--------------'
      write(6,*)inputmat
      write(6,*)'matlhs--------------'
      write(6,*)matlhs
      write(6,*)'invrsmat-------------'
      write(6,*)invrsmat

      coef=matmul(invrsmat,matlhs)
      write(6,*)'coef---------'
      write(6,*)coef

      profil2=matmul(mat2dv,coef)
      write(6,*)'profil2----------'
      write(6,*)profil2
      
!     write(6,*)'-----------------------------------'
      do k=kmax,1,-1
      write(66,*)depth(k),profil(k),profil2(k)
      enddo
     
      

      ! Inverse Mat
!     do i=1,n
!     write(6,*)invrsmat(i,:)
!     enddo


      matprod=matmul(inputmat,invrsmat)

      ! Product Mat:
      write(6,*)'-------------------'
      do i=1,n
       write(6,*)matprod(i,:)
      enddo

      end

!............................................................

subroutine s_matrix_reverser(matrix,inverse,n,errorflag)
IMPLICIT NONE
!Declarations
INTEGER :: n
INTEGER :: errorflag !Return error status. ­1 for error, 0 for normal
REAL*8, DIMENSION(n,n) :: matrix !Input matrix
REAL*8, DIMENSION(n,n) :: inverse !Inverted matrix
LOGICAL :: FLAG = .TRUE.
INTEGER :: i, j, k, l
REAL*8 :: m
REAL*8, DIMENSION(n,2*n) :: augmatrix !augmented matrix

!Augment input matrix with an identity matrix
DO i=1,n
DO j=1,2*n
IF (j<=n) THEN
augmatrix(i,j) = matrix(i,j)
ELSE IF ((i+n) == j) THEN
augmatrix(i,j) = 1
Else
augmatrix(i,j) = 0
ENDIF
END DO
END DO

!Reduce augmented matrix to upper traingular form
DO k=1,n-1
IF (augmatrix(k,k) == 0) THEN
FLAG = .FALSE.
DO i = k+1, n
IF (augmatrix(i,k) /= 0) THEN
DO j = 1,2*n
augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
END DO
FLAG = .TRUE.
EXIT
ENDIF

IF (FLAG .EQV. .FALSE.) THEN
PRINT*, "Matrix is non ­ invertible"
inverse=0
errorflag=-1
return
ENDIF
END DO
ENDIF

DO j = k+1, n
m = augmatrix(j,k)/augmatrix(k,k)
DO i = k, 2*n
augmatrix(j,i)=augmatrix(j,i)-m*augmatrix(k,i)
END DO
END DO
END DO

!Test for invertibility

DO i=1,n
IF (augmatrix(i,i) == 0) THEN
PRINT*, "Matrix is non ­ invertible"
inverse=0
errorflag=-1
return
ENDIF
END DO

!Make diagonal elements as 1

do i=1,n
m=augmatrix(i,i)
DO j=i,(2*n)
augmatrix(i,j)=(augmatrix(i,j)/m)
END DO
END DO

!Reduced right side half of augmented matrix to identity matrix
DO k =n-1,1,-1
DO i=1,k
m=augmatrix(i,k+1)
DO j=k,(2*n)
augmatrix(i,j)=augmatrix(i,j)-augmatrix(k+1,j)*m
END DO
END DO
END DO

!store answer

DO i =1, n
DO j = 1, n
inverse(i,j) = augmatrix(i,j+n)
END DO
END DO
errorflag = 0

end subroutine s_matrix_reverser
