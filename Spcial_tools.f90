!一维数组线性插值
SUBROUTINE Linear_Interpolation(array1,array2)
    IMPLICIT NONE
    
    REAL,ALLOCATABLE::array1(:),array2(:)
    INTEGER::size1,size2,i,j,k
    
    size1=SIZE(array1)
    size2=SIZE(array2)
    k=(size2-1)/(size1-1)
    IF(size1/=size2)THEN
        DO i=1,size2-2
            array2(i)=array1(i/k+1)+(MOD(i,k)-1)*(array1(i/k+2)-array1(i/k+1))/k
        END DO
        array2(size2-1)=array1(size1)-(array1(size1)-array1(size1-1))/k
        array2(size2)=array1(size1)
    ELSE
        array2=array1
    END IF
    
END
!二维数组线性插值
SUBROUTINE Linear_Interpolation_2d(array1,array2,direction)
    IMPLICIT NONE
    
    REAL,ALLOCATABLE::array1(:,:),array2(:,:)
    INTEGER::direction
    INTEGER::size1,size2,size3,i,j,k
    
    IF(direction==1)THEN
        size1=SIZE(array1(:,1))
        size2=SIZE(array2(:,1))
        size3=SIZE(array1(1,:))
        k=(size2-1)/(size1-1)
        DO j=1,size3
            IF(size1/=size2)THEN
                DO i=1,size2-2
                    array2(i,j)=array1(i/k+1,j)+(MOD(i,k)-1)*(array1(i/k+2,j)-array1(i/k+1,j))/k
                END DO
                array2(size2-1,j)=array1(size1,j)-(array1(size1,j)-array1(size1-1,j))/k
                array2(size2,j)=array1(size1,j)
            ELSE
                array2=array1
            END IF
        END DO
    ELSE IF(direction==2)THEN
        size1=SIZE(array1(1,:))
        size2=SIZE(array2(1,:))
        size3=SIZE(array1(:,1))
        k=(size2-1)/(size1-1)
        DO j=1,size3
            IF(size1/=size2)THEN
                DO i=1,size2-2
                    array2(j,i)=array1(j,i/k+1)+(MOD(i,k)-1)*(array1(j,i/k+2)-array1(j,i/k+1))/k
                END DO
                array2(j,size2-1)=array1(j,size1)-(array1(j,size1)-array1(j,size1-1))/k
                array2(j,size2)=array1(j,size1)
            ELSE
                array2=array1
            END IF
        END DO
    END IF
    
END