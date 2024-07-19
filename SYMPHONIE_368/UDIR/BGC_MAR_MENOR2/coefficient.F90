       beta_=   ratio_bionegdif*0.25 
      alpha_=1.-ratio_bionegdif*0.5 


!...........................
! Constante:
!      anyv3d(i,j,k,12)=x3*(omega_w(i,j,k+1,1)-omega_w(i,j,k,1))
       anyv3d(i,j,k,12)=x3*(anyv3d(i,j,k+1,id_webio)-anyv3d(i,j,k,id_webio)) &
                          *mask_t(i,j,kmax) !01-03-19

!...........................
! En facteur de B(k+3):
        anyv3d(i,j,k,19)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(                                  &
                                  -anyv1d(k+1,4))*beta_  & 
                                                ) & !ooo>

! En facteur de B(k+2):
        anyv3d(i,j,k,18)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(                                  &
                      +anyv1d(k,4)-anyv1d(k+1,3))*beta_  & 
                     +(-anyv1d(k+1,4)           )*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k+1):
        anyv3d(i,j,k,17)=   x3*mask_t(i,j,kmax)*( & !ooo>
                    -(            -anyv1d(k+1,4)        &
                      +anyv1d(k,3)-anyv1d(k+1,2))*beta_  & 
                     +(anyv1d(k,4)-anyv1d(k+1,3))*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k  ):
        anyv3d(i,j,k,16)=1.+x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,4)-anyv1d(k+1,3)        &
                      +anyv1d(k,2)-anyv1d(k+1,1))*beta_  & 
                     +(anyv1d(k,3)-anyv1d(k+1,2))*alpha_ & 
                                                ) & !ooo>
       
! En facteur de B(k-1):
        anyv3d(i,j,k,15)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,3)-anyv1d(k+1,2)        &
                      +anyv1d(k,1)              )*beta_  & 
                     +(anyv1d(k,2)-anyv1d(k+1,1))*alpha_ & 
                                                ) & !ooo>

! En facteur de B(k-2):
        anyv3d(i,j,k,14)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k,2)-anyv1d(k+1,1)        &
                                                )*beta_  & 
                     +(anyv1d(k,1)              )*alpha_ & 
                                                ) & !ooo>
                        
! En facteur de B(k-3):
        anyv3d(i,j,k,13)=  +x3*mask_t(i,j,kmax)*( & !ooo>
                    -( anyv1d(k, 1)                     &
                                                )*beta_  & 
                                                ) & !ooo>
