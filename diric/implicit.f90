!#define my_mod_solide

module ludecomp

   use decomp_2d, only : mytype

   implicit none

   private
   public :: ludecomp7

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Interface for septadiag LU decomposition used in implicit mode
   ! ludecomp7_0 is called when ncly=0
   ! ludecomp7_12 is called when ncly=1 or 2
   !
   interface ludecomp7
      module procedure ludecomp7_12
      module procedure ludecomp7_0
   end interface ludecomp7

   contains

   !*******************************************************************
   !Decomposition LU matrice septadiag
   subroutine ludecomp7_12(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ny)
   ! 
   !*******************************************************************
      use decomp_2d, only : mytype

      implicit none

      integer :: j
      integer, intent(in) :: ny
      real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
      real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
      real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm

      ggm=0.;hhm=0.;ssm=0.
      vvm=0.;wwm=0.;zzm=0.
  
      ggm(1)=aam(1)
      hhm(1)=bbm(1)
      ssm(1)=ccm(1)
  
      zzm(2)=ddm(2)/ggm(1)
      ggm(2)=aam(2)-zzm(2)*hhm(1)
      hhm(2)=bbm(2)-zzm(2)*ssm(1)
      ssm(2)=ccm(2)-zzm(2)*rrm(1)

      wwm(3)=eem(3)/ggm(1)
      zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
      ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
      hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
      ssm(3)=ccm(3)-zzm(3)*rrm(2)
  
      do j=4,ny
         vvm(j)=qqm(j)/ggm(j-3)
         wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
         zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
         ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
         hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
         ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
      enddo

   end subroutine ludecomp7_12
   !

   !*******************************************************************
   !Decomposition LU matrice cyclique septadiag
   subroutine ludecomp7_0(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,ny)
   ! 
   !*******************************************************************
      use decomp_2d, only : mytype

      implicit none

      integer :: j,k
      integer, intent(in) :: ny
      real(mytype),dimension(ny), intent(in)  :: aam,bbm,ccm,ddm,eem,rrm,qqm
      real(mytype),dimension(ny), intent(out) :: ggm,hhm,ssm
      real(mytype),dimension(ny), intent(out) :: vvm,wwm,zzm
      real(mytype),dimension(ny), intent(out) :: l1m,l2m,l3m
      real(mytype),dimension(ny), intent(out) :: u1m,u2m,u3m

      ! a=diag, b=diag+1, c=diag+2, r=diag+3
      !         d=diag-1, e=diag-2, q=diag-3
      ggm=0.;hhm=0.;ssm=0. ! U part (diag, diag+1, diag+2, diag+3=rrm)
      u1m=0.;u2m=0.;u3m=0. ! U cyclic part (last col, last-1, last-2)
      vvm=0.;wwm=0.;zzm=0. ! L part (diag=1, diag-3, diag-2, diag-1)
      l1m=0.;l2m=0.;l3m=0. ! L cyclic part (last row, last-1, last-2)

      ggm(1)=aam(1)
      hhm(1)=bbm(1)
      ssm(1)=ccm(1)
      u1m(1)=ddm(1)
      u2m(1)=eem(1)
      u3m(1)=qqm(1)
      l1m(1)=bbm(ny)/ggm(1)
      l2m(1)=ccm(ny-1)/ggm(1)
      l3m(1)=rrm(ny-2)/ggm(1)

      zzm(2)=ddm(2)/ggm(1)
      ggm(2)=aam(2)-zzm(2)*hhm(1)
      hhm(2)=bbm(2)-zzm(2)*ssm(1)
      ssm(2)=ccm(2)-zzm(2)*rrm(1)
      u1m(2)=eem(2)-zzm(2)*u1m(1)
      u2m(2)=qqm(2)-zzm(2)*u2m(1)
      u3m(2)=      -zzm(2)*u3m(1)
      l1m(2)=(ccm(ny)  -hhm(1)*l1m(1))/ggm(2)
      l2m(2)=(rrm(ny-1)-hhm(1)*l2m(1))/ggm(2)
      l3m(2)=(         -hhm(1)*l3m(1))/ggm(2)

      wwm(3)=eem(3)/ggm(1)
      zzm(3)=(ddm(3)-wwm(3)*hhm(1))/ggm(2)
      ggm(3)=aam(3)-wwm(3)*ssm(1)-zzm(3)*hhm(2)
      hhm(3)=bbm(3)-wwm(3)*rrm(1)-zzm(3)*ssm(2)
      ssm(3)=ccm(3)              -zzm(3)*rrm(2)
      u1m(3)=qqm(3)-wwm(3)*u1m(1)-zzm(3)*u1m(2)
      u2m(3)=      -wwm(3)*u2m(1)-zzm(3)*u2m(2)
      u3m(3)=      -wwm(3)*u3m(1)-zzm(3)*u3m(2)
      l1m(3)=(rrm(ny)-ssm(1)*l1m(1)-hhm(2)*l1m(2))/ggm(3)
      l2m(3)=(       -ssm(1)*l2m(1)-hhm(2)*l2m(2))/ggm(3)
      l3m(3)=(       -ssm(1)*l3m(1)-hhm(2)*l3m(2))/ggm(3)
  
      do j=4,ny-3
         vvm(j)=qqm(j)/ggm(j-3)
         wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
         zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
         ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
         hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
         ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
         u1m(j)= - vvm(j)*u1m(j-3) - wwm(j)*u1m(j-2) - zzm(j)*u1m(j-1)
         u2m(j)= - vvm(j)*u2m(j-3) - wwm(j)*u2m(j-2) - zzm(j)*u2m(j-1)
         u3m(j)= - vvm(j)*u3m(j-3) - wwm(j)*u3m(j-2) - zzm(j)*u3m(j-1)
         l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1))/ggm(j)
         l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1))/ggm(j)
         l3m(j)= - (rrm(j-3)*l3m(j-3)+ssm(j-2)*l3m(j-2)+hhm(j-1)*l3m(j-1))/ggm(j)
      enddo

      do j=ny-2,ny
         vvm(j)=qqm(j)/ggm(j-3)
         wwm(j)=(eem(j)-vvm(j)*hhm(j-3))/ggm(j-2)
         zzm(j)=(ddm(j)-wwm(j)*hhm(j-2)-vvm(j)*ssm(j-3))/ggm(j-1)
         ggm(j)=aam(j)-vvm(j)*rrm(j-3)-wwm(j)*ssm(j-2)-zzm(j)*hhm(j-1)
         hhm(j)=bbm(j)-wwm(j)*rrm(j-2)-zzm(j)*ssm(j-1)
         ssm(j)=ccm(j)-zzm(j)*rrm(j-1)
      enddo
      j=ny-2
      u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
                - l3m(j-1)*rrm(j-1)
      u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
                - l3m(j-2)*rrm(j-2) - l3m(j-1)*ssm(j-1)
      u3m(j)= - vvm(j)*u3m(j-3)-wwm(j)*u3m(j-2)-zzm(j)*u3m(j-1) &
                - l3m(j-3)*rrm(j-3) - l3m(j-2)*ssm(j-2) - l3m(j-1)*hhm(j-1)
      do k=1,j-1
         u1m(j)=u1m(j) - u1m(k)*l3m(k)
         u2m(j)=u2m(j) - u2m(k)*l3m(k)
         u3m(j)=u3m(j) - u3m(k)*l3m(k)
      enddo
      l1m(j)= - ( rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
         + vvm(ny)*u3m(j-1)+wwm(ny)*u3m(j) ) / ( ggm(j)+u3m(j) )
      l2m(j)= - (rrm(j-3)*l2m(j-3)+ssm(j-2)*l2m(j-2)+hhm(j-1)*l2m(j-1) &
         + vvm(ny-1)*u3m(j-2)+wwm(ny-1)*u3m(j-1)+zzm(ny-1)*u3m(j) ) / (ggm(j)+u3m(j))
      do k=1,j-1
         l1m(j)=l1m(j) - l1m(k)*u3m(k) / (ggm(j)+u3m(j))
         l2m(j)=l2m(j) - l2m(k)*u3m(k) / (ggm(j)+u3m(j))
      enddo
      j=ny-1
      u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
                 - l2m(j-2)*rrm(j-2) - l2m(j-1)*ssm(j-1)
      u2m(j)= - vvm(j)*u2m(j-3)-wwm(j)*u2m(j-2)-zzm(j)*u2m(j-1) &
                 - l2m(j-3)*rrm(j-3) - l2m(j-2)*ssm(j-2) - l2m(j-1)*hhm(j-1)
      do k=1,j-1
         u1m(j)=u1m(j) - u1m(k)*l2m(k)
         u2m(j)=u2m(j) - u2m(k)*l2m(k)
      enddo
      l1m(j)= - (rrm(j-3)*l1m(j-3)+ssm(j-2)*l1m(j-2)+hhm(j-1)*l1m(j-1) &
         + vvm(ny)*u2m(j-2)+wwm(ny)*u2m(j-1)+zzm(ny)*u2m(j) ) / (ggm(j)+u2m(j));
      do k=1,j-1
         l1m(j)=l1m(j) - l1m(k)*u2m(k) / (ggm(j)+u2m(j));
      enddo
      j=ny
      u1m(j)= - vvm(j)*u1m(j-3)-wwm(j)*u1m(j-2)-zzm(j)*u1m(j-1) &
                 - l1m(j-3)*rrm(j-3)-l1m(j-2)*ssm(j-2)-l1m(j-1)*hhm(j-1);
      do k=1,j-1
         u1m(j)=u1m(j) - u1m(k)*l1m(k)
      enddo

end subroutine ludecomp7_0
   !

end module ludecomp

module matinv

   use decomp_2d, only : mytype

   implicit none

   private
   public :: septinv, xseptinv, zseptinv

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Interface for septadiag matrix inversion used in implicit mode
   ! septinv_0 is called when ncly=0
   ! septinv_12 is called when ncly=1 or 2
   !
   interface septinv
      module procedure septinv_12
      module procedure septinv_0
   end interface septinv
   !
   ! xseptinv is called for x inversion
   !
   interface xseptinv
      module procedure xseptinv_12
      module procedure xseptinv_0
   end interface xseptinv
   !
   ! zseptinv is called for z inversion
   !
   interface zseptinv
      module procedure zseptinv_12
      module procedure zseptinv_0
   end interface zseptinv

   contains

   !********************************************************************
   !INVERSE SEPTADIAG
   subroutine septinv_12(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,nx,ny,nz)
   ! 
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none
  
      integer :: j
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
      real(mytype),dimension(ny),intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(ny),intent(in) :: vvm,wwm,zzm

      !Descente
      xsol(:,1,:)=bbb(:,1,:)
      xsol(:,2,:)=bbb(:,2,:)-zzm(2)*xsol(:,1,:)
      xsol(:,3,:)=bbb(:,3,:)-zzm(3)*xsol(:,2,:) &
                            -wwm(3)*xsol(:,1,:)

      do j=4,ny
         xsol(:,j,:)=bbb(:,j,:)-vvm(j)*xsol(:,j-3,:) &
                               -wwm(j)*xsol(:,j-2,:) &
                               -zzm(j)*xsol(:,j-1,:)
      enddo

      !Remontee
      xsol(:,ny  ,:)= xsol(:,ny  ,:)/ggm(ny)
      xsol(:,ny-1,:)=(xsol(:,ny-1,:)-hhm(ny-1)*xsol(:,ny  ,:))/ggm(ny-1)
      xsol(:,ny-2,:)=(xsol(:,ny-2,:)-hhm(ny-2)*xsol(:,ny-1,:) &
                                    -ssm(ny-2)*xsol(:,ny  ,:))/ggm(ny-2)

      do j=ny-3,1,-1
         xsol(:,j,:)=(xsol(:,j,:)-hhm(j)*xsol(:,j+1,:) &
                                 -ssm(j)*xsol(:,j+2,:) &
                                 -rrm(j)*xsol(:,j+3,:))/ggm(j)
      enddo

   end subroutine septinv_12
   !
   subroutine xseptinv_12(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,nx,ny,nz)
   !
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none
  
      integer :: i
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
      real(mytype),dimension(nx),intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(nx),intent(in) :: vvm,wwm,zzm

      !Descente
      xsol(1,:,:)=bbb(1,:,:)
      xsol(2,:,:)=bbb(2,:,:)-zzm(2)*xsol(1,:,:)
      xsol(3,:,:)=bbb(3,:,:)-zzm(3)*xsol(2,:,:) &
                            -wwm(3)*xsol(1,:,:)

      do i=4,nx
         xsol(i,:,:)=bbb(i,:,:)-vvm(i)*xsol(i-3,:,:) &
                               -wwm(i)*xsol(i-2,:,:) &
                               -zzm(i)*xsol(i-1,:,:)
      enddo

      !Remontee
      xsol(nx  ,:,:)= xsol(nx  ,:,:)/ggm(nx)
      xsol(nx-1,:,:)=(xsol(nx-1,:,:)-hhm(nx-1)*xsol(nx  ,:,:))/ggm(nx-1)
      xsol(nx-2,:,:)=(xsol(nx-2,:,:)-hhm(nx-2)*xsol(nx-1,:,:) &
                                    -ssm(nx-2)*xsol(nx  ,:,:))/ggm(nx-2)

      do i=nx-3,1,-1
         xsol(i,:,:)=(xsol(i,:,:)-hhm(i)*xsol(i+1,:,:) &
                                 -ssm(i)*xsol(i+2,:,:) &
                                 -rrm(i)*xsol(i+3,:,:))/ggm(i)
      enddo

   end subroutine xseptinv_12
   !
   subroutine zseptinv_12(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,nx,ny,nz)
   ! 
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none
  
      integer :: k
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz),intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz),intent(in) :: bbb
      real(mytype),dimension(nz),intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(nz),intent(in) :: vvm,wwm,zzm

      !Descente
      xsol(:,:,1)=bbb(:,:,1)
      xsol(:,:,2)=bbb(:,:,2)-zzm(2)*xsol(:,:,1)
      xsol(:,:,3)=bbb(:,:,3)-zzm(3)*xsol(:,:,2) &
                            -wwm(3)*xsol(:,:,1)

      do k=4,nz
         xsol(:,:,k)=bbb(:,:,k)-vvm(k)*xsol(:,:,k-3) &
                               -wwm(k)*xsol(:,:,k-2) &
                               -zzm(k)*xsol(:,:,k-1)
      enddo

      !Remontee
      xsol(:,:,nz  )= xsol(:,:,nz  )/ggm(nz)
      xsol(:,:,nz-1)=(xsol(:,:,nz-1)-hhm(nz-1)*xsol(:,:,nz  ))/ggm(nz-1)
      xsol(:,:,nz-2)=(xsol(:,:,nz-2)-hhm(nz-2)*xsol(:,:,nz-1) &
                                    -ssm(nz-2)*xsol(:,:,nz  ))/ggm(nz-2)

      do k=nz-3,1,-1
         xsol(:,:,k)=(xsol(:,:,k)-hhm(k)*xsol(:,:,k+1) &
                                 -ssm(k)*xsol(:,:,k+2) &
                                 -rrm(k)*xsol(:,:,k+3))/ggm(k)
      enddo

   end subroutine zseptinv_12

   !********************************************************************
   !INVERSE SEPTADIAG CYCLIQUE
   subroutine septinv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)
   ! 
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none

      integer :: j,jj
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
      real(mytype),dimension(ny), intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(ny), intent(in) :: vvm,wwm,zzm
      real(mytype),dimension(ny), intent(in) :: l1m,l2m,l3m
      real(mytype),dimension(ny), intent(in) :: u1m,u2m,u3m

      !Descente avec matrice triangulaire inf
      xsol(:,1,:)=bbb(:,1,:)
      xsol(:,2,:)=bbb(:,2,:)-zzm(2)*xsol(:,1,:)
      xsol(:,3,:)=bbb(:,3,:)-zzm(3)*xsol(:,2,:) &
                            -wwm(3)*xsol(:,1,:)
      do j=4,ny-3
         xsol(:,j,:)=bbb(:,j,:)-vvm(j)*xsol(:,j-3,:) &
                               -wwm(j)*xsol(:,j-2,:) &
                               -zzm(j)*xsol(:,j-1,:)
      enddo
      j=ny-2
      xsol(:,j,:)=bbb(:,j,:)-vvm(j)*xsol(:,j-3,:) &
                            -wwm(j)*xsol(:,j-2,:) &
                            -zzm(j)*xsol(:,j-1,:)
      do jj=1,j-1
         xsol(:,j,:)=xsol(:,j,:)-l3m(jj)*xsol(:,jj,:)
      enddo
      j=ny-1
      xsol(:,j,:)=bbb(:,j,:)-vvm(j)*xsol(:,j-3,:) &
                            -wwm(j)*xsol(:,j-2,:) &
                            -zzm(j)*xsol(:,j-1,:)
      do jj=1,j-1
         xsol(:,j,:)=xsol(:,j,:)-l2m(jj)*xsol(:,jj,:)
      enddo
      j=ny
      xsol(:,j,:)=bbb(:,j,:)-vvm(j)*xsol(:,j-3,:) &
                            -wwm(j)*xsol(:,j-2,:) &
                            -zzm(j)*xsol(:,j-1,:)
      do jj=1,j-1
         xsol(:,j,:)=xsol(:,j,:)-l1m(jj)*xsol(:,jj,:)
      enddo
 
      !Remontee avec matrice triangulaire sup.
      xsol(:,ny  ,:)= xsol(:,ny  ,:)/(ggm(ny)+u1m(ny))
      xsol(:,ny-1,:)=(xsol(:,ny-1,:)-(hhm(ny-1)+u1m(ny-1))*xsol(:,ny  ,:))/(ggm(ny-1)+u2m(ny-1))
      xsol(:,ny-2,:)=(xsol(:,ny-2,:)-(hhm(ny-2)+u2m(ny-2))*xsol(:,ny-1,:) &
                                    -(ssm(ny-2)+u1m(ny-2))*xsol(:,ny  ,:))/(ggm(ny-2)+u3m(ny-2))

      do j=ny-3,1,-1
         xsol(:,j,:)=(xsol(:,j,:)-hhm(j)*xsol(:,j+1 ,:) &
                                 -ssm(j)*xsol(:,j+2 ,:) &
                                 -rrm(j)*xsol(:,j+3 ,:) & 
                                 -u3m(j)*xsol(:,ny-2,:) &
                                 -u2m(j)*xsol(:,ny-1,:) &
                                 -u1m(j)*xsol(:,ny  ,:) )/ggm(j)
      enddo

   end subroutine septinv_0
   !
   subroutine xseptinv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)
   ! 
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none

      integer :: i,ii
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
      real(mytype),dimension(nx), intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(nx), intent(in) :: vvm,wwm,zzm
      real(mytype),dimension(nx), intent(in) :: l1m,l2m,l3m
      real(mytype),dimension(nx), intent(in) :: u1m,u2m,u3m

      !Descente avec matrice triangulaire inf
      xsol(1,:,:)=bbb(1,:,:)
      xsol(2,:,:)=bbb(2,:,:)-zzm(2)*xsol(1,:,:)
      xsol(3,:,:)=bbb(3,:,:)-zzm(3)*xsol(2,:,:) &
                            -wwm(3)*xsol(1,:,:)
      do i=4,nx-3
         xsol(i,:,:)=bbb(i,:,:)-vvm(i)*xsol(i-3,:,:) &
                               -wwm(i)*xsol(i-2,:,:) &
                               -zzm(i)*xsol(i-1,:,:)
      enddo
      i=nx-2
      xsol(i,:,:)=bbb(i,:,:)-vvm(i)*xsol(i-3,:,:) &
                            -wwm(i)*xsol(i-2,:,:) &
                            -zzm(i)*xsol(i-1,:,:)
      do ii=1,i-1
         xsol(i,:,:)=xsol(i,:,:)-l3m(ii)*xsol(ii,:,:)
      enddo
      i=nx-1
      xsol(i,:,:)=bbb(i,:,:)-vvm(i)*xsol(i-3,:,:) &
                            -wwm(i)*xsol(i-2,:,:) &
                            -zzm(i)*xsol(i-1,:,:)
      do ii=1,i-1
         xsol(i,:,:)=xsol(i,:,:)-l2m(ii)*xsol(ii,:,:)
      enddo
      i=nx
      xsol(i,:,:)=bbb(i,:,:)-vvm(i)*xsol(i-3,:,:) &
                            -wwm(i)*xsol(i-2,:,:) &
                            -zzm(i)*xsol(i-1,:,:)
      do ii=1,i-1
         xsol(i,:,:)=xsol(i,:,:)-l1m(ii)*xsol(ii,:,:)
      enddo

      !Remontee avec matrice triangulaire sup.
      xsol(nx  ,:,:)= xsol(nx  ,:,:)/(ggm(nx)+u1m(nx))
      xsol(nx-1,:,:)=(xsol(nx-1,:,:)-(hhm(nx-1)+u1m(nx-1))*xsol(nx  ,:,:))/(ggm(nx-1)+u2m(nx-1))
      xsol(nx-2,:,:)=(xsol(nx-2,:,:)-(hhm(nx-2)+u2m(nx-2))*xsol(nx-1,:,:) &
                                    -(ssm(nx-2)+u1m(nx-2))*xsol(nx  ,:,:))/(ggm(nx-2)+u3m(nx-2))

      do i=nx-3,1,-1
         xsol(i,:,:)=(xsol(i,:,:)-hhm(i)*xsol(i+1 ,:,:) &
                                 -ssm(i)*xsol(i+2 ,:,:) &
                                 -rrm(i)*xsol(i+3 ,:,:) &
                                 -u3m(i)*xsol(nx-2,:,:) &
                                 -u2m(i)*xsol(nx-1,:,:) &
                                 -u1m(i)*xsol(nx  ,:,:) )/ggm(i)
      enddo

   end subroutine xseptinv_0
   !
   subroutine zseptinv_0(xsol,bbb,ggm,hhm,ssm,rrm,vvm,wwm,zzm,l1m,l2m,l3m,u1m,u2m,u3m,nx,ny,nz)
   ! 
   !********************************************************************
      use decomp_2d, only : mytype

      implicit none

      integer :: k,kk
      integer, intent(in) :: nx,ny,nz
      real(mytype),dimension(nx,ny,nz), intent(out) :: xsol
      real(mytype),dimension(nx,ny,nz), intent(in) :: bbb
      real(mytype),dimension(nz), intent(in) :: ggm,hhm,ssm,rrm
      real(mytype),dimension(nz), intent(in) :: vvm,wwm,zzm
      real(mytype),dimension(nz), intent(in) :: l1m,l2m,l3m
      real(mytype),dimension(nz), intent(in) :: u1m,u2m,u3m

      !Descente avec matrice triangulaire inf
      xsol(:,:,1)=bbb(:,:,1)
      xsol(:,:,2)=bbb(:,:,2)-zzm(2)*xsol(:,:,1)
      xsol(:,:,3)=bbb(:,:,3)-zzm(3)*xsol(:,:,2) &
                            -wwm(3)*xsol(:,:,1)
      do k=4,nz-3
         xsol(:,:,k)=bbb(:,:,k)-vvm(k)*xsol(:,:,k-3) &
                               -wwm(k)*xsol(:,:,k-2) &
                               -zzm(k)*xsol(:,:,k-1)
      enddo
      k=nz-2
      xsol(:,:,k)=bbb(:,:,k)-vvm(k)*xsol(:,:,k-3) &
                            -wwm(k)*xsol(:,:,k-2) &
                            -zzm(k)*xsol(:,:,k-1)
      do kk=1,k-1
         xsol(:,:,k)=xsol(:,:,k)-l3m(kk)*xsol(:,:,kk)
      enddo
      k=nz-1
      xsol(:,:,k)=bbb(:,:,k)-vvm(k)*xsol(:,:,k-3) &
                            -wwm(k)*xsol(:,:,k-2) &
                            -zzm(k)*xsol(:,:,k-1)
      do kk=1,k-1
         xsol(:,:,k)=xsol(:,:,k)-l2m(kk)*xsol(:,:,kk)
      enddo
      k=nz
      xsol(:,:,k)=bbb(:,:,k)-vvm(k)*xsol(:,:,k-3) &
                            -wwm(k)*xsol(:,:,k-2) &
                            -zzm(k)*xsol(:,:,k-1)
      do kk=1,k-1
         xsol(:,:,k)=xsol(:,:,k)-l1m(kk)*xsol(:,:,kk)
      enddo

      !Remontee avec matrice triangulaire sup.
      xsol(:,:,nz  )= xsol(:,:,nz  )/(ggm(nz)+u1m(nz))
      xsol(:,:,nz-1)=(xsol(:,:,nz-1)-(hhm(nz-1)+u1m(nz-1))*xsol(:,:,nz  ))/(ggm(nz-1)+u2m(nz-1))
      xsol(:,:,nz-2)=(xsol(:,:,nz-2)-(hhm(nz-2)+u2m(nz-2))*xsol(:,:,nz-1) &
                                    -(ssm(nz-2)+u1m(nz-2))*xsol(:,:,nz  ))/(ggm(nz-2)+u3m(nz-2))

      do k=nz-3,1,-1
         xsol(:,:,k)=(xsol(:,:,k)-hhm(k)*xsol(:,:,k+1 ) &
                                 -ssm(k)*xsol(:,:,k+2 ) &
                                 -rrm(k)*xsol(:,:,k+3 ) & 
                                 -u3m(k)*xsol(:,:,nz-2) &
                                 -u2m(k)*xsol(:,:,nz-1) &
                                 -u1m(k)*xsol(:,:,nz  ) )/ggm(k)
      enddo

   end subroutine zseptinv_0

end module matinv

!********************************************************************
!
subroutine  inttimp (ux1,uy1,uz1,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1,px1,py1,pz1,&
     td1,te1,tf1,tg1,th1,ti1,di1,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,di2,&
     ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,di3)
! 
!********************************************************************
USE param
USE variables
USE decomp_2d
use derivY
use matinv

implicit none

integer :: ijk,nxyz,i,j,k,code
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: td1,te1,tf1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux2,uy2,uz2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ux3,uy3,uz3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: td3,te3,tf3,di3
real(mytype) :: x, y, z

td1=0.;te1=0.;tf1=0.;tg1=0.;th1=0.;ti1=0.;di1=0.
ta2=0.;tb2=0.;tc2=0.;td2=0.;te2=0.;tf2=0.;di2=0.
ta3=0.;tb3=0.;tc3=0.;td3=0.;te3=0.;tf3=0.;di3=0.

nxyz=xsize(1)*xsize(2)*xsize(3)

if ((nscheme.eq.2).or.(nscheme.eq.3)) then
   if (nrank.eq.0) print *,'RK scheme not ready, abort computation'
   call decomp_2d_finalize
   CALL MPI_FINALIZE(code)
   call exit(0)
endif

!!!!!!!!!!!!!
!! CN2+AB2 !!
!!!!!!!!!!!!!
if (nscheme.eq.1) then
   if ((itime.eq.1).and.(ilit.eq.0)) then

      !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES
      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               !uhat
               td1(i,j,k)= gdt(itr)*ta1(i,j,k)
               te1(i,j,k)= gdt(itr)*tb1(i,j,k) 
               tf1(i,j,k)= gdt(itr)*tc1(i,j,k)
               
               !save (N+Lxz)(n-1)
               gx(i,j,k)=ta1(i,j,k)
               gy(i,j,k)=tb1(i,j,k)
               gz(i,j,k)=tc1(i,j,k)
               
            enddo
         enddo
      enddo

   else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               !uhat
               td1(i,j,k)= adt(itr)*ta1(i,j,k)+bdt(itr)*gx(i,j,k)-px1(i,j,k)
               te1(i,j,k)= adt(itr)*tb1(i,j,k)+bdt(itr)*gy(i,j,k)-py1(i,j,k)
               tf1(i,j,k)= adt(itr)*tc1(i,j,k)+bdt(itr)*gz(i,j,k)-pz1(i,j,k)
               
               !save (N+Lxz)(n-1)
               gx(i,j,k)=ta1(i,j,k)
               gy(i,j,k)=tb1(i,j,k)
               gz(i,j,k)=tc1(i,j,k)  
            enddo
         enddo
      enddo

   endif
endif

!!!!!!!!!!!!!!!
!!! CN2+AB3 !!!
!!!!!!!!!!!!!!!
if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then

      if (nrank==0) print *,'start AB3 with Euler',itime
      !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               !uhat
               td1(i,j,k)= dt*ta1(i,j,k)
               te1(i,j,k)= dt*tb1(i,j,k) 
               tf1(i,j,k)= dt*tc1(i,j,k)
               
               !save (N+Lxz)(n-1)
               gx(i,j,k)=ta1(i,j,k)
               gy(i,j,k)=tb1(i,j,k)
               gz(i,j,k)=tc1(i,j,k)
               
            enddo
         enddo
      enddo

   else if ((itime.eq.2).and.(ilit.eq.0)) then

      if (nrank==0) print *,'then with AB2',itime
      !THEN AB2 + IMPLICIT CN SCHEMES
      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               !uhat
               td1(i,j,k)= 1.5*dt*ta1(i,j,k)-0.5*dt*gx(i,j,k)-px1(i,j,k)
               te1(i,j,k)= 1.5*dt*tb1(i,j,k)-0.5*dt*gy(i,j,k)-py1(i,j,k)
               tf1(i,j,k)= 1.5*dt*tc1(i,j,k)-0.5*dt*gz(i,j,k)-pz1(i,j,k)

               !save (N+Lxz)(n-1)&(n-2)
               hx(i,j,k)=gx(i,j,k)
               hy(i,j,k)=gy(i,j,k)
               hz(i,j,k)=gz(i,j,k)
               gx(i,j,k)=ta1(i,j,k)
               gy(i,j,k)=tb1(i,j,k)
               gz(i,j,k)=tc1(i,j,k)  
            enddo
         enddo
      enddo
         
   else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               !uhat
               td1(i,j,k)= adt(itr)*ta1(i,j,k)+bdt(itr)*gx(i,j,k) &
                    +cdt(itr)*hx(i,j,k)-px1(i,j,k)
               te1(i,j,k)= adt(itr)*tb1(i,j,k)+bdt(itr)*gy(i,j,k) &
                    +cdt(itr)*hy(i,j,k)-py1(i,j,k)
               tf1(i,j,k)= adt(itr)*tc1(i,j,k)+bdt(itr)*gz(i,j,k) &
                    +cdt(itr)*hz(i,j,k)-pz1(i,j,k)

               !save (N+Lxz)(n-1)
               hx(i,j,k)=gx(i,j,k)
               hy(i,j,k)=gy(i,j,k)
               hz(i,j,k)=gz(i,j,k)
               gx(i,j,k)=ta1(i,j,k)
               gy(i,j,k)=tb1(i,j,k)
               gz(i,j,k)=tc1(i,j,k)  
            enddo
         enddo
      enddo

   endif
endif

call transpose_x_to_y(ux1,ux2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)
call transpose_y_to_z(ux2,ux3)
call transpose_y_to_z(uy2,uy3)
call transpose_y_to_z(uz2,uz3)

! DIFFUSIVE TERMS IN Z
iimplicit=0
call derzz (ta3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tb3,uy3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
call derzz (tc3,uz3,di3,sz,sfz ,ssz ,swz ,zsize(1),zsize(2),zsize(3),0)
call transpose_z_to_y(ta3,ta2)
call transpose_z_to_y(tb3,tb2)
call transpose_z_to_y(tc3,tc2)
!DIFFUSIVE TERMS IN Y
call deryy (td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
call deryy (te2,uy2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
call deryy (tf2,uz2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
if (istret.ne.0) then
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      td2(i,j,k)=td2(i,j,k)*pp2y(j)
   enddo
   enddo
   enddo
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      te2(i,j,k)=te2(i,j,k)*pp2y(j)
   enddo
   enddo
   enddo
   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      tf2(i,j,k)=tf2(i,j,k)*pp2y(j)
   enddo
   enddo
   enddo
endif
ta2 = td2 + ta2;
tb2 = te2 + tb2;
tc2 = tf2 + tc2;
call transpose_y_to_x(ta2,ta1)
call transpose_y_to_x(tb2,tb1)
call transpose_y_to_x(tc2,tc1)
!DIFFUSIVE TERMS IN X
call derxx (tg1,ux1,di1,sx,sfx ,ssx ,swx ,xsize(1),xsize(2),xsize(3),0)
call derxx (th1,uy1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
call derxx (ti1,uz1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
iimplicit=1
! Total RHS with diffusion + current velocity
ta1=td1+dt*xnu*(ta1+tg1)+ux1
tb1=te1+dt*xnu*(tb1+th1)+uy1
tc1=tf1+dt*xnu*(tc1+ti1)+uz1

! X pencil inversion
call xmultmatrix7(td1,ta1,ux1,0)
call xmultmatrix7(te1,tb1,uy1,1)
call xmultmatrix7(tf1,tc1,uz1,1)
! rhs complet
ta1(:,:,:)=td1(:,:,:)+ta1(:,:,:)
tb1(:,:,:)=te1(:,:,:)+tb1(:,:,:)
tc1(:,:,:)=tf1(:,:,:)+tc1(:,:,:)
!
ux1=0.; uy1=0.; uz1=0.
if (nclx.eq.0) then
   call xseptinv(ux1,ta1,xggm0,xhhm0,xssm0,xrrm0,xvvm0,xwwm0,xzzm0,xl1m,xl2m,xl3m,xu1m,xu2m,xu3m,xsize(1),xsize(2),xsize(3))
   call xseptinv(uy1,tb1,xggm0,xhhm0,xssm0,xrrm0,xvvm0,xwwm0,xzzm0,xl1m,xl2m,xl3m,xu1m,xu2m,xu3m,xsize(1),xsize(2),xsize(3))
   call xseptinv(uz1,tc1,xggm0,xhhm0,xssm0,xrrm0,xvvm0,xwwm0,xzzm0,xl1m,xl2m,xl3m,xu1m,xu2m,xu3m,xsize(1),xsize(2),xsize(3))
elseif (nclx.eq.1) then
   call septinv(ux1,ta1,xggm10,xhhm10,xssm10,xrrm10,xvvm10,xwwm10,xzzm10,xsize(1),xsize(2),xsize(3))
   call septinv(uy1,tb1,xggm11,xhhm11,xssm11,xrrm11,xvvm11,xwwm11,xzzm11,xsize(1),xsize(2),xsize(3))
   call septinv(uz1,tc1,xggm11,xhhm11,xssm11,xrrm11,xvvm11,xwwm11,xzzm11,xsize(1),xsize(2),xsize(3))
elseif (nclx.eq.2) then

endif

! Z pencil inversion
call transpose_x_to_y(ux1,ta2)
call transpose_x_to_y(uy1,tb2)
call transpose_x_to_y(uz1,tc2)
call transpose_y_to_z(ta2,ta3)
call transpose_y_to_z(tb2,tb3)
call transpose_y_to_z(tc2,tc3)
!
call zmultmatrix7(td3,ta3,ux3,1)
call zmultmatrix7(te3,tb3,uy3,1)
call zmultmatrix7(tf3,tc3,uz3,0)
! rhs complet
ta3(:,:,:)=td3(:,:,:)+ta3(:,:,:)
tb3(:,:,:)=te3(:,:,:)+tb3(:,:,:)
tc3(:,:,:)=tf3(:,:,:)+tc3(:,:,:)
!
ux3=0.; uy3=0.; uz3=0.
if (nclz.eq.0) then
   call zseptinv(ux3,ta3,zggm0,zhhm0,zssm0,zrrm0,zvvm0,zwwm0,zzzm0,zl1m,zl2m,zl3m,zu1m,zu2m,zu3m,zsize(1),zsize(2),zsize(3))
   call zseptinv(uy3,tb3,zggm0,zhhm0,zssm0,zrrm0,zvvm0,zwwm0,zzzm0,zl1m,zl2m,zl3m,zu1m,zu2m,zu3m,zsize(1),zsize(2),zsize(3))
   call zseptinv(uz3,tc3,zggm0,zhhm0,zssm0,zrrm0,zvvm0,zwwm0,zzzm0,zl1m,zl2m,zl3m,zu1m,zu2m,zu3m,zsize(1),zsize(2),zsize(3))
elseif (nclz.eq.1) then
   call septinv(ux3,ta3,zggm11,zhhm11,zssm11,zrrm11,zvvm11,zwwm11,zzzm11,zsize(1),zsize(2),zsize(3))
   call septinv(uy3,tb3,zggm11,zhhm11,zssm11,zrrm11,zvvm11,zwwm11,zzzm11,zsize(1),zsize(2),zsize(3))
   call septinv(uz3,tc3,zggm10,zhhm10,zssm10,zrrm10,zvvm10,zwwm10,zzzm10,zsize(1),zsize(2),zsize(3))
elseif (nclz.eq.2) then

endif

! Y pencil inversion
call transpose_z_to_y(ux3,ta2)
call transpose_z_to_y(uy3,tb2)
call transpose_z_to_y(uz3,tc2)
!
call multmatrix7(td2,ta2,ux2,1)
call multmatrix7(te2,tb2,uy2,0)
call multmatrix7(tf2,tc2,uz2,1)
! rhs complet
ta2(:,:,:)=td2(:,:,:)+ta2(:,:,:)
tb2(:,:,:)=te2(:,:,:)+tb2(:,:,:)
tc2(:,:,:)=tf2(:,:,:)+tc2(:,:,:)
! CL Burgraff
if (itype.eq.7) then
  do k=1,ysize(3)
    z=(k-1+ystart(3)-1)*dz
  j=ysize(2)
    y=(j-1+ystart(2)-1)*dy
  do i=1,ysize(1)
    x=(i-1+ystart(1)-1)*dx
    ta2(i,ysize(2),k)= 8.*(x**4-2.*x**3+x**2)*(4.*y**3-2.*y)
  enddo
  enddo
endif
!
ux2=0.; uy2=0.; uz2=0.;
if (ncly.eq.0) then
   call septinv(ux2,ta2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
   call septinv(uy2,tb2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
   call septinv(uz2,tc2,ggm0,hhm0,ssm0,rrm0,vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ysize(1),ysize(2),ysize(3))
elseif (ncly.eq.1) then
   call septinv(ux2,ta2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
   call septinv(uy2,tb2,ggm10,hhm10,ssm10,rrm10,vvm10,wwm10,zzm10,ysize(1),ysize(2),ysize(3))
   call septinv(uz2,tc2,ggm11,hhm11,ssm11,rrm11,vvm11,wwm11,zzm11,ysize(1),ysize(2),ysize(3))
elseif (ncly.eq.2) then
   call septinv(ux2,ta2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
   call septinv(uy2,tb2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
   call septinv(uz2,tc2,ggm,hhm,ssm,rrm,vvm,wwm,zzm,ysize(1),ysize(2),ysize(3))
endif

call transpose_y_to_x(ux2,ux1)
call transpose_y_to_x(uy2,uy1)
call transpose_y_to_x(uz2,uz1)

! IF ab2 OR ab3
if ((nscheme.eq.1).or.(nscheme.eq.4)) then
   ! IF NOT (first time step with no restart)
   if ( (ilit.eq.1).or.(itime.gt.1) ) then
      ux1(:,:,:)=ux1(:,:,:)+px1(:,:,:)
      uy1(:,:,:)=uy1(:,:,:)+py1(:,:,:)
      uz1(:,:,:)=uz1(:,:,:)+pz1(:,:,:)
   endif
endif

return
end subroutine inttimp

!********************************************************************
!PREMULT FOR INTTIMP WITH SEPTADIAG
subroutine multmatrix7(td2,ta2,ux2,npaire)
! 
!********************************************************************
USE param, only: xcst, ncly, istret
USE variables, only: sy,sfy,sfyp,ssy,ssyp,swy,swyp,pp2y
USE derivY, only: alsajy,alsa2y,alsa3y,alsaty,alsamy
USE decomp_2d, only: mytype, ysize

implicit none

integer, intent(in) :: npaire
integer :: j,jm1,jp1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2

   if (istret.ne.0) then
      do j=1,ysize(2)
         ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
      enddo
   endif

   if (ncly.eq.0) then

      j=1; jm1=ysize(2)
      td2(:,j,:) = alsajy*ta2(:,jm1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      j=ysize(2); jp1=1
      td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,jp1,:)

   elseif (ncly.eq.1) then

      if (npaire.eq.0) then
         td2(:,1,:) = ta2(:,1,:)
      else
         td2(:,1,:) = 2.*alsajy*ta2(:,2,:) + ta2(:,1,:)
      endif
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      if (npaire.eq.0) then
         td2(:,ysize(2),:) = ta2(:,ysize(2),:)
      else
         td2(:,ysize(2),:) = 2.*alsajy*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:)
      endif

   elseif (ncly.eq.2) then

      td2(:,1,:) = 0.
      j=2
      td2(:,j,:) = alsa2y*ta2(:,j-1,:) + ta2(:,j,:) + alsa2y*ta2(:,j+1,:)
      j=3
      td2(:,j,:) = alsa3y*ta2(:,j-1,:) + ta2(:,j,:) + alsa3y*ta2(:,j+1,:)
      do j=4,ysize(2)-3
         td2(:,j,:) = alsajy*ta2(:,j-1,:) + ta2(:,j,:) + alsajy*ta2(:,j+1,:)
      enddo
      j=ysize(2)-2
      td2(:,j,:) = alsaty*ta2(:,j-1,:) + ta2(:,j,:) + alsaty*ta2(:,j+1,:)
      j=ysize(2)-1
      td2(:,j,:) = alsamy*ta2(:,j-1,:) + ta2(:,j,:) + alsamy*ta2(:,j+1,:)
      td2(:,ysize(2),:) = 0.

   endif
   ta2=td2

   if ( (ncly.eq.1).and.(npaire==0) ) then
      call deryy(td2,ux2,di2,sy,sfy,ssy,swy,ysize(1),ysize(2),ysize(3),0)
   else
      call deryy(td2,ux2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)
   endif

   do j=1,ysize(2)
      td2(:,j,:) = - xcst*td2(:,j,:)
   enddo

   if (ncly.eq.2) then
      td2(:,1,:) = 0.
      td2(:,ysize(2),:) = 0.
   endif

end subroutine multmatrix7
!
subroutine xmultmatrix7(td1,ta1,ux1,npaire)
! 
!********************************************************************
USE param, only: xcst, nclx
USE variables, only: sx,sfx,sfxp,ssx,ssxp,swx,swxp
USE derivX, only: alsaix,alsa2x,alsa3x,alsatx,alsamx
USE decomp_2d, only: mytype, xsize

implicit none

integer, intent(in) :: npaire
integer :: i,im1,ip1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: ux1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)), intent(inout) :: td1,ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1

   if (nclx.eq.0) then

      i=1; im1=xsize(1)
      td1(i,:,:) = alsaix*ta1(im1,:,:) + ta1(i,:,:) + alsaix*ta1(i+1,:,:)
      do i=2,xsize(1)-1
         td1(i,:,:) = alsaix*ta1(i-1,:,:) + ta1(i,:,:) + alsaix*ta1(i+1,:,:)
      enddo
      i=xsize(1); ip1=1
      td1(i,:,:) = alsaix*ta1(i-1,:,:) + ta1(i,:,:) + alsaix*ta1(ip1,:,:)

   elseif (nclx.eq.1) then

      if (npaire.eq.0) then
         td1(1,:,:) = ta1(1,:,:)
      else
         td1(1,:,:) = 2.*alsaix*ta1(2,:,:) + ta1(1,:,:)
      endif
      do i=2,xsize(1)-1
         td1(i,:,:) = alsaix*ta1(i-1,:,:) + ta1(i,:,:) + alsaix*ta1(i+1,:,:)
      enddo
      if (npaire.eq.0) then
         td1(xsize(1),:,:) = ta1(xsize(1),:,:)
      else
         td1(xsize(1),:,:) = 2.*alsaix*ta1(xsize(1)-1,:,:) + ta1(xsize(1),:,:)
      endif
   
   elseif (nclx.eq.2) then

      td1(1,:,:) = 0.
      i=2
      td1(i,:,:) = alsa2x*ta1(i-1,:,:) + ta1(i,:,:) + alsa2x*ta1(i+1,:,:)
      i=3
      td1(i,:,:) = alsa3x*ta1(i-1,:,:) + ta1(i,:,:) + alsa3x*ta1(i+1,:,:)
      do i=4,xsize(1)-3
         td1(i,:,:) = alsaix*ta1(i-1,:,:) + ta1(i,:,:) + alsaix*ta1(i+1,:,:)
      enddo
      i=xsize(1)-2
      td1(i,:,:) = alsatx*ta1(i-1,:,:) + ta1(i,:,:) + alsatx*ta1(i+1,:,:)
      i=xsize(1)-1
      td1(i,:,:) = alsamx*ta1(i-1,:,:) + ta1(i,:,:) + alsamx*ta1(i+1,:,:)
      td1(xsize(1),:,:) = 0.

   endif
   ta1=td1

   if ( (nclx.eq.1).and.(npaire==0) ) then
      call derxx(td1,ux1,di1,sx,sfx,ssx,swx,xsize(1),xsize(2),xsize(3),0)
   else
      call derxx(td1,ux1,di1,sx,sfxp,ssxp,swxp,xsize(1),xsize(2),xsize(3),1)
   endif

   do i=1,xsize(1)
      td1(i,:,:) = - xcst*td1(i,:,:)
   enddo

   if (nclx.eq.2) then
      td1(1,:,:) = 0.
      td1(xsize(1),:,:) = 0.
   endif

end subroutine xmultmatrix7
!
subroutine zmultmatrix7(td3,ta3,ux3,npaire)
! 
!********************************************************************
USE param, only: xcst, nclz
USE variables, only: sz,sfz,sfzp,ssz,sszp,swz,swzp
USE derivZ, only: alsakz,alsa2z,alsa3z,alsamz,alsatz
USE decomp_2d, only: mytype, zsize

implicit none

integer, intent(in) :: npaire
integer :: k,km1,kp1
real(mytype),dimension(zsize(1),zsize(2),zsize(3)), intent(inout) :: ux3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)), intent(inout) :: td3,ta3
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: di3

   if (nclz.eq.0) then

      k=1; km1=zsize(3)
      td3(:,:,k) = alsakz*ta3(:,:,km1) + ta3(:,:,k) + alsakz*ta3(:,:,k+1)
      do k=2,zsize(3)-1
         td3(:,:,k) = alsakz*ta3(:,:,k-1) + ta3(:,:,k) + alsakz*ta3(:,:,k+1)
      enddo
      k=zsize(3); kp1=1
      td3(:,:,k) = alsakz*ta3(:,:,k-1) + ta3(:,:,k) + alsakz*ta3(:,:,kp1)

   elseif (nclz.eq.1) then

      if (npaire.eq.0) then
         td3(:,:,1) = ta3(:,:,1)
      else
         td3(:,:,1) = 2.*alsakz*ta3(:,:,2) + ta3(:,:,1)
      endif
      do k=2,zsize(3)-1
         td3(:,:,k) = alsakz*ta3(:,:,k-1) + ta3(:,:,k) + alsakz*ta3(:,:,k+1)
      enddo
      if (npaire.eq.0) then
         td3(:,:,zsize(3)) = ta3(:,:,zsize(3))
      else
         td3(:,:,zsize(3)) = 2.*alsakz*ta3(:,:,zsize(3)-1) + ta3(:,:,zsize(3))
      endif

   elseif (nclz.eq.2) then

      td3(:,:,1) = 0.
      k=2
      td3(:,:,k) = alsa2z*ta3(:,:,k-1) + ta3(:,:,k) + alsa2z*ta3(:,:,k+1)
      k=3
      td3(:,:,k) = alsa3z*ta3(:,:,k-1) + ta3(:,:,k) + alsa3z*ta3(:,:,k+1)      
      do k=4,zsize(3)-3
         td3(:,:,k) = alsakz*ta3(:,:,k-1) + ta3(:,:,k) + alsakz*ta3(:,:,k+1)
      enddo
      k=zsize(3)-2
      td3(:,:,k) = alsatz*ta3(:,:,k-1) + ta3(:,:,k) + alsatz*ta3(:,:,k+1)
      k=zsize(3)-1
      td3(:,:,k) = alsamz*ta3(:,:,k-1) + ta3(:,:,k) + alsamz*ta3(:,:,k+1)
      td3(:,:,zsize(3)) = 0.

   endif
   ta3=td3

   if ( (nclz.eq.1).and.(npaire==0) ) then
      call derzz(td3,ux3,di3,sz,sfz,ssz,swz,zsize(1),zsize(2),zsize(3),0)
   else
      call derzz(td3,ux3,di3,sz,sfzp,sszp,swzp,zsize(1),zsize(2),zsize(3),1)
   endif

   do k=1,zsize(3)
      td3(:,:,k) = - xcst*td3(:,:,k)
   enddo

   if (nclz.eq.2) then
      td3(:,:,1) = 0.
      td3(:,:,zsize(3)) = 0.
   endif

end subroutine zmultmatrix7

!********************************************************************
!
subroutine implicit_schemes()
!
!********************************************************************

USE param
USE derivX
USE derivY
USE derivZ
USE variables
USE var
use ludecomp

implicit none

integer  :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MATRICE M2 TIME IMPLICIT !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!FOR VELOCITY  !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!

xcst = xnu*dt*0.5

!!!!!!!!!!!!!!!!!!!!!!
!CAS SEPTADIAGONALE !!
!!!!!!!!!!!!!!!!!!!!!!

!!! NCL = 2, dirichlet imposé, fonction nulle à la paroi
!
!DIAG
aam(1     )=as1y
aam(ny    )=asny
aam(2     )=-2.*as2y
aam(ny-1  )=-2.*asmy
aam(3     )=-2.*(as3y+bs3y)
aam(ny-2  )=-2.*(asty+bsty)
aam(4:ny-3)=-2.*(asjy+bsjy+csjy)
if (istret==0) then
   aam = 1.-xcst*aam
else
   aam = 1./pp2y -xcst*aam
endif
!CL sur aam
aam(1 )=1.
aam(ny)=1.
!
!DIAG SUP 1
bbm(1     )=bs1y
bbm(ny    )=bsny
bbm(2     )=as2y
bbm(ny-1  )=asmy
bbm(3     )=as3y
bbm(ny-2  )=asty
bbm(4:ny-3)=asjy
bbm = -xcst*bbm
if (istret==0) then
  bbm(2     )=bbm(2     )+alsa2y
  bbm(ny-1  )=bbm(ny-1  )+alsamy
  bbm(3     )=bbm(3     )+alsa3y
  bbm(ny-2  )=bbm(ny-2  )+alsaty
  bbm(4:ny-3)=bbm(4:ny-3)+alsajy
else
  bbm(2     )=bbm(2     )+alsa2y/pp2y(3)
  bbm(ny-1  )=bbm(ny-1  )+alsamy/pp2y(ny)
  bbm(3     )=bbm(3     )+alsa3y/pp2y(4)
  bbm(ny-2  )=bbm(ny-2  )+alsaty/pp2y(ny-1)
  bbm(4:ny-3)=bbm(4:ny-3)+alsajy/pp2y(5:ny-2)
endif
!CL sur bbm
bbm(1 )=0.
bbm(ny)=0.
!
!DIAG SUP 2
ccm(1     )=cs1y
ccm(ny    )=csny
ccm(2     )=0.
ccm(ny-1  )=0.
ccm(3     )=bs3y
ccm(ny-2  )=bsty
ccm(4:ny-3)=bsjy
ccm = -xcst*ccm
!CL sur ccm
ccm(1 )=0.
ccm(ny)=0.
!
!DIAG SUP 3
rrm(1     )=ds1y
rrm(ny    )=dsny
rrm(2     )=0. 
rrm(ny-1  )=0.
rrm(3     )=0.
rrm(ny-2  )=0.
rrm(4:ny-3)=csjy
rrm = -xcst*rrm
!CL sur rrm
rrm(1 )=0.
rrm(ny)=0.
!
!DIAG INF 1
if (istret==0) then
  ddm=bbm
else
  ddm(1     )=bs1y
  ddm(ny    )=bsny
  ddm(2     )=as2y
  ddm(ny-1  )=asmy
  ddm(3     )=as3y
  ddm(ny-2  )=asty
  ddm(4:ny-3)=asjy
  ddm = -xcst*ddm
  ddm(2     )=ddm(2     )+alsa2y/pp2y(1)
  ddm(ny-1  )=ddm(ny-1  )+alsamy/pp2y(ny-2)
  ddm(3     )=ddm(3     )+alsa3y/pp2y(2)
  ddm(ny-2  )=ddm(ny-2  )+alsaty/pp2y(ny-3)
  ddm(4:ny-3)=ddm(4:ny-3)+alsajy/pp2y(3:ny-4)
  !CL sur ddm
  ddm(1 )=0.
  ddm(ny)=0.
endif
!
!DIAG INF 2
eem=ccm
!
!DIAG INF 3
qqm=rrm

!!! NCL = 1, npaire=0, dirichlet imposé, fonction impaire
!
! DIAG
   xaam10(1     )=0.
   xaam10(nx    )=0.
   xaam10(2     )=-2.*asix-3.*bsix-2.*csix
   xaam10(nx-1  )=xaam10(2)
   xaam10(3:nx-2)=-2.*(asix+bsix+csix)
   xaam10 = 1. - xcst*xaam10
aam10(1     )=0.
aam10(ny    )=0.
aam10(2     )=-2.*asjy-3.*bsjy-2.*csjy
aam10(ny-1  )=aam10(2)
aam10(3:ny-2)=-2.*(asjy+bsjy+csjy)
if (istret==0) then
   aam10 = 1. - xcst*aam10
else
   aam10 = 1./pp2y - xcst*aam10
endif
   zaam10(1     )=0.
   zaam10(nz    )=0.
   zaam10(2     )=-2.*askz-3.*bskz-2.*cskz
   zaam10(nz-1  )=zaam10(2)
   zaam10(3:nz-2)=-2.*(askz+bskz+cskz)
   zaam10 = 1. - xcst*zaam10
!
!DIAG SUP 1
   xbbm10(1     )=0.
   xbbm10(nx    )=0.
   xbbm10(2     )=asix-csix
   xbbm10(nx-1  )=asix
   xbbm10(3     )=asix
   xbbm10(nx-2  )=asix-csix
   xbbm10(4:nx-3)=asix
   xbbm10 = alsaix - xcst*xbbm10
   xbbm10(1 )=0.
   xbbm10(nx)=0.
bbm10(1     )=0.
bbm10(ny    )=0.
bbm10(2     )=asjy-csjy
bbm10(ny-1  )=asjy
bbm10(3     )=asjy
bbm10(ny-2  )=asjy-csjy
bbm10(4:ny-3)=asjy
if (istret==0) then
   bbm10 = alsajy - xcst*bbm10
else
   bbm10(2:ny-1) = alsajy/pp2y(3:ny) - xcst*bbm10(2:ny-1)
endif
!CL sur bbm10
bbm10(1 )=0.
bbm10(ny)=0.
   zbbm10(1     )=0.
   zbbm10(nz    )=0.
   zbbm10(2     )=askz-cskz
   zbbm10(nz-1  )=askz
   zbbm10(3     )=askz
   zbbm10(nz-2  )=askz-cskz
   zbbm10(4:nz-3)=askz
   zbbm10 = alsakz - xcst*zbbm10
   zbbm10(1 )=0.
   zbbm10(nz)=0.
!
!DIAG SUP 2
   xccm10(1     )=0.
   xccm10(nx    )=0.
   xccm10(2     )=bsix
   xccm10(nx-1  )=0.
   xccm10(3     )=bsix
   xccm10(nx-2  )=bsix
   xccm10(4:nx-3)=bsix
   xccm10 = -xcst*xccm10
ccm10(1     )=0.
ccm10(ny    )=0.
ccm10(2     )=bsjy
ccm10(ny-1  )=0.
ccm10(3     )=bsjy
ccm10(ny-2  )=bsjy
ccm10(4:ny-3)=bsjy
ccm10 = -xcst*ccm10
   zccm10(1     )=0.
   zccm10(nz    )=0.
   zccm10(2     )=bskz
   zccm10(nz-1  )=0.
   zccm10(3     )=bskz
   zccm10(nz-2  )=bskz
   zccm10(4:nz-3)=bskz
   zccm10 = -xcst*zccm10
!
!DIAG SUP 3
   xrrm10(1     )=0.
   xrrm10(nx    )=0.
   xrrm10(2     )=csix
   xrrm10(nx-1  )=0.
   xrrm10(3     )=csix
   xrrm10(nx-2  )=0.
   xrrm10(4:nx-3)=csix
   xrrm10 = -xcst*xrrm10
rrm10(1     )=0.
rrm10(ny    )=0.
rrm10(2     )=csjy
rrm10(ny-1  )=0.
rrm10(3     )=csjy
rrm10(ny-2  )=0.
rrm10(4:ny-3)=csjy
rrm10 = -xcst*rrm10
   zrrm10(1     )=0.
   zrrm10(nz    )=0.
   zrrm10(2     )=cskz
   zrrm10(nz-1  )=0.
   zrrm10(3     )=cskz
   zrrm10(nz-2  )=0.
   zrrm10(4:nz-3)=cskz
   zrrm10 = -xcst*zrrm10
!
!DIAG INF 1
   xddm10(1     )=0.
   xddm10(nx    )=0.
   xddm10(2     )=asix
   xddm10(nx-1  )=asix-csix
   xddm10(3     )=asix-csix
   xddm10(nx-2  )=asix
   xddm10(4:nx-3)=asix
   xddm10 = alsajy - xcst*xddm10
   xddm10(1 )=0.
   xddm10(nx)=0.
ddm10(1     )=0.
ddm10(ny    )=0.
ddm10(2     )=asjy
ddm10(ny-1  )=asjy-csjy
ddm10(3     )=asjy-csjy
ddm10(ny-2  )=asjy
ddm10(4:ny-3)=asjy
if (istret==0) then
   ddm10 = alsajy - xcst*ddm10
else
   ddm10(2:ny-1) = alsajy/pp2y(1:ny-2) - xcst*ddm10(2:ny-1)
endif
!CL sur ddm10
ddm10(1 )=0.
ddm10(ny)=0.
   zddm10(1     )=0.
   zddm10(nz    )=0.
   zddm10(2     )=askz
   zddm10(nz-1  )=askz-cskz
   zddm10(3     )=askz-cskz
   zddm10(nz-2  )=askz
   zddm10(4:nz-3)=askz
   zddm10 = alsakz - xcst*zddm10
   zddm10(1 )=0.
   zddm10(nz)=0.
!
!DIAG INF 2
   xeem10(1     )=0.
   xeem10(nx    )=0.
   xeem10(2     )=0.
   xeem10(nx-1  )=bsix
   xeem10(3     )=bsix
   xeem10(nx-2  )=bsix
   xeem10(4:nx-3)=bsix
   xeem10 = -xcst*xeem10
eem10(1     )=0.
eem10(ny    )=0.
eem10(2     )=0.
eem10(ny-1  )=bsjy
eem10(3     )=bsjy
eem10(ny-2  )=bsjy
eem10(4:ny-3)=bsjy
eem10 = -xcst*eem10
   zeem10(1     )=0.
   zeem10(nz    )=0.
   zeem10(2     )=0.
   zeem10(nz-1  )=bskz
   zeem10(3     )=bskz
   zeem10(nz-2  )=bskz
   zeem10(4:nz-3)=bskz
   zeem10 = -xcst*zeem10
!
!DIAG INF 3
   xqqm10(1     )=0.
   xqqm10(nx    )=0.
   xqqm10(2     )=0.
   xqqm10(nx-1  )=csix
   xqqm10(3     )=0.
   xqqm10(nx-2  )=csix
   xqqm10(4:nx-3)=csix
   xqqm10 = -xcst*xqqm10
qqm10(1     )=0.
qqm10(ny    )=0.
qqm10(2     )=0.
qqm10(ny-1  )=csjy
qqm10(3     )=0.
qqm10(ny-2  )=csjy
qqm10(4:ny-3)=csjy
qqm10 = -xcst*qqm10
   zqqm10(1     )=0.
   zqqm10(nz    )=0.
   zqqm10(2     )=0.
   zqqm10(nz-1  )=cskz
   zqqm10(3     )=0.
   zqqm10(nz-2  )=cskz
   zqqm10(4:nz-3)=cskz
   zqqm10 = -xcst*zqqm10

!!! NCL = 1, npaire=1, neumann imposé, fonction paire
!
! DIAG
   xaam11(1     )=-2.*(asix+bsix+csix)
   xaam11(nx    )=xaam11(1)
   xaam11(2     )=-2.*asix-bsix-2.*csix
   xaam11(nx-1  )=xaam11(2)
   xaam11(3:nx-2)=-2.*(asix+bsix+csix)
   xaam11 = 1. - xcst*xaam11
aam11(1     )=-2.*(asjy+bsjy+csjy)
aam11(ny    )=aam11(1)
aam11(2     )=-2.*asjy-bsjy-2.*csjy
aam11(ny-1  )=aam11(2)
aam11(3:ny-2)=-2.*(asjy+bsjy+csjy)
if (istret==0) then
   aam11 = 1. - xcst*aam11
else
   aam11 = 1./pp2y - xcst*aam11
endif
   zaam11(1     )=-2.*(askz+bskz+cskz)
   zaam11(nz    )=zaam11(1)
   zaam11(2     )=-2.*askz-bskz-2.*cskz
   zaam11(nz-1  )=zaam11(2)
   zaam11(3:nz-2)=-2.*(askz+bskz+cskz)
   zaam11 = 1. - xcst*zaam11
!
!DIAG SUP 1
   xbbm11(1     )=2.*asix
   xbbm11(nx    )=0.
   xbbm11(2     )=asix+csix
   xbbm11(nx-1  )=asix
   xbbm11(3     )=asix
   xbbm11(nx-2  )=asix+csix
   xbbm11(4:nx-3)=asix
   xbbm11 = alsaix - xcst*xbbm11
   xbbm11(1 )=xbbm11(1)+alsaix
   xbbm11(nx)=0.
bbm11(1     )=2.*asjy
bbm11(ny    )=0.
bbm11(2     )=asjy+csjy
bbm11(ny-1  )=asjy
bbm11(3     )=asjy
bbm11(ny-2  )=asjy+csjy
bbm11(4:ny-3)=asjy
if (istret==0) then
   bbm11 = alsajy - xcst*bbm11
else
   bbm11(1:ny-1) = alsajy/pp2y(2:ny) - xcst*bbm11(1:ny-1)
endif
!CL sur bbm11
if (istret.eq.0) then
   bbm11(1 )=bbm11(1)+alsajy
else
   bbm11(1 )=bbm11(1)+alsajy/pp2y(2)
endif
bbm11(ny)=0.
   zbbm11(1     )=2.*askz
   zbbm11(nz    )=0.
   zbbm11(2     )=askz+cskz
   zbbm11(nz-1  )=askz
   zbbm11(3     )=askz
   zbbm11(nz-2  )=askz+cskz
   zbbm11(4:nz-3)=askz
   zbbm11 = alsakz - xcst*zbbm11
   zbbm11(1 )=zbbm11(1)+alsakz
   zbbm11(nz)=0.
!
!DIAG SUP 2
   xccm11(1     )=2.*bsix
   xccm11(nx    )=0.
   xccm11(2     )=bsix
   xccm11(nx-1  )=0.
   xccm11(3     )=bsix
   xccm11(nx-2  )=bsix
   xccm11(4:nx-3)=bsix
   xccm11 = -xcst*xccm11
ccm11(1     )=2.*bsjy
ccm11(ny    )=0.
ccm11(2     )=bsjy
ccm11(ny-1  )=0.
ccm11(3     )=bsjy
ccm11(ny-2  )=bsjy
ccm11(4:ny-3)=bsjy
ccm11 = -xcst*ccm11
   zccm11(1     )=2.*bskz
   zccm11(nz    )=0.
   zccm11(2     )=bskz
   zccm11(nz-1  )=0.
   zccm11(3     )=bskz
   zccm11(nz-2  )=bskz
   zccm11(4:nz-3)=bskz
   zccm11 = -xcst*zccm11
!
!DIAG SUP 3
   xrrm11(1     )=2.*csix
   xrrm11(nx    )=0.
   xrrm11(2     )=csix
   xrrm11(nx-1  )=0.
   xrrm11(3     )=csix
   xrrm11(nx-2  )=0.
   xrrm11(4:nx-3)=csix
   xrrm11 = -xcst*xrrm11
rrm11(1     )=2.*csjy
rrm11(ny    )=0.
rrm11(2     )=csjy
rrm11(ny-1  )=0.
rrm11(3     )=csjy
rrm11(ny-2  )=0.
rrm11(4:ny-3)=csjy
rrm11 = -xcst*rrm11
   zrrm11(1     )=2.*cskz
   zrrm11(nz    )=0.
   zrrm11(2     )=cskz
   zrrm11(nz-1  )=0.
   zrrm11(3     )=cskz
   zrrm11(nz-2  )=0.
   zrrm11(4:nz-3)=cskz
   zrrm11 = -xcst*zrrm11
!
!DIAG INF 1
   xddm11(1     )=0.
   xddm11(nx    )=2.*asix
   xddm11(2     )=asix
   xddm11(nx-1  )=asix+csix
   xddm11(3     )=asix+csix
   xddm11(nx-2  )=asix
   xddm11(4:nx-3)=asix
   xddm11 = alsaix - xcst*xddm11
   xddm11(1 )=0.
   xddm11(nx)=xddm11(nx)+alsaix
ddm11(1     )=0.
ddm11(ny    )=2.*asjy
ddm11(2     )=asjy
ddm11(ny-1  )=asjy+csjy
ddm11(3     )=asjy+csjy
ddm11(ny-2  )=asjy
ddm11(4:ny-3)=asjy
if (istret==0) then
   ddm11 = alsajy - xcst*ddm11
else
   ddm11(2:ny) = alsajy/pp2y(1:ny-1) - xcst*ddm11(2:ny)
endif
!CL sur ddm11
ddm11(1 )=0.
if (istret.eq.0) then
   ddm11(ny)=ddm11(ny)+alsajy
else
   ddm11(ny)=ddm11(ny)+alsajy/pp2y(ny-1)!a1
endif
   zddm11(1     )=0.
   zddm11(nz    )=2.*askz
   zddm11(2     )=askz
   zddm11(nz-1  )=askz+cskz
   zddm11(3     )=askz+cskz
   zddm11(nz-2  )=askz
   zddm11(4:nz-3)=askz
   zddm11 = alsakz - xcst*zddm11
   zddm11(1 )=0.
   zddm11(nz)=zddm11(nz)+alsakz
!
!DIAG INF 2
   xeem11(1     )=0.
   xeem11(nx    )=2.*bsix
   xeem11(2     )=0.
   xeem11(nx-1  )=bsix
   xeem11(3     )=bsix
   xeem11(nx-2  )=bsix
   xeem11(4:nx-3)=bsix
   xeem11 = -xcst*xeem11
eem11(1     )=0.
eem11(ny    )=2.*bsjy
eem11(2     )=0.
eem11(ny-1  )=bsjy
eem11(3     )=bsjy
eem11(ny-2  )=bsjy
eem11(4:ny-3)=bsjy
eem11 = -xcst*eem11
   zeem11(1     )=0.
   zeem11(nz    )=2.*bskz
   zeem11(2     )=0.
   zeem11(nz-1  )=bskz
   zeem11(3     )=bskz
   zeem11(nz-2  )=bskz
   zeem11(4:nz-3)=bskz
   zeem11 = -xcst*zeem11
!
!DIAG INF 3
   xqqm11(1     )=0.
   xqqm11(nx    )=2.*csix
   xqqm11(2     )=0.
   xqqm11(nx-1  )=csix
   xqqm11(3     )=0.
   xqqm11(nx-2  )=csix
   xqqm11(4:nx-3)=csix
   xqqm11 = -xcst*xqqm11
qqm11(1     )=0.
qqm11(ny    )=2.*csjy
qqm11(2     )=0.
qqm11(ny-1  )=csjy
qqm11(3     )=0.
qqm11(ny-2  )=csjy
qqm11(4:ny-3)=csjy
qqm11 = -xcst*qqm11
   zqqm11(1     )=0.
   zqqm11(nz    )=2.*cskz
   zqqm11(2     )=0.
   zqqm11(nz-1  )=cskz
   zqqm11(3     )=0.
   zqqm11(nz-2  )=cskz
   zqqm11(4:nz-3)=cskz
   zqqm11 = -xcst*zqqm11

!!! NXL = 0
!DIAG
   xaam0= 1.-xcst*(-2.*(asix+bsix+csix))
if (istret==0) then
   aam0 = 1.-xcst*(-2.*(asjy+bsjy+csjy))
else
   aam0 = 1./pp2y-xcst*(-2.*(asjy+bsjy+csjy))
endif
   zaam0= 1.-xcst*(-2.*(askz+bskz+cskz))
!
!DIAG SUP 1
   xbbm0 = alsaix-xcst*asix
if (istret==0) then
   bbm0 = alsajy-xcst*asjy
else
   bbm0(1:ny-1) = alsajy/pp2y(2:ny) -xcst*asjy
   bbm0(ny) = alsajy/pp2y(1) -xcst*asjy
endif
   zbbm0 = alsakz-xcst*askz
!
!DIAG SUP 2
   xccm0= -xcst*bsix
ccm0 = -xcst*bsjy
   zccm0= -xcst*bskz
!
!DIAG SUP 3
   xrrm0= -xcst*csix
rrm0 = -xcst*csjy
   zrrm0= -xcst*cskz
!
!DIAG INF 1
   xddm0=xbbm0
if (istret==0) then
  ddm0=bbm0
else
   ddm0(1) = alsajy/pp2y(ny) -xcst*asjy
   ddm0(2:ny) = alsajy/pp2y(1:ny-1) -xcst*asjy
endif
   zddm0=zbbm0
!
!DIAG INF 2
   xeem0=xccm0
eem0=ccm0
   zeem0=zccm0
!
!DIAG INF 3
   xqqm0=xrrm0
qqm0=rrm0
   zqqm0=zrrm0

call ludecomp7(aam,bbm,ccm,ddm,eem,qqm,ggm,hhm,ssm,rrm,&
     vvm,wwm,zzm,ny)
call ludecomp7(xaam10,xbbm10,xccm10,xddm10,xeem10,xqqm10,xggm10,xhhm10,xssm10,xrrm10,&
     xvvm10,xwwm10,xzzm10,nx)
call ludecomp7(aam10,bbm10,ccm10,ddm10,eem10,qqm10,ggm10,hhm10,ssm10,rrm10,&
     vvm10,wwm10,zzm10,ny)
call ludecomp7(zaam10,zbbm10,zccm10,zddm10,zeem10,zqqm10,zggm10,zhhm10,zssm10,zrrm10,&
     zvvm10,zwwm10,zzzm10,nz)
call ludecomp7(xaam11,xbbm11,xccm11,xddm11,xeem11,xqqm11,xggm11,xhhm11,xssm11,xrrm11,&
     xvvm11,xwwm11,xzzm11,nx)
call ludecomp7(aam11,bbm11,ccm11,ddm11,eem11,qqm11,ggm11,hhm11,ssm11,rrm11,&
     vvm11,wwm11,zzm11,ny)
call ludecomp7(zaam11,zbbm11,zccm11,zddm11,zeem11,zqqm11,zggm11,zhhm11,zssm11,zrrm11,&
     zvvm11,zwwm11,zzzm11,nz)
call ludecomp7(xaam0,xbbm0,xccm0,xddm0,xeem0,xqqm0,xggm0,xhhm0,xssm0,xrrm0,&
     xvvm0,xwwm0,xzzm0,xl1m,xl2m,xl3m,xu1m,xu2m,xu3m,nx)
call ludecomp7(aam0,bbm0,ccm0,ddm0,eem0,qqm0,ggm0,hhm0,ssm0,rrm0,&
     vvm0,wwm0,zzm0,l1m,l2m,l3m,u1m,u2m,u3m,ny)
call ludecomp7(zaam0,zbbm0,zccm0,zddm0,zeem0,zqqm0,zggm0,zhhm0,zssm0,zrrm0,&
     zvvm0,zwwm0,zzzm0,zl1m,zl2m,zl3m,zu1m,zu2m,zu3m,nz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FIN MATRICE M2 TIME IMPLICIT    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine implicit_schemes

!************************************************************
!
subroutine scalarimp(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3)
!
!************************************************************
!IMPLICIT TIME INTEGRATION FOR D2/DY2
!
!************************************************************
USE param
USE variables
USE decomp_2d
USE derivY
USE MPI
use matinv
#ifdef my_mod_solide
use conjugate_ht, only : temp_bot,temp_top,ny_sol_bot,cl_bot,cl_top,update_temp_solide
#endif

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi1,phis1,phiss1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: di1,tg1,th1,ti1,td1
!real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phisauv

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: uy2,uz2,phi2,ta2,tb2,tc2,td2,di2

real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: uz3,phi3,di3,ta3,tb3
integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nxyz,code

!parametres CL
real(mytype) :: x,y,z,r,lambda,phislbda,adiab,tjet,liss

#ifdef my_mod_solide
real(mytype),dimension(ysize(1),2,ysize(3)) :: mytmptemp
#endif

tg1=0.;th1=0.;ti1=0.;td1=0.
ta2=0.;tb2=0.;tc2=0.;td2=0.
ta3=0.;tb3=0.

nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)
nvect3=zsize(1)*zsize(2)*zsize(3)

if ((nscheme.eq.2).or.(nscheme.eq.3)) then
   if (nrank.eq.0) print *,'RK scheme not ready, abort computation'
   call decomp_2d_finalize
   CALL MPI_FINALIZE(code)
   call exit(0)
endif

!X PENCILS
tg1=ux1*phi1
call derx (th1,tg1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! npaire=0 pour phi
!call derx (th1,phi1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1) ! npaire=1
!call derx (th1,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) ! npaire=0
!call derxxt (tg1,phi1,di1,sx,sfxpt,ssxpt,swxpt,xsize(1),xsize(2),xsize(3),1) ! npaire=1
call derxxt (tg1,phi1,di1,sx,sfxt,ssxt,swxt,xsize(1),xsize(2),xsize(3),0) ! npaire=0
call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(uy1,uy2)
call transpose_x_to_y(uz1,uz2)

!Y PENCILS
tc2=phi2*uy2
call dery (tb2,tc2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=0 pour phi
!call dery (tb2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=1
!call dery (tb2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0
 
if (istret.ne.0) then 
  
!   call dery (tc2,phi2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) ! npaire=1
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0

   do k=1,ysize(3)
   do j=1,ysize(2)
   do i=1,ysize(1)
      ta2(i,j,k)=-pp4y(j)*tc2(i,j,k)
   enddo
   enddo
   enddo

endif

call transpose_y_to_z(phi2,phi3)
call transpose_y_to_z(uz2,uz3)

!Z PENCILS
ta3=uz3*phi3
call derz (tb3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! npaire=0 pour phi
!call derz (tb3,phi3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) ! npaire=1
!call derz (tb3,phi3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0) ! npaire=0
!call derzzt (ta3,phi3,di3,sz,sfzpt,sszpt,swzpt,zsize(1),zsize(2),zsize(3),1) ! npaire=1
call derzzt (ta3,phi3,di3,sz,sfzt,sszt,swzt,zsize(1),zsize(2),zsize(3),0) ! npaire=0

call transpose_z_to_y(ta3,tc2)
call transpose_z_to_y(tb3,td2)

!Y PENCILS ADD TERMS
do k=1,ysize(3)
   do j=1,ysize(2)
      do i=1,ysize(1)
         tc2(i,j,k)=tc2(i,j,k)+ta2(i,j,k) !SECOND DERIVATIVE
!         td2(i,j,k)=uz2(i,j,k)*td2(i,j,k)+uy2(i,j,k)*tb2(i,j,k) !FIRST DERIVATIVE
         td2(i,j,k)=td2(i,j,k)+tb2(i,j,k) !FIRST DERIVATIVE  
      enddo                                                        
   enddo
enddo

call transpose_y_to_x(tc2,ti1)
call transpose_y_to_x(td2,td1) 


!X PENCILS ADD TERMS
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         tg1(i,j,k)= tg1(i,j,k)+ti1(i,j,k) !SECOND DERIVATIVE
!         th1(i,j,k)= ux1(i,j,k)*th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
         th1(i,j,k)= th1(i,j,k)+td1(i,j,k) !FIRST DERIVATIVE
      enddo
   enddo
enddo

do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         tg1(i,j,k)= (xnu/sc)*tg1(i,j,k)-th1(i,j,k) 
      enddo
   enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Kasagi Source term !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!
tg1=tg1+ux1*xnu/sc

! TIME ADVANCEMENT EXPLICIT AB + IMPLICIT CN2 (d2/dy2)
nxyz=xsize(1)*xsize(2)*xsize(3)  

!!!!!!!!!!!!!
!! CN2+AB2 !!
!!!!!!!!!!!!!
if (nscheme.eq.1) then
   if (itime.eq.1.and.ilit.eq.0) then !start with EXPLICIT EULER + CN2

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               td1(i,j,k)= gdt(itr)*tg1(i,j,k)
               phis1(i,j,k)= tg1(i,j,k)
            enddo
         enddo
      enddo
   
   else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               td1(i,j,k)= adt(itr)*tg1(i,j,k)+bdt(itr)*phis1(i,j,k)
               phis1(i,j,k)= tg1(i,j,k)          
            enddo
         enddo
      enddo

   endif

endif


!!!!!!!!!!!!!!!!!!!!
!CN2+AB3         !!!
!!!!!!!!!!!!!!!!!!!!
if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'Temperature start with Euler',itime
      !START WITH EXPLICIT EULER + IMPLICIT CN SCHEMES

      do k=1,xsize(3)
         do j=1,xsize(2)
            do i=1,xsize(1)
               td1(i,j,k)= dt*tg1(i,j,k)
               phis1(i,j,k)= tg1(i,j,k)          
            enddo
         enddo
      enddo
   
   else !CONTINUE WITH EXPLICIT AB2 + IMPLICIT CN SCHEMES
      if  ((itime.eq.2).and.(ilit.eq.0)) then
         if (nrank==0) print *,'then with AB2',itime


         do k=1,xsize(3)
            do j=1,xsize(2)
               do i=1,xsize(1)
                  td1(i,j,k)= 1.5*dt*tg1(i,j,k)-0.5*dt*phis1(i,j,k)
                  phiss1(i,j,k)=phis1(i,j,k)
                  phis1(i,j,k)= tg1(i,j,k)          
               enddo
            enddo
         enddo

      else !FINALLY EXPLICIT AB3 + IMPLICIT CN SCHEMES

         do k=1,xsize(3)
            do j=1,xsize(2)
               do i=1,xsize(1)
                  td1(i,j,k)= adt(itr)*tg1(i,j,k)+bdt(itr)*phis1(i,j,k)&
                       +cdt(itr)*phiss1(i,j,k)
                  phiss1(i,j,k)=phis1(i,j,k)
                  phis1(i,j,k)= tg1(i,j,k)          
               enddo
            enddo
         enddo

      endif
   endif
endif

!Y-PENCIL FOR MATRIX INVERSION
call transpose_x_to_y(phi1,phi2)
call transpose_x_to_y(td1,ta2)

#ifdef my_mod_solide
mytmptemp(:,1,:)=phi2(:,1,:)
mytmptemp(:,2,:)=phi2(:,ysize(2),:)
#endif

!ta2: A.T_hat
!td2:(A+xcstB).Tn
call multmatrix7T(td2,ta2,phi2)

!right hand side
ta2(:,:,:) = ta2(:,:,:) + td2(:,:,:)
if (ncly==2) then
#ifdef my_mod_solide
   ta2(:,1,:) = 0.5*(mytmptemp(:,1,:)+temp_bot(:,ny_sol_bot+3,:))
   ta2(:,ysize(2),:) = 0.5*(mytmptemp(:,2,:)+temp_top(:,1,:))
#else
   ta2(:,1       ,:)=g_0
   ta2(:,ysize(2),:)=g_n
#endif
endif

!Inversion systeme lineaire Mx=b: (A-xcst.B)u^n+1=uhat+(A+xcst.B)u^n
phi2=0.
if (ncly==0) then
   call septinv(phi2,ta2,ggmt0,hhmt0,ssmt0,rrmt0,vvmt0,wwmt0,zzmt0,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt,ysize(1),ysize(2),ysize(3))
elseif (ncly==1) then
  call septinv(phi2,ta2,ggmt1,hhmt1,ssmt1,rrmt1,vvmt1,wwmt1,zzmt1,ysize(1),ysize(2),ysize(3))
elseif (ncly==2) then
  call septinv(phi2,ta2,ggmt,hhmt,ssmt,rrmt,vvmt,wwmt,zzmt,ysize(1),ysize(2),ysize(3))
endif

#ifdef my_mod_solide
   call dery (tc2,phi2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0) ! npaire=0
   cl_bot(:,1,:)=tc2(:,1       ,:)
   cl_top(:,1,:)=tc2(:,ysize(2),:)
   call update_temp_solide()
#endif

call transpose_y_to_x(phi2,phi1)

end subroutine scalarimp

!********************************************************************
!PREMULT FOR INTTIMP WITH SEPTADIAG
subroutine multmatrix7T(td2,ta2,ux2)
! 
!********************************************************************
USE param
USE variables
USE derivY
USE decomp_2d

implicit none

integer :: i,j,k,code
!real(mytype) :: xcst modification module param
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: ux2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)), intent(inout) :: td2,ta2
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: di2

!xcst= xnu*gdt(itr)*0.5

!A.uhat
   if (istret.ne.0) then
      do j=1,ysize(2)
         ta2(:,j,:)=ta2(:,j,:)/pp2y(j)
      enddo
   endif

   if (ncly.eq.0) then

      td2(:,1,:) = alsajyt*ta2(:,2,:) + ta2(:,1,:) + alsajyt*ta2(:,ysize(2),:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
      enddo
      td2(:,ysize(2),:) = alsajyt*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) + alsajyt*ta2(:,1,:)
      ta2=td2

   elseif (ncly.eq.1) then

!      td2(:,1,:) = 2.*alsajyt*ta2(:,2,:) + ta2(:,1,:) ! npaire=1
      td2(:,1,:) = ta2(:,1,:) ! npaire=0
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
      enddo
!      td2(:,ysize(2),:) = 2.*alsajyt*ta2(:,ysize(2)-1,:) + ta2(:,ysize(2),:) ! npaire=1
      td2(:,ysize(2),:) = ta2(:,ysize(2),:) ! npaire=0
      ta2=td2

   elseif (ncly.eq.2) then

      td2(:,1,:) = 0.
      td2(:,2,:) = alsa2y*ta2(:,1,:) + ta2(:,2,:) + alsa2y*ta2(:,3,:)
      td2(:,3,:) = alsa3y*ta2(:,2,:) + ta2(:,3,:) + alsa3y*ta2(:,4,:)
      do j=4,ysize(2)-3
         td2(:,j,:) = alsajyt*ta2(:,j-1,:) + ta2(:,j,:) + alsajyt*ta2(:,j+1,:)
      enddo
      td2(:,ysize(2)-2,:) = alsaty*ta2(:,ysize(2)-3,:) + ta2(:,ysize(2)-2,:) + alsaty*ta2(:,ysize(2)-1,:)
      td2(:,ysize(2)-1,:) = alsamy*ta2(:,ysize(2)-2,:) + ta2(:,ysize(2)-1,:) + alsamy*ta2(:,ysize(2),:)
      td2(:,ysize(2),:) = 0.
      ta2=td2

   endif


!(A+nu*dt.B).un
   if (ncly.eq.0) then

      call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)

   elseif (ncly.eq.1) then

!      call deryyt(td2,ux2,di2,sy,sfypt,ssypt,swypt,ysize(1),ysize(2),ysize(3),1) ! npaire=1
      call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0) ! npaire=0

   elseif (ncly.eq.2) then

      call deryyt(td2,ux2,di2,sy,sfyt,ssyt,swyt,ysize(1),ysize(2),ysize(3),0)

   endif

   if (istret.ne.0) then
      do j=1,ysize(2)
         ux2(:,j,:)=ux2(:,j,:)/pp2y(j)
      enddo
   endif

   if (ncly.eq.0) then

      td2(:,1,:) = alsajyt*ux2(:,2,:) + ux2(:,1,:) + alsajyt*ux2(:,ysize(2),:) &
                    + xcst_pr*td2(:,1,:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
      enddo
      td2(:,ysize(2),:) = alsajyt*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) + alsajyt*ux2(:,1,:) &
                           + xcst_pr*td2(:,ysize(2),:)

   elseif (ncly.eq.1) then

!      td2(:,1,:) = 2.*alsajyt*ux2(:,2,:) + ux2(:,1,:) & ! npaire=1
      td2(:,1,:) = ux2(:,1,:) & ! npaire=0
                 + xcst_pr*td2(:,1,:)
      do j=2,ysize(2)-1
         td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
      enddo
!      td2(:,ysize(2),:) = 2.*alsajyt*ux2(:,ysize(2)-1,:) + ux2(:,ysize(2),:) & ! npaire=1
      td2(:,ysize(2),:) = ux2(:,ysize(2),:) & ! npaire=0
                        + xcst_pr*td2(:,ysize(2),:)

   elseif (ncly.eq.2) then

      td2(:,1,:) = 0.
      td2(:,2,:) = alsa2y*ux2(:,1,:) + ux2(:,2,:) + alsa2y*ux2(:,3,:) &
                 + xcst_pr*td2(:,2,:)
      td2(:,3,:) = alsa3y*ux2(:,2,:) + ux2(:,3,:) + alsa3y*ux2(:,4,:) &
                 + xcst_pr*td2(:,3,:)
      do j=4,ysize(2)-3
         td2(:,j,:) = alsajyt*ux2(:,j-1,:) + ux2(:,j,:) + alsajyt*ux2(:,j+1,:) &
                    + xcst_pr*td2(:,j,:)
      enddo
      td2(:,ysize(2)-2,:) = alsaty*ux2(:,ysize(2)-3,:) + ux2(:,ysize(2)-2,:) + alsaty*ux2(:,ysize(2)-1,:) &
                          + xcst_pr*td2(:,ysize(2)-2,:)
      td2(:,ysize(2)-1,:) = alsamy*ux2(:,ysize(2)-2,:) + ux2(:,ysize(2)-1,:) + alsamy*ux2(:,ysize(2),:) &
                          + xcst_pr*td2(:,ysize(2)-1,:)
      td2(:,ysize(2),:) = 0.

   endif

end subroutine multmatrix7T

!********************************************************************
!
subroutine derxxt(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivX 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
real(mytype), dimension(ny,nz) :: sx
real(mytype),  dimension(nx):: sfx,ssx,swx 

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=asixt*(ux(2,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx  ,j,k))&
           +bsixt*(ux(3,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-1,j,k))&
           +csixt*(ux(4,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-2,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=asixt*(ux(3,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(1   ,j,k))&
           +bsixt*(ux(4,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx  ,j,k))&
           +csixt*(ux(5,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx-1,j,k))
      rx(2,j,k)=0.
      tx(3,j,k)=asixt*(ux(4,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(2 ,j,k))&
           +bsixt*(ux(5,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(1 ,j,k))&
           +csixt*(ux(6,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(nx,j,k))
      rx(3,j,k)=0.
      do i=4,nx-3
         tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-1,j,k))&
              +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-2,j,k))&
              +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-3,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-2,j,k)=asixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-3,j,k))&
           +bsixt*(ux(nx  ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-4,j,k))&
           +csixt*(ux(1   ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-5,j,k))
      rx(nx-2,j,k)=0.
      tx(nx-1,j,k)=asixt*(ux(nx  ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-2,j,k))&
           +bsixt*(ux(1   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-3,j,k))&
           +csixt*(ux(2   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-4,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=asixt*(ux(1 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-1,j,k))&
           +bsixt*(ux(2 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-2,j,k))&
           +csixt*(ux(3 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-3,j,k))
      rx(nx  ,j,k)=alsaixt
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
      sx(j,k)=(   tx(1,j,k)-alsaixt*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alsaixt*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo
endif

if (nclx==1) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=asixt*(ux(2,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(2,j,k))&
              +bsixt*(ux(3,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(3,j,k))&
              +csixt*(ux(4,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(4,j,k))
         tx(2,j,k)=asixt*(ux(3,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(1,j,k))&
              +bsixt*(ux(4,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(2,j,k))&
              +csixt*(ux(5,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(3,j,k))
         tx(3,j,k)=asixt*(ux(4,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))&
              +bsixt*(ux(5,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(1,j,k))&
              +csixt*(ux(6,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-1,j,k))&
                 +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-2,j,k))&
                 +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         tx(nx-2,j,k)=asixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-3,j,k))&
              +bsixt*(ux(nx  ,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-4,j,k))&
              +csixt*(ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-5,j,k))
         tx(nx-1,j,k)=asixt*(ux(nx  ,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-2,j,k))&
              +bsixt*(ux(nx-1,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-3,j,k))&
              +csixt*(ux(nx-2,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-4,j,k))
         tx(nx  ,j,k)=asixt*(ux(nx-1,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-1,j,k))&
              +bsixt*(ux(nx-2,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-2,j,k))&
              +csixt*(ux(nx-3,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-3,j,k))
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=0.
         tx(2,j,k)=asixt*(ux(3,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(1,j,k))&
              +bsixt*(ux(4,j,k)-ux(2,j,k)&
              -ux(2,j,k)-ux(2,j,k))&
              +csixt*(ux(5,j,k)-ux(2,j,k)&
              -ux(2,j,k)-ux(3,j,k))
         tx(3,j,k)=asixt*(ux(4,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))&
              +bsixt*(ux(5,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(1,j,k))&
              +csixt*(ux(6,j,k)-ux(3,j,k)&
              -ux(3,j,k)-ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-1,j,k))&
                 +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-2,j,k))&
                 +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         tx(nx-2,j,k)=asixt*( ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-3,j,k))&
              +bsixt*( ux(nx  ,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-4,j,k))&
              +csixt*(-ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-5,j,k))
         tx(nx-1,j,k)=asixt*( ux(nx  ,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-2,j,k))&
              +bsixt*(-ux(nx-1,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-3,j,k))&
              +csixt*(-ux(nx-2,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-4,j,k))
         tx(nx  ,j,k)=0.
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
      enddo
      enddo
   endif
endif
 
if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
           +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
      tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
           -ux(2,j,k)+ux(1,j,k))
      tx(3,j,k)=as3x*(ux(4,j,k)-ux(3,j,k)&
           -ux(3,j,k)+ux(2,j,k))&
           +bs3x*(ux(5,j,k)-ux(3,j,k)&
           -ux(3,j,k)+ux(1,j,k))
      do i=4,nx-3
         tx(i,j,k)=asixt*(ux(i+1,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-1,j,k))&
              +bsixt*(ux(i+2,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-2,j,k))&
              +csixt*(ux(i+3,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-3,j,k))
      enddo
      tx(nx-2,j,k)=astx*(ux(nx-1,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-3,j,k))&
           +bstx*(ux(nx  ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-4,j,k))
      tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-2,j,k))
      tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
           +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
   enddo
   enddo
endif

return  
end subroutine derxxt

!********************************************************************
!
subroutine deryyt(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivY 

implicit none

integer :: nx,ny,nz,npaire,i,j,k 
real(mytype), dimension(nx,ny,nz) :: ty,uy,ry
real(mytype), dimension(nx,nz) :: sy
real(mytype), dimension(ny) :: sfy,ssy,swy

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=asjyt*(uy(i,2,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny,k))&
           +bsjyt*(uy(i,3,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny-1,k))&
           +csjyt*(uy(i,4,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny-2,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,1,k))&
           +bsjyt*(uy(i,4,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,ny,k))&
           +csjyt*(uy(i,5,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,ny-1,k))
      ry(i,2,k)=0.
      ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,2,k))&
           +bsjyt*(uy(i,5,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,1,k))&
           +csjyt*(uy(i,6,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,ny,k))
      ry(i,3,k)=0.
   enddo
   enddo
   do k=1,nz
   do j=4,ny-3
   do i=1,nx
      ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-1,k))&
           +bsjyt*(uy(i,j+2,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-2,k))&
           +csjyt*(uy(i,j+3,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-3,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-2,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-3,k))&
           +bsjyt*(uy(i,ny  ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-4,k))&
           +csjyt*(uy(i,1   ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-5,k))
      ry(i,ny-2,k)=0.
      ty(i,ny-1,k)=asjyt*(uy(i,ny  ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-2,k))&
           +bsjyt*(uy(i,1   ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-3,k))&
           +csjyt*(uy(i,2   ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-4,k))
      ry(i,ny-1,k)=0.
      ty(i,ny  ,k)=asjyt*(uy(i,1 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-1,k))&
           +bsjyt*(uy(i,2 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-2,k))&
           +csjyt*(uy(i,3 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-3,k))
      ry(i,ny  ,k)=alsajyt
   enddo
   enddo
if (iimplicit.eq.1) return
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*ssy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      ry(i,ny,k)=ry(i,ny,k)*swy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      ry(i,j,k)=(ry(i,j,k)-sfy(j)*ry(i,j+1,k))*swy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(   ty(i,1,k)-alsajyt*ty(i,ny,k))/&
           (1.+ry(i,1,k)-alsajyt*ry(i,ny,k))
   enddo
   enddo   
   do k=1,nz
   do j=1,ny     
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-sy(i,k)*ry(i,j,k)
   enddo
   enddo
   enddo
endif

if (ncly==1) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=asjyt*(uy(i,2,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,2,k))&
              +bsjyt*(uy(i,3,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,3,k))&
              +csjyt*(uy(i,4,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,4,k))
         ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,1,k))&
              +bsjyt*(uy(i,4,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,2,k))&
              +csjyt*(uy(i,5,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,3,k))
         ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))&
              +bsjyt*(uy(i,5,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,1,k))&
              +csjyt*(uy(i,6,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=4,ny-3         
      do i=1,nx   
         ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-1,k))&
              +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-2,k))&
              +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-2,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-3,k))&
              +bsjyt*(uy(i,ny  ,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-4,k))&
              +csjyt*(uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-5,k))
         ty(i,ny-1,k)=asjyt*(uy(i,ny  ,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-2,k))&
              +bsjyt*(uy(i,ny-1,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-3,k))&
              +csjyt*(uy(i,ny-2,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-4,k))
         ty(i,ny  ,k)=asjyt*(uy(i,ny-1,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-1,k))&
              +bsjyt*(uy(i,ny-2,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-2,k))&
              +csjyt*(uy(i,ny-3,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-3,k))
      enddo
      enddo
if (iimplicit.eq.1) return
      do k=1,nz
      do j=2,ny         
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      enddo
      enddo   
      do k=1,nz
      do j=ny-1,1,-1         
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=0.
         ty(i,2,k)=asjyt*(uy(i,3,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,1,k))&
              +bsjyt*(uy(i,4,k)-uy(i,2,k)&
              -uy(i,2,k)-uy(i,2,k))&
              +csjyt*(uy(i,5,k)-uy(i,2,k)&
              -uy(i,2,k)-uy(i,3,k))
         ty(i,3,k)=asjyt*(uy(i,4,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))&
              +bsjyt*(uy(i,5,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,1,k))&
              +csjyt*(uy(i,6,k)-uy(i,3,k)&
              -uy(i,3,k)-uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=4,ny-3
      do i=1,nx   
         ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-1,k))&
              +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-2,k))&
              +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-2,k)=asjyt*( uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-3,k))&
              +bsjyt*( uy(i,ny  ,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-4,k))&
              +csjyt*(-uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-5,k))
         ty(i,ny-1,k)=asjyt*( uy(i,ny  ,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-2,k))&
              +bsjyt*(-uy(i,ny-1,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-3,k))&
              +csjyt*(-uy(i,ny-2,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-4,k))
         ty(i,ny  ,k)=0.
      enddo
      enddo   
if (iimplicit.eq.1) return
      do k=1,nz
      do j=2,ny
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      enddo
      enddo   
      do k=1,nz
      do j=ny-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      enddo
      enddo
      enddo
   endif
endif

if (ncly==2) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=as1y*uy(i,1,k)+bs1y*uy(i,2,k)&
           +cs1y*uy(i,3,k)+ds1y*uy(i,4,k)
      ty(i,2,k)=as2y*(uy(i,3,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,1,k))
      ty(i,3,k)=as3y*(uy(i,4,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,2,k))&
           +bs3y*(uy(i,5,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,1,k))
   enddo
   enddo
   do k=1,nz
   do j=4,ny-3
   do i=1,nx   
      ty(i,j,k)=asjyt*(uy(i,j+1,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-1,k))&
           +bsjyt*(uy(i,j+2,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-2,k))&
           +csjyt*(uy(i,j+3,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-3,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-2,k)=asty*(uy(i,ny-1,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-3,k))&
           +bsty*(uy(i,ny  ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-4,k))
      ty(i,ny-1,k)=asmy*(uy(i,ny  ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-2,k))
      ty(i,ny  ,k)=asny*uy(i,ny  ,k)+bsny*uy(i,ny-1,k)&
           +csny*uy(i,ny-2,k)+dsny*uy(i,ny-3,k)
   enddo
   enddo
if (iimplicit.eq.1) return  
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
   enddo
   enddo
   enddo
endif

return  
end subroutine deryyt

!********************************************************************
!
subroutine derzzt(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire) 
!
!********************************************************************

USE param 
USE derivZ 

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: sz 
real(mytype), dimension(nz) :: sfz,ssz,swz

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=askzt*(uz(i,j,2)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz  ))&
           +bskzt*(uz(i,j,3)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-1))&
           +cskzt*(uz(i,j,4)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-2))
      rz(i,j,1)=-1.
      tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,1 ))&
           +bskzt*(uz(i,j,4)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz))&
           +cskzt*(uz(i,j,5)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz-1))
      rz(i,j,2)=0.
      tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,2 ))&
           +bskzt*(uz(i,j,5)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,1 ))&
           +cskzt*(uz(i,j,6)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,nz))
      rz(i,j,3)=0.
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-1))&
           +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-2))&
           +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-3))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=askzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-3))&
           +bskzt*(uz(i,j,nz  )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-4))&
           +cskzt*(uz(i,j,1   )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-5))
      rz(i,j,nz-2)=0.
      tz(i,j,nz-1)=askzt*(uz(i,j,nz  )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-2))&
           +bskzt*(uz(i,j,1   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-3))&
           +cskzt*(uz(i,j,2   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-4))
      rz(i,j,nz-1)=0.
      tz(i,j,nz  )=askzt*(uz(i,j,1 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-1))&
           +bskzt*(uz(i,j,2 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-2))&
           +cskzt*(uz(i,j,3 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-3))
      rz(i,j,nz  )=alsakz
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      rz(i,j,nz)=rz(i,j,nz)*swz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if (nclz==1) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=askzt*(uz(i,j,2)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,2))&
              +bskzt*(uz(i,j,3)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,3))&
              +cskzt*(uz(i,j,4)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,4))
         tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,1))&
              +bskzt*(uz(i,j,4)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,2))&
              +cskzt*(uz(i,j,5)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,3))
         tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))&
              +bskzt*(uz(i,j,5)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,1))&
              +cskzt*(uz(i,j,6)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-1))&
              +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-2))&
              +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=askzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-3))&
              +bskzt*(uz(i,j,nz  )-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-4))&
              +cskzt*(uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-5))
         tz(i,j,nz-1)=askzt*(uz(i,j,nz  )-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-2))&
              +bskzt*(uz(i,j,nz-1)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-3))&
              +cskzt*(uz(i,j,nz-2)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-4))
         tz(i,j,nz  )=askzt*(uz(i,j,nz-1)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-1))&
              +bskzt*(uz(i,j,nz-2)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-2))&
              +cskzt*(uz(i,j,nz-3)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-3))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=0.
         tz(i,j,2)=askzt*(uz(i,j,3)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,1))&
              +bskzt*(uz(i,j,4)-uz(i,j,2)&
              -uz(i,j,2)-uz(i,j,2))&
              +cskzt*(uz(i,j,5)-uz(i,j,2)&
              -uz(i,j,2)-uz(i,j,3))
         tz(i,j,3)=askzt*(uz(i,j,4)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))&
              +bskzt*(uz(i,j,5)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,1))&
              +cskzt*(uz(i,j,6)-uz(i,j,3)&
              -uz(i,j,3)-uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-1))&
              +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-2))&
              +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=askzt*( uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-3))&
              +bskzt*( uz(i,j,nz  )-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-4))&
              +cskzt*(-uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-5))
         tz(i,j,nz-1)=askzt*( uz(i,j,nz  )-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-2))&
              +bskzt*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-3))&
              +cskzt*(-uz(i,j,nz-2)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-4))
         tz(i,j,nz  )=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
endif

if (nclz==2) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
           +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
      tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
           -uz(i,j,2)+uz(i,j,1))
      tz(i,j,3)=as3z*(uz(i,j,4)-uz(i,j,3)&
           -uz(i,j,3)+uz(i,j,2))&
           +bs3z*(uz(i,j,5)-uz(i,j,3)&
           -uz(i,j,3)+uz(i,j,1))
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=askzt*(uz(i,j,k+1)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-1))&
           +bskzt*(uz(i,j,k+2)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-2))&
           +cskzt*(uz(i,j,k+3)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-3))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=astz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-3))&
           +bstz*(uz(i,j,nz  )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-4))
      tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-2))
      tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
           +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
   enddo
   enddo
   enddo
endif

return  
end subroutine derzzt

!********************************************************************
!
subroutine scalar_schemes(fpi2t)
!
!********************************************************************
   use decomp_2d, only : mytype
   use param
   use derivX
   use derivY
   use derivZ
   use variables
   use ludecomp

   implicit none

   real(mytype), intent(in) :: fpi2t
   integer :: i,j,k
   real(mytype), dimension(ny) :: aamt10,bbmt10,ccmt10,ddmt10,eemt10,ggmt10,hhmt10,wwmt10,zzmt10
   real(mytype), dimension(ny) :: rrmt10,qqmt10,vvmt10,ssmt10
   real(mytype), dimension(ny) :: aamt11,bbmt11,ccmt11,ddmt11,eemt11,ggmt11,hhmt11,wwmt11,zzmt11
   real(mytype), dimension(ny) :: rrmt11,qqmt11,vvmt11,ssmt11

   xcst_pr=xnu*dt*0.5/sc

   alsaixt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
   asixt  =((6.-9.*alsaixt)/4.)/dx2
   bsixt  =((-3.+24.*alsaixt)/5.)/(4.*dx2)
   csixt  =((2.-11.*alsaixt)/20.)/(9.*dx2)

   alsajyt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
   asjyt  =((6.-9.*alsajyt)/4.)/dy2
   bsjyt  =((-3.+24.*alsajyt)/5.)/(4.*dy2)
   csjyt  =((2.-11.*alsajyt)/20.)/(9.*dy2)

   alsakzt=(45.*fpi2t*pi*pi-272.)/(2.*(45.*fpi2t*pi*pi-208.))
   askzt  =((6.-9.*alsakzt)/4.)/dz2
   bskzt  =((-3.+24.*alsakzt)/5.)/(4.*dz2)
   cskzt  =((2.-11.*alsakzt)/20.)/(9.*dz2)

   if (nclx.eq.0) then
      sfxt(1)   =alsaixt
      sfxt(2)   =alsaixt
      sfxt(nx-2)=alsaixt
      sfxt(nx-1)=alsaixt
      sfxt(nx)  =0.
      scxt(1)   =2.
      scxt(2)   =1.
      scxt(nx-2)=1.
      scxt(nx-1)=1.
      scxt(nx  )=1.+alsaixt*alsaixt
      sbxt(1)   =alsaixt
      sbxt(2)   =alsaixt
      sbxt(nx-2)=alsaixt
      sbxt(nx-1)=alsaixt
      sbxt(nx  )=0.
      do i=3,nx-3
         sfxt(i)=alsaixt
         scxt(i)=1.
         sbxt(i)=alsaixt
      enddo
   endif

   if (nclx.eq.1) then
      sfxt(1)   =alsaixt+alsaixt
      sfxt(2)   =alsaixt
      sfxt(nx-2)=alsaixt
      sfxt(nx-1)=alsaixt
      sfxt(nx)  =0.
      scxt(1)   =1.
      scxt(2)   =1.
      scxt(nx-2)=1.
      scxt(nx-1)=1.
      scxt(nx  )=1.
      sbxt(1)   =alsaixt
      sbxt(2)   =alsaixt
      sbxt(nx-2)=alsaixt
      sbxt(nx-1)=alsaixt+alsaixt
      sbxt(nx  )=0.
      do i=3,nx-3
         sfxt(i)=alsaixt
         scxt(i)=1.
         sbxt(i)=alsaixt
      enddo
   endif

   if (nclx.eq.2) then
      sfxt(1)   =alsa1x
      sfxt(2)   =alsa2x
      sfxt(3)   =alsa3x
      sfxt(nx-3)=alsaixt
      sfxt(nx-2)=alsatx
      sfxt(nx-1)=alsamx
      sfxt(nx)  =0.
      scxt(1)   =1.
      scxt(2)   =1.
      scxt(3)   =1.
      scxt(nx-3)=1.
      scxt(nx-2)=1.
      scxt(nx-1)=1.
      scxt(nx  )=1.
      sbxt(1)   =alsa2x 
      sbxt(2)   =alsa3x
      sbxt(3)   =alsaixt
      sbxt(nx-3)=alsatx
      sbxt(nx-2)=alsamx
      sbxt(nx-1)=alsanx
      sbxt(nx  )=0.
      do i=4,nx-4
         sfxt(i)=alsaixt
         scxt(i)=1.
         sbxt(i)=alsaixt
      enddo
   endif

   if (ncly.eq.0) then
      sfyt(1)   =alsajyt
      sfyt(2)   =alsajyt
      sfyt(ny-2)=alsajyt
      sfyt(ny-1)=alsajyt
      sfyt(ny)  =0.
      scyt(1)   =2.
      scyt(2)   =1.
      scyt(ny-2)=1.
      scyt(ny-1)=1.
      scyt(ny  )=1.+alsajyt*alsajyt
      sbyt(1)   =alsajyt
      sbyt(2)   =alsajyt
      sbyt(ny-2)=alsajyt
      sbyt(ny-1)=alsajyt
      sbyt(ny  )=0.
      do j=3,ny-3
         sfyt(j)=alsajyt
         scyt(j)=1.
         sbyt(j)=alsajyt
      enddo
   endif

   if (ncly.eq.1) then
      sfyt(1)   =alsajyt+alsajyt
      sfyt(2)   =alsajyt
      sfyt(ny-2)=alsajyt
      sfyt(ny-1)=alsajyt
      sfyt(ny)  =0.
      scyt(1)   =1.
      scyt(2)   =1.
      scyt(ny-2)=1.
      scyt(ny-1)=1.
      scyt(ny  )=1.
      sbyt(1)   =alsajyt
      sbyt(2)   =alsajyt
      sbyt(ny-2)=alsajyt
      sbyt(ny-1)=alsajyt+alsajyt
      sbyt(ny  )=0.
      do j=3,ny-3
         sfyt(j)=alsajyt
         scyt(j)=1.
         sbyt(j)=alsajyt
      enddo
   endif

   if (ncly.eq.2) then
      sfyt(1)   =alsa1y
      sfyt(2)   =alsa2y
      sfyt(3)   =alsa3y
      sfyt(ny-3)=alsajyt
      sfyt(ny-2)=alsaty
      sfyt(ny-1)=alsamy
      sfyt(ny)  =0.
      scyt(1)   =1.
      scyt(2)   =1.
      scyt(3)   =1.
      scyt(ny-3)=1.
      scyt(ny-2)=1.
      scyt(ny-1)=1.
      scyt(ny  )=1.
      sbyt(1)   =alsa2y 
      sbyt(2)   =alsa3y
      sbyt(3)   =alsajyt
      sbyt(ny-3)=alsaty
      sbyt(ny-2)=alsamy
      sbyt(ny-1)=alsany
      sbyt(ny  )=0.
      do j=4,ny-4
         sfyt(j)=alsajyt
         scyt(j)=1.
         sbyt(j)=alsajyt
      enddo
   endif

   if (nclz.eq.0) then
      sfzt(1)   =alsakzt
      sfzt(2)   =alsakzt
      sfzt(nz-2)=alsakzt
      sfzt(nz-1)=alsakzt
      sfzt(nz)  =0.
      sczt(1)   =2.
      sczt(2)   =1.
      sczt(nz-2)=1.
      sczt(nz-1)=1.
      sczt(nz  )=1.+alsakzt*alsakzt
      sbzt(1)   =alsakzt
      sbzt(2)   =alsakzt
      sbzt(nz-2)=alsakzt
      sbzt(nz-1)=alsakzt
      sbzt(nz  )=0.
      do k=3,nz-3
         sfzt(k)=alsakzt
         sczt(k)=1.
         sbzt(k)=alsakzt
      enddo
   endif

   if (nclz.eq.1) then
      sfzt(1)   =alsakzt+alsakzt
      sfzt(2)   =alsakzt
      sfzt(nz-2)=alsakzt
      sfzt(nz-1)=alsakzt
      sfzt(nz)  =0.
      sczt(1)   =1.
      sczt(2)   =1.
      sczt(nz-2)=1.
      sczt(nz-1)=1.
      sczt(nz  )=1.
      sbzt(1)   =alsakzt
      sbzt(2)   =alsakzt
      sbzt(nz-2)=alsakzt
      sbzt(nz-1)=alsakzt+alsakzt
      sbzt(nz  )=0.
      do k=3,nz-3
         sfzt(k)=alsakzt
         sczt(k)=1.
         sbzt(k)=alsakzt
      enddo
   endif
   
   if (nclz.eq.2) then
      sfzt(1)   =alsa1z
      sfzt(2)   =alsa2z
      sfzt(3)   =alsa3z
      sfzt(nz-3)=alsakzt
      sfzt(nz-2)=alsatz
      sfzt(nz-1)=alsamz
      sfzt(nz)  =0.
      sczt(1)   =1.
      sczt(2)   =1.
      sczt(3)   =1.
      sczt(nz-3)=1.
      sczt(nz-2)=1.
      sczt(nz-1)=1.
      sczt(nz  )=1.
      sbzt(1)   =alsa2z 
      sbzt(2)   =alsa3z
      sbzt(3)   =alsakzt
      sbzt(nz-3)=alsatz
      sbzt(nz-2)=alsamz
      sbzt(nz-1)=alsanz
      sbzt(nz  )=0.
      do k=4,nz-4
         sfzt(k)=alsakzt
         sczt(k)=1.
         sbzt(k)=alsakzt
      enddo
   endif

   if (nclx.eq.1) then
      sfxt(1)=0.
   endif
   if (ncly.eq.1) then
      sfyt(1)=0.
   endif
   if (nclz.eq.1) then
      sfzt(1)=0.
   endif

   call prepare (sbxt,scxt,sfxt,ssxt,swxt,nx)
   call prepare (sbyt,scyt,sfyt,ssyt,swyt,ny)
   call prepare (sbzt,sczt,sfzt,sszt,swzt,nz)

   call prepare (sbxt,scxt,sfxpt,ssxpt,swxpt,nx)
   call prepare (sbyt,scyt,sfypt,ssypt,swypt,ny)
   call prepare (sbzt,sczt,sfzpt,sszpt,swzpt,nz)

   if (nclx.eq.1) then
      sbxt(nx-1)=0.
      call prepare (sbxt,scxt,sfxt,ssxt,swxt,nx)
   endif
   if (ncly.eq.1) then
      sbyt(ny-1)=0.
      call prepare (sbyt,scyt,sfyt,ssyt,swyt,ny)
   endif
   if (nclz.eq.1) then
      sbzt(nz-1)=0.
      call prepare (sbzt,sczt,sfzt,sszt,swzt,nz)
   endif

!!!!!!!!!!!!!!!!!!!!!!
!CAS SEPTADIAGONALE !!
!!!!!!!!!!!!!!!!!!!!!!

!!! NCL = 2, dirichlet imposé, fonction nulle à la paroi
!
!DIAG
aamt(1     )=as1y
aamt(ny    )=asny
aamt(2     )=-2.*as2y
aamt(ny-1  )=-2.*asmy
aamt(3     )=-2.*(as3y+bs3y)
aamt(ny-2  )=-2.*(asty+bsty)
aamt(4:ny-3)=-2.*(asjyt+bsjyt+csjyt)
if (istret==0) then
   aamt = 1.-xcst_pr*aamt
else
   aamt = 1./pp2y-xcst_pr*aamt
endif
!CL sur aamt
if (istret==0) then
   aamt(1 )=alpha_0+beta_0*(11./6./dy)
   aamt(ny)=alpha_n+beta_n*(11./6./dy)
else
   aamt(1 )=alpha_0+beta_0*ppy(1 )*(11./6./dy)
   aamt(ny)=alpha_n+beta_n*ppy(ny)*(11./6./dy)
endif
!
!DIAG SUP 1
bbmt(1     )=bs1y
bbmt(ny    )=bsny
bbmt(2     )=as2y
bbmt(ny-1  )=asmy
bbmt(3     )=as3y
bbmt(ny-2  )=asty
bbmt(4:ny-3)=asjyt
bbmt = -xcst_pr*bbmt
if (istret==0) then
  bbmt(2     )=bbmt(2     )+alsa2y
  bbmt(ny-1  )=bbmt(ny-1  )+alsamy
  bbmt(3     )=bbmt(3     )+alsa3y
  bbmt(ny-2  )=bbmt(ny-2  )+alsaty
  bbmt(4:ny-3)=bbmt(4:ny-3)+alsajyt
else
  bbmt(2     )=bbmt(2     )+alsa2y/pp2y(3)
  bbmt(ny-1  )=bbmt(ny-1  )+alsamy/pp2y(ny)
  bbmt(3     )=bbmt(3     )+alsa3y/pp2y(4)
  bbmt(ny-2  )=bbmt(ny-2  )+alsaty/pp2y(ny-1)
  bbmt(4:ny-3)=bbmt(4:ny-3)+alsajyt/pp2y(5:ny-2)
endif
!CL sur bbmt
if (istret==0) then
   bbmt(1 )=beta_0*(-18./6./dy)
   bbmt(ny)=0.
else
   bbmt(1 )=beta_0*ppy(1)*(-18./6./dy)
   bbmt(ny)=0.
endif
!
!DIAG SUP 2
ccmt(1     )=cs1y
ccmt(ny    )=csny
ccmt(2     )=0.
ccmt(ny-1  )=0.
ccmt(3     )=bs3y
ccmt(ny-2  )=bsty
ccmt(4:ny-3)=bsjyt
ccmt = -xcst_pr*ccmt
!CL sur ccmt
if (istret==0) then
   ccmt(1 )=beta_0*(9./6./dy)
   ccmt(ny)=0.
else
   ccmt(1 )=beta_0*ppy(1)*(9./6./dy)
   ccmt(ny)=0.
endif
!
!DIAG SUP 3
rrmt(1     )=ds1y
rrmt(ny    )=dsny
rrmt(2     )=0. 
rrmt(ny-1  )=0.
rrmt(3     )=0.
rrmt(ny-2  )=0.
rrmt(4:ny-3)=csjyt
rrmt = -xcst_pr*rrmt
!CL sur rrmt
if (istret==0) then
   rrmt(1 )=beta_0*(-2./6./dy)
   rrmt(ny)=0.
else
   rrmt(1 )=beta_0*ppy(1)*(-2./6./dy)
   rrmt(ny)=0.
endif
!
!DIAG INF 1
if (istret==0) then
  ddmt=bbmt
else
  ddmt(1     )=bs1y
  ddmt(ny    )=bsny
  ddmt(2     )=as2y
  ddmt(ny-1  )=asmy
  ddmt(3     )=as3y
  ddmt(ny-2  )=asty
  ddmt(4:ny-3)=asjyt
  ddmt = -xcst_pr*ddmt
  ddmt(2     )=ddmt(2     )+alsa2y/pp2y(1)
  ddmt(ny-1  )=ddmt(ny-1  )+alsamy/pp2y(ny-2)
  ddmt(3     )=ddmt(3     )+alsa3y/pp2y(2)
  ddmt(ny-2  )=ddmt(ny-2  )+alsaty/pp2y(ny-3)
  ddmt(4:ny-3)=ddmt(4:ny-3)+alsajyt/pp2y(3:ny-4)
endif
!CL sur ddmt
if (istret==0) then
   ddmt(1 )=0.
   ddmt(ny)=beta_n*(-18./6./dy)
else
   ddmt(1 )=0.
   ddmt(ny)=beta_n*ppy(ny)*(-18./6./dy)
endif
!
!DIAG INF 2
eemt=ccmt
!CL sur eemt
if (istret==0) then
   eemt(1)=0.
   eemt(ny)=beta_n*(9./6./dy)
else
   eemt(1)=0.
   eemt(ny)=beta_n*ppy(ny)*(9./6./dy)
endif
!
!DIAG INF 3
qqmt=rrmt
!CL sur qqmt
if (istret==0) then
   qqmt(1)=0.
   qqmt(ny)=beta_n*(-2./6./dy)
else
   qqmt(1)=0.
   qqmt(ny)=beta_n*ppy(ny)*(-2./6./dy)
endif

!!! NCL = 1, npaire=0, dirichlet imposé, fonction impaire
!
! DIAG
aamt10(1     )=0.
aamt10(ny    )=0.
aamt10(2     )=-2.*asjyt-3.*bsjyt-2.*csjyt
aamt10(ny-1  )=aamt10(2)
aamt10(3:ny-2)=-2.*(asjyt+bsjyt+csjyt)
if (istret==0) then
   aamt10 = 1. - xcst_pr*aamt10
else
   aamt10 = 1./pp2y - xcst_pr*aamt10
endif
!
!DIAG SUP 1
bbmt10(1     )=0.
bbmt10(ny    )=0.
bbmt10(2     )=asjyt-csjyt
bbmt10(ny-1  )=asjyt
bbmt10(3     )=asjyt
bbmt10(ny-2  )=asjyt-csjyt
bbmt10(4:ny-3)=asjyt
if (istret==0) then
   bbmt10 = alsajyt - xcst_pr*bbmt10
else
   bbmt10(2:ny-1) = alsajyt/pp2y(3:ny) - xcst_pr*bbmt10(2:ny-1)
endif
!CL sur bbm10t
bbmt10(1 )=0.
bbmt10(ny)=0.
!
!DIAG SUP 2
ccmt10(1     )=0.
ccmt10(ny    )=0.
ccmt10(2     )=bsjyt
ccmt10(ny-1  )=0.
ccmt10(3     )=bsjyt
ccmt10(ny-2  )=bsjyt
ccmt10(4:ny-3)=bsjyt
ccmt10 = -xcst_pr*ccmt10
!
!DIAG SUP 3
rrmt10(1     )=0.
rrmt10(ny    )=0.
rrmt10(2     )=csjyt
rrmt10(ny-1  )=0.
rrmt10(3     )=csjyt
rrmt10(ny-2  )=0.
rrmt10(4:ny-3)=csjyt
rrmt10 = -xcst_pr*rrmt10
!
!DIAG INF 1
ddmt10(1     )=0.
ddmt10(ny    )=0.
ddmt10(2     )=asjyt
ddmt10(ny-1  )=asjyt-csjyt
ddmt10(3     )=asjyt-csjyt
ddmt10(ny-2  )=asjyt
ddmt10(4:ny-3)=asjyt
if (istret==0) then
   ddmt10 = alsajyt - xcst_pr*ddmt10
else
   ddmt10(2:ny-1) = alsajyt/pp2y(1:ny-2) - xcst_pr*ddmt10(2:ny-1)
endif
!CL sur ddmt10
ddmt10(1 )=0.
ddmt10(ny)=0.
!
!DIAG INF 2
eemt10(1     )=0.
eemt10(ny    )=0.
eemt10(2     )=0.
eemt10(ny-1  )=bsjyt
eemt10(3     )=bsjyt
eemt10(ny-2  )=bsjyt
eemt10(4:ny-3)=bsjyt
eemt10 = -xcst_pr*eemt10
!
!DIAG INF 3
qqmt10(1     )=0.
qqmt10(ny    )=0.
qqmt10(2     )=0.
qqmt10(ny-1  )=csjyt
qqmt10(3     )=0.
qqmt10(ny-2  )=csjyt
qqmt10(4:ny-3)=csjyt
qqmt10 = -xcst_pr*qqmt10

!!! NCL = 1, npaire=1, neumann imposé, fonction paire
!
! DIAG
aamt11(1     )=-2.*(asjyt+bsjyt+csjyt)
aamt11(ny    )=aamt11(1)
aamt11(2     )=-2.*asjyt-bsjyt-2.*csjyt
aamt11(ny-1  )=aamt11(2)
aamt11(3:ny-2)=-2.*(asjyt+bsjyt+csjyt)
if (istret==0) then
   aamt11 = 1. - xcst_pr*aamt11
else
   aamt11 = 1./pp2y - xcst_pr*aamt11
endif
!
!DIAG SUP 1
bbmt11(1     )=2.*asjyt
bbmt11(ny    )=0.
bbmt11(2     )=asjyt+csjyt
bbmt11(ny-1  )=asjyt
bbmt11(3     )=asjyt
bbmt11(ny-2  )=asjyt+csjyt
bbmt11(4:ny-3)=asjyt
if (istret==0) then
   bbmt11 = alsajyt - xcst_pr*bbmt11
else
   bbmt11(1:ny-1) = alsajyt/pp2y(2:ny) - xcst_pr*bbmt11(1:ny-1)
endif
!CL sur bbm11t
if (istret==0) then
   bbmt11(1 )=bbmt11(1)+alsajyt
else
   bbmt11(1 )=bbmt11(1)+alsajyt/pp2y(2)
endif
bbmt11(ny)=0.
!
!DIAG SUP 2
ccmt11(1     )=2.*bsjyt
ccmt11(ny    )=0.
ccmt11(2     )=bsjyt
ccmt11(ny-1  )=0.
ccmt11(3     )=bsjyt
ccmt11(ny-2  )=bsjyt
ccmt11(4:ny-3)=bsjyt
ccmt11 = -xcst_pr*ccmt11
!
!DIAG SUP 3
rrmt11(1     )=2.*csjyt
rrmt11(ny    )=0.
rrmt11(2     )=csjyt
rrmt11(ny-1  )=0.
rrmt11(3     )=csjyt
rrmt11(ny-2  )=0.
rrmt11(4:ny-3)=csjyt
rrmt11 = -xcst_pr*rrmt11
!
!DIAG INF 1
ddmt11(1     )=0.
ddmt11(ny    )=2.*asjyt
ddmt11(2     )=asjyt
ddmt11(ny-1  )=asjyt+csjyt
ddmt11(3     )=asjyt+csjyt
ddmt11(ny-2  )=asjyt
ddmt11(4:ny-3)=asjyt
if (istret==0) then
   ddmt11 = alsajyt - xcst_pr*ddmt11
else
   ddmt11(2:ny) = alsajyt/pp2y(1:ny-1) - xcst_pr*ddmt11(2:ny)
endif
!CL sur ddmt11
ddmt11(1 )=0.
if (istret==0) then
   ddmt11(ny)=ddmt11(ny)+alsajyt!a1
else
   ddmt11(ny)=ddmt11(ny)+alsajyt/pp2y(ny-1)!a1
endif
!
!DIAG INF 2
eemt11(1     )=0.
eemt11(ny    )=2.*bsjyt
eemt11(2     )=0.
eemt11(ny-1  )=bsjyt
eemt11(3     )=bsjyt
eemt11(ny-2  )=bsjyt
eemt11(4:ny-3)=bsjyt
eemt11 = -xcst_pr*eemt11
!
!DIAG INF 3
qqmt11(1     )=0.
qqmt11(ny    )=2.*csjyt
qqmt11(2     )=0.
qqmt11(ny-1  )=csjyt
qqmt11(3     )=0.
qqmt11(ny-2  )=csjyt
qqmt11(4:ny-3)=csjyt
qqmt11 = -xcst_pr*qqmt11

!!! NXL = 0
!DIAG
if (istret==0) then
   aamt0 = 1.-xcst_pr*(-2.*(asjyt+bsjyt+csjyt))
else
   aamt0 = 1./pp2y-xcst_pr*(-2.*(asjyt+bsjyt+csjyt))
endif
!
!DIAG SUP 1
if (istret==0) then
   bbmt0 = alsajyt-xcst_pr*asjyt
else
   bbmt0(1:ny-1) = alsajyt/pp2y(2:ny)-xcst_pr*asjyt
   bbmt0(ny) = alsajyt/pp2y(1)-xcst_pr*asjyt
endif
!
!DIAG SUP 2
ccmt0 = -xcst_pr*bsjyt
!
!DIAG SUP 3
rrmt0 = -xcst_pr*csjyt
!
!DIAG INF 1
if (istret==0) then
   ddmt0=bbmt0
else
   ddmt0(1)=alsajyt/pp2y(ny) -xcst_pr*asjyt
   ddmt0(2:ny)=alsajyt/pp2y(1:ny-1) -xcst_pr*asjyt
endif
!
!DIAG INF 2
eemt0=ccmt0
!
!DIAG INF 3
qqmt0=rrmt0

call ludecomp7(aamt,bbmt,ccmt,ddmt,eemt,qqmt,ggmt,hhmt,ssmt,rrmt,&
     vvmt,wwmt,zzmt,ny)
call ludecomp7(aamt10,bbmt10,ccmt10,ddmt10,eemt10,qqmt10,ggmt10,hhmt10,ssmt10,rrmt10,&
     vvmt10,wwmt10,zzmt10,ny)
call ludecomp7(aamt11,bbmt11,ccmt11,ddmt11,eemt11,qqmt11,ggmt11,hhmt11,ssmt11,rrmt11,&
     vvmt11,wwmt11,zzmt11,ny)
call ludecomp7(aamt0,bbmt0,ccmt0,ddmt0,eemt0,qqmt0,ggmt0,hhmt0,ssmt0,rrmt0,&
     vvmt0,wwmt0,zzmt0,l1mt,l2mt,l3mt,u1mt,u2mt,u3mt,ny)

! npaire=0
   aamt1=aamt10; bbmt1=bbmt10; ccmt1=ccmt10; ddmt1=ddmt10
   eemt1=eemt10; ggmt1=ggmt10; hhmt1=hhmt10; wwmt1=wwmt10; zzmt1=zzmt10
   rrmt1=rrmt10; qqmt1=qqmt10; vvmt1=vvmt10; ssmt1=ssmt10
!! npaire=1
!   aamt1=aamt11; bbmt1=bbmt11; ccmt1=ccmt11; ddmt1=ddmt11
!   eemt1=eemt11; ggmt1=ggmt11; hhmt1=hhmt11; wwmt1=wwmt11; zzmt1=zzmt11
!   rrmt1=rrmt11; qqmt1=qqmt11; vvmt1=vvmt11; ssmt1=ssmt11

end subroutine scalar_schemes
