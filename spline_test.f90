!gfortran -o xspline_test -I/Users/slazerso/bin/libstell_dir spline_test.f90 ~/bin/libstell.a
PROGRAM SPLINE_TEST
   USE EZspline
   USE EZspline_obj
   IMPLICIT NONE

   INTEGER :: i, j, ier, nhigh
   DOUBLE PRECISION :: u,v, ftemp, fact, factu, factv, factuv
   DOUBLE PRECISION, DIMENSION(2) :: F_grad
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x1, x2
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f, fu, fv, fuv

   TYPE(EZspline2_r8) :: F_spl, F_spl_mod

   INTEGER, PARAMETER :: nx1 = 8
   INTEGER, PARAMETER :: nx2 = 8
   INTEGER, PARAMETER :: bcs1(2) = (/-1,-1/)
   DOUBLE PRECISION, PARAMETER :: pi2 = 6.283185482025146D+00


   ALLOCATE(x1(nx1),x2(nx2))
   ALLOCATE(f(nx1,nx2),fu(nx1,nx2),fv(nx1,nx2),fuv(nx1,nx2))

   do i = 1, nx1
      do j = 1, nx2
         u = 2.0*(DBLE(i-1)/DBLE(nx1-1)-0.5)
         v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
         CALL F_VAL(u,v,fact,factu,factv,factuv)
         f(i,j) = fact
         fu(i,j) = factu
         fv(i,j) = factv
         fuv(i,j) = factuv
      enddo
   enddo

   CALL EZspline_init(F_spl,nx1,nx2,bcs1,bcs1,ier)
   CALL EZspline_init(F_spl_mod,nx1,nx2,bcs1,bcs1,ier)

   do i = 1, nx1
      u = 2.0*(DBLE(i-1)/DBLE(nx1-1)-0.5)
      F_spl%x1(i)=u
      F_spl_mod%x1(i)=u
   enddo

   do j = 1, nx2
      v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
      F_spl%x2(j)=v
      F_spl_mod%x2(j)=v
   enddo

   F_spl%isHermite = 0
   F_spl_mod%isHermite = 0

   CALL EZspline_setup(F_spl,f,ier)
   CALL EZspline_setup(F_spl_mod,f,ier)

   !F_spl_mod%fspl(2,:,:)=fu
   !F_spl_mod%fspl(3,:,:)=fv
   !F_spl_mod%fspl(4,:,:)=fuv

   nhigh  = nx1*16
   do i = 1, nhigh
      do j = 1, nhigh
         u = 2.0*(DBLE(i-1)/DBLE(nhigh-1)-0.5)
         v = 2.0*(DBLE(j-1)/DBLE(nhigh-1)-0.5)
         CALL F_VAL(u,v,fact,factu,factv,factuv)
         CALL EZspline_interp(F_spl,u,v,ftemp,ier)
         CALL EZspline_gradient(F_spl,u,v,F_grad,ier)
         WRITE(6,'(2(1X,I6.6),8(1X,E22.12))') i,j,u,v,fact,factu,factv,ftemp,F_grad(1),F_grad(2)
      enddo
   enddo
   CALL FLUSH(6)




   DEALLOCATE(x1,x2)
   DEALLOCATE(f,fu,fv,fuv)


CONTAINS 

   SUBROUTINE F_VAL(u,v,f,fu,fv,fuv)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: u,v
      DOUBLE PRECISION, INTENT(out) :: f,fu,fv,fuv
      INTEGER :: m, n
      DOUBLE PRECISION, PARAMETER :: A = 1.0
      DOUBLE PRECISION, PARAMETER :: B = 0.25
      DOUBLE PRECISION, PARAMETER :: pi2 = 6.283185482025146D+00
      INTEGER, PARAMETER :: MPOL = 1
      INTEGER, PARAMETER :: NTOR = 1

      f = 0; fu = 0; fv = 0; fuv=0;
      do m = 1, MPOL
         do n = 1, NTOR
            f   = f   + A + B*cos(pi2*m*u+pi2*n*v) + 8.0*(v-2.0)*(u-2.0)
            fu  = fu  - B*m*pi2*sin(pi2*m*u+pi2*n*v) + 8.0*(v-2.0)
            fv  = fv  - B*n*pi2*sin(pi2*m*u+pi2*n*v) + 8.0*(u-2.0)
            fuv = fuv + B*m*pi2*n*pi2*cos(pi2*m*u+pi2*n*v) + 8.0
         enddo
      enddo
   RETURN
   END SUBROUTINE F_VAL

END PROGRAM SPLINE_TEST
