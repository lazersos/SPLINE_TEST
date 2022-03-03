!gfortran -o xspline_test -I/Users/slazerso/bin/libstell_dir spline_test.f90 ~/bin/libstell.a
PROGRAM SPLINE_TEST
   USE EZspline
   USE EZspline_obj
   IMPLICIT NONE

   INTEGER :: i, j, ier, nhigh
   DOUBLE PRECISION :: u,v, ftemp, fact, factu, factv, factuv, F1d_grad
   DOUBLE PRECISION, DIMENSION(2) :: F2d_grad
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x1, x2, f1d, fu1d
   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f2d, fu2d, fv2d, fuv2d

   TYPE(EZspline1_r8) :: F1d_spl
   TYPE(EZspline2_r8) :: F2d_spl, F2d_spl_mod

   INTEGER, PARAMETER :: nx1 = 32
   INTEGER, PARAMETER :: nx2 = 32
   INTEGER, PARAMETER :: bcs1(2) = (/-1,-1/)
   DOUBLE PRECISION, PARAMETER :: pi2 = 6.283185482025146D+00


   ALLOCATE(x1(nx1),x2(nx2),f1d(nx1),fu1d(nx1))
   ALLOCATE(f2d(nx1,nx2),fu2d(nx1,nx2),fv2d(nx1,nx2),fuv2d(nx1,nx2))

   do i = 1, nx1
      do j = 1, nx2
         u = 2.0*(DBLE(i-1)/DBLE(nx1-1)-0.5)
         v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
         CALL F_VAL_STEP(u,v,fact,factu,factv,factuv)
         f2d(i,j) = fact
         fu2d(i,j) = factu
         fv2d(i,j) = factv
         fuv2d(i,j) = factuv
      enddo
   end do
   j = nx2/2
   v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
   do i = 1, nx1
      u = 2.0*(DBLE(i-1)/DBLE(nx1-1)-0.5)
      CALL F_VAL_STEP(u,v,fact,factu,factv,factuv)
      f1d(i) = fact
      fu1d(i) = factu
   enddo

   CALL EZspline_init(F1d_spl,nx1,bcs1,ier)
   CALL EZspline_init(F2d_spl,nx1,nx2,bcs1,bcs1,ier)
   CALL EZspline_init(F2d_spl_mod,nx1,nx2,bcs1,bcs1,ier)

   do i = 1, nx1
      u = 2.0*(DBLE(i-1)/DBLE(nx1-1)-0.5)
      F2d_spl%x1(i)=u
      F2d_spl_mod%x1(i)=u
      F1d_spl%x1(i)=u
   enddo

   do j = 1, nx2
      v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
      F2d_spl%x2(j)=v
      F2d_spl_mod%x2(j)=v
   enddo

   F1d_spl%isHermite = 1
   F2d_spl%isHermite = 1
   F2d_spl_mod%isHermite = 1

   CALL EZspline_setup(F1d_spl,f1d,ier,EXACT_DIM=.TRUE.)
   CALL EZspline_setup(F2d_spl,f2d,ier,EXACT_DIM=.TRUE.)
   CALL EZspline_setup(F2d_spl_mod,f2d,ier,EXACT_DIM=.TRUE.)

   !F_spl_mod%fspl(2,:,:)=fu
   !F_spl_mod%fspl(3,:,:)=fv
   !F_spl_mod%fspl(4,:,:)=fuv

   nhigh  = nx1*16
   do i = 1, nhigh
      do j = 1, nhigh
         u = 2.0*(DBLE(i-1)/DBLE(nhigh-1)-0.5)
         v = 2.0*(DBLE(j-1)/DBLE(nhigh-1)-0.5)
         CALL F_VAL_STEP(u,v,fact,factu,factv,factuv)
         CALL EZspline_interp(F2d_spl,u,v,ftemp,ier)
         CALL EZspline_gradient(F2d_spl,u,v,F2d_grad,ier)
         WRITE(326,'(2(1X,I6.6),8(1X,E22.12))') i,j,u,v,fact,factu,factv,ftemp,F2d_grad(1),F2d_grad(2)
      enddo
   enddo
   CALL FLUSH(326)

   j=nx2/2
   v = 2.0*(DBLE(j-1)/DBLE(nx2-1)-0.5)
   do i = 1, nhigh
      u = 2.0*(DBLE(i-1)/DBLE(nhigh-1)-0.5)
      CALL F_VAL_STEP(u,v,fact,factu,factv,factuv)
      CALL EZspline_interp(F1d_spl,u,ftemp,ier)
      CALL EZspline_gradient(F1d_spl,u,F1d_grad,ier)
      WRITE(327,'(2(1X,I6.6),8(1X,E22.12))') i,j,u,v,fact,factu,factv,ftemp,F1d_grad
   end do
   CALL FLUSH(327)




   DEALLOCATE(x1,x2,f1d,fu1d)
   DEALLOCATE(f2d,fu2d,fv2d,fuv2d)


CONTAINS 

   SUBROUTINE F2_VAL(u,v,f,fu,fv,fuv)
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
            f   = f   + A + B*cos(pi2*m*u+pi2*n*v)
            fu  = fu  - B*m*pi2*sin(pi2*m*u+pi2*n*v) 
            fv  = fv  - B*n*pi2*sin(pi2*m*u+pi2*n*v) 
            fuv = fuv + B*m*pi2*n*pi2*cos(pi2*m*u+pi2*n*v)
         enddo
      enddo
   RETURN
   END SUBROUTINE F2_VAL

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

   SUBROUTINE F_VAL_STEP(u,v,f,fu,fv,fuv)
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
      IF (u==0) fu = fu + 100
      IF (u==0) fuv = fuv + 100
      IF (u>=0) f = f+100
   RETURN
   END SUBROUTINE F_VAL_STEP

END PROGRAM SPLINE_TEST
