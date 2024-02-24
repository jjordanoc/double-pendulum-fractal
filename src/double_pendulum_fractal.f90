module double_pendulum_fractal
   use iso_fortran_env
   use fplot_core
   use forcolormap, only : colormaps_list
   use omp_lib
   implicit none
   private

   public :: main_loop
contains
   subroutine main_loop
      logical, parameter :: read_from_file = .false.
      real(real128), parameter :: dt = 0.1
      integer iter; real(real128) y(5), E0
      real(real64), parameter :: PI=4.Q0*DATAN(1.D0)
      integer, parameter :: iterations = 100000
      integer(int32), parameter :: n = 500

      ! Plot Variables
      type(surface_plot) :: plt
      type(surface_plot_data) :: d1
      class(plot_axis), pointer :: xAxis, yAxis, zAxis
      type(custom_colormap) :: map
      type(cmap) :: colors

      ! Problem variables
      real(real64), dimension(n, n, 2), target :: space
      real(real64), pointer, dimension(:,:) :: th1_0, th2_0
      real(real64), dimension(n, n) :: tflip
      integer(int32) :: i, j

      ! Set up the colormap
      call colors%set("glasgow", 0.0d0, 10.0d0)
      call map%set_colormap(colors)

      ! Initialize the plot
      call plt%initialize(term=GNUPLOT_TERMINAL_QT)
      call plt%set_colormap(map)

      ! Set the orientation of the plot
      call plt%set_elevation(0.0d0)
      call plt%set_azimuth(0.0d0)

      ! Define titles
      call plt%set_title("Double pendulum fractal")

      xAxis => plt%get_x_axis()
      call xAxis%set_title("Theta 1")

      yAxis => plt%get_y_axis()
      call yAxis%set_title("Theta 2")

      zAxis => plt%get_z_axis()
      call zAxis%set_title("Time to flip")

      ! Define the data
      space = meshgrid(linspace(-PI, PI, n), linspace(-PI, PI, n))
      th1_0 => space(:,:,1)
      th2_0 => space(:,:,2)
      
      !$OMP PARALLEL PRIVATE(y, iter) SHARED(th1_0, th2_0, tflip)
      !$omp do
      do i = 1, n
         do j = 1, n
            y(1) = 0.0 ! t = 0
            y(2) = th1_0(i,j) ! th1 = th1_0
            y(3) = th2_0(i,j) ! th2 = th2_0
            y(4) = 0.0 ! pth1 = 0
            y(5) = 0.0 ! pht2 = 0
            !main evolution loop
            do iter = 0,iterations
               ! pendulum flips when angle is greater than pi in absolute value
               if (abs(y(2)) > PI .or. abs(y(3)) > PI) then
                  exit
               end if
               call gl8(y, dt)
            end do
            tflip(i,j) = y(1)
            write (*,*) 'Final time: ', i, j, y(1)
         end do
      end do
      !$OMP END DO
      !$OMP END PARALLEL
      call d1%define_data(th1_0, th2_0, log(tflip + 1))
      call plt%set_use_map_view(.true.)
      call plt%set_show_gridlines(.false.)
      call plt%push(d1)
      call plt%draw()
   end subroutine main_loop

   ! function energy(vec)
   !    real energy, vec(5)
   !    ! H = p^2/2 + x^4/4
   !    real, parameter :: m = 1.0, g = 9.81, l = 1.0
   !    real :: th1, th2
   !    th1 = vec(2)
   !    th2 = vec(3)
   !    energy = (-1.0/2.0)*m*g*l*(3*cos(th1)+cos(th2))
   ! end function energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine evalf(y, dydx) ! Here you put the odes to solve
      real(real128) y(5), dydx(5)
      ! y = [t, th1, th2, pth1, pth2]
      real(real128), parameter :: m = 1.0, g = 9.81, l = 1.0
      real(real128) :: th1, th2, pth1, pth2
      real(real128) :: ml_2
      ml_2 = m * l**2
      th1 = y(2)
      th2 = y(3)
      pth1 = y(4)
      pth2 = y(5)
      ! dt/dt = 1
      dydx(1) = 1.0
      ! dth1/dt = 6/ml^2 * (2pth1 - 3cos(th1-th2)pth2)/(16 - 9cos^2(th1-th2))
      dydx(2) = (6.0/ml_2) * (2.0*pth1-3.0*cos(th1-th2)*pth2)/(16.0-9.0*(cos(th1-th2)**2))
      ! dth2/dt = 6/ml^2 * (2pth2 - 3cos(th1-th2)pth1)/(16 - 9cos^2(th1-th2))
      dydx(3) = (6.0/ml_2) * (8.0*pth2-3.0*cos(th1-th2)*pth1)/(16.0-9.0*(cos(th1-th2)**2))
      ! dpth1/dt = -1/2ml_2 * (3*g/l*sin(th1)+dth1/dt*dth2/dt*sin(th1-th2))
      dydx(4) = (-1.0/2.0)*ml_2*(3.0*(g/l)*sin(th1)+dydx(2)*dydx(3)*sin(th1-th2))
      ! dpth2/dt = -1/2ml_2 * (3*g/l*sin(th2)+dth1/dt*dth2/dt*sin(th1-th2))
      dydx(5) = (-1.0/2.0)*ml_2*((g/l)*sin(th2)-dydx(2)*dydx(3)*sin(th1-th2))
   end subroutine evalf




! 8th order implicit Gauss-Legendre integrator
! only change n -> the number of equations
   subroutine gl8(y, dt)
      integer, parameter :: s = 4, n = 5
      real(real128) y(n), g(n,s), dt; integer i, k, j

      ! Butcher tableau for 8th order Gauss-Legendre method
      real(real128), parameter :: a(s,s) = reshape((/ &
         0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
         0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
         0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
         -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
         0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
         0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
         0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
         0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
      real(real128), parameter ::   b(s) = (/ &
         0.173927422568726928686531974610999704Q0,   0.326072577431273071313468025389000296Q0,  &
         0.326072577431273071313468025389000296Q0,   0.173927422568726928686531974610999704Q0  /)

      ! iterate trial steps
      g = 0.0; do k = 1,16
         g = matmul(g,a)
         do i = 1,s
            call evalf(y + g(:,i)*dt, g(:,i))
         end do
      end do

      ! update the solution
      y = y + matmul(g,b)*dt
   end subroutine gl8
end module double_pendulum_fractal
