module double_pendulum_fractal
   use iso_fortran_env
   use fplot_core
   use forcolormap, only : colormaps_list
   implicit none
   private

   public :: main_loop
contains
   subroutine main_loop
      real, parameter :: dt = 0.1
      integer iter; real y(5), E0
      character(len=*), parameter :: OUT_FILE = 'data.txt' ! Output file.
      character (len=*), parameter :: format = '(6g24.16)'
      real(real64), parameter :: PI=4.D0*DATAN(1.D0)
      integer, parameter :: iterations = 100000

      ! Parameters
      integer(int32), parameter :: n = 25

      ! Local Variables
      ! real(real64), dimension(n, n, 2), target :: space
      real(real64), dimension(n) :: th1_0_space, th2_0_space
      real(real64), dimension(n**2) :: th1_0, th2_0, tflip
      type(plot_3d) :: plt
      type(plot_data_3d) :: d1
      class(plot_axis), pointer :: xAxis, yAxis, zAxis
      type(custom_colormap) :: map
      type(cmap) :: colors
      integer(int32) :: i, j, k

      open (1, file=OUT_FILE, status='replace')

      ! Set up the colormap
      call colors%set("glasgow", -8.0d0, 8.0d0)
      call map%set_colormap(colors)

      ! Define the data
      ! space = meshgrid(linspace(-3.0d0, 3.0d0, n), linspace(-3.0d0, 3.0d0, n))
      th1_0_space = linspace(-3.0d0, 3.0d0, n)
      th2_0_space = linspace(-3.0d0, 3.0d0, n)

      ! Initialize the plot
      call plt%initialize(term=GNUPLOT_TERMINAL_QT)
      call plt%set_colormap(map)

      ! Set the orientation of the plot
      call plt%set_elevation(20.0d0)
      call plt%set_azimuth(30.0d0)

      ! Establish lighting
      ! call plt%set_use_lighting(.true.)

      ! Define titles
      call plt%set_title("Double pendulum fractal")

      xAxis => plt%get_x_axis()
      call xAxis%set_title("Theta 1")

      yAxis => plt%get_y_axis()
      call yAxis%set_title("Theta 2")

      zAxis => plt%get_z_axis()
      call zAxis%set_title("Time to flip")

      k = 1
      do i = 1, n
         do j = 1, n
            ! condiciones iniciales
            y(1) = 0.0 ! t = 0
            y(2) = th1_0_space(i) ! th1 = th1_0
            y(3) = th2_0_space(j) ! th2 = th2_0
            y(4) = 0.0 ! pth1 = 0
            y(5) = 0.0 ! pht2 = 0
            !main evolution loop
            do iter = 0,iterations
               ! output positions, errors in energy

               ! write (1,format)  y(1), y(2), E0-energy(y), y(3)
               ! pendulum flips when angle is greater than pi in absolute value
               if (abs(y(2)) > PI .or. abs(y(3)) > PI) then
                  exit
               end if
               call gl8(y, dt)
            end do
            th1_0(k) = th1_0_space(i)
            th2_0(k) = th2_0_space(j)
            tflip(k) = y(1)
            write (1,*)  th1_0(k), th2_0(k), tflip(k)
            k = k + 1
            write (*,*) 'Final time: ', i, j, y(1)
         end do
      end do
      call d1%define_data(th1_0, th2_0, log(tflip + 1))
      ! call plt%set_use_map_view(.true.)
      call plt%push(d1)
      call plt%draw()
      ! E0 = energy(y)

      close(1)
   end subroutine main_loop

   function energy(vec)
      real energy, vec(5)
      ! H = p^2/2 + x^4/4
      real, parameter :: m = 1.0, g = 9.81, l = 1.0
      real :: th1, th2
      th1 = vec(2)
      th2 = vec(3)
      energy = (-1.0/2.0)*m*g*l*(3*cos(th1)+cos(th2))
   end function energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! implicit Gauss-Legendre methods; symplectic with arbitrary Hamiltonian, A-stable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! equations of motion simple anharmonic oscillator

   subroutine evalf(y, dydx)
      real y(5), dydx(5)                       ! Do I change it?
      !YES, here you put the odes to solve

      ! y = [t, th1, th2, pth1, pth2]
      ! real, parameter :: m = 1.0, l = 1.0
      real, parameter :: m = 1.0, g = 9.81, l = 1.0
      real :: th1, th2, pth1, pth2
      real :: ml_2
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
      real y(n), g(n,s), dt; integer i, k, j

      ! Butcher tableau for 8th order Gauss-Legendre method

      ! esta creando una matriz de sxs creando un array y luego reshapeandolo
      real, parameter :: a(s,s) = reshape((/ &
         0.869637112843634643432659873054998518Q-1, -0.266041800849987933133851304769531093Q-1, &
         0.126274626894047245150568805746180936Q-1, -0.355514968579568315691098184956958860Q-2, &
         0.188118117499868071650685545087171160Q0,   0.163036288715636535656734012694500148Q0,  &
         -0.278804286024708952241511064189974107Q-1,  0.673550059453815551539866908570375889Q-2, &
         0.167191921974188773171133305525295945Q0,   0.353953006033743966537619131807997707Q0,  &
         0.163036288715636535656734012694500148Q0,  -0.141906949311411429641535704761714564Q-1, &
         0.177482572254522611843442956460569292Q0,   0.313445114741868346798411144814382203Q0,  &
         0.352676757516271864626853155865953406Q0,   0.869637112843634643432659873054998518Q-1 /), (/s,s/))
      real, parameter ::   b(s) = (/ &
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