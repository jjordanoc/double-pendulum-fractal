module example_plot_1
   use, intrinsic :: iso_fortran_env
   use fplot_core
!    use forcolormap, only : colormaps_list
   implicit none
contains
   subroutine plot1
      ! Parameters
      integer(int32), parameter :: n = 1000

      ! Local Variables
      real(real64), dimension(n) :: x, y1, y2
      type(plot_2d) :: plt
      type(plot_data_2d) :: d1, d2
      class(plot_axis), pointer :: xAxis, yAxis
      type(legend), pointer :: leg

      ! Initialize the plot object
      call plt%initialize(term=GNUPLOT_TERMINAL_QT)

      ! Define titles
      call plt%set_title("Example Plot")
      call plt%set_font_size(14)

      xAxis => plt%get_x_axis()
      call xAxis%set_title("X Axis")

      yAxis => plt%get_y_axis()
      call yAxis%set_title("Y Axis")

      ! Establish legend properties
      leg => plt%get_legend()
      call leg%set_is_visible(.true.)
      call leg%set_draw_inside_axes(.false.)
      call leg%set_horizontal_position(LEGEND_CENTER)
      call leg%set_vertical_position(LEGEND_BOTTOM)
      call leg%set_draw_border(.false.)

      ! Define the data, and then add it to the plot
      x = linspace(0.0d0, 10.0d0, n)
      y1 = sin(5.0d0 * x)
      y2 = 2.0d0 * cos(2.0d0 * x)

      call d1%define_data(x, y1)
      call d2%define_data(x, y2)

      ! Define properties for each data set
      call d1%set_name("Data Set 1")
      call d1%set_draw_markers(.true.)
      call d1%set_marker_frequency(10)
      call d1%set_marker_style(MARKER_EMPTY_CIRCLE)
      call d1%set_marker_scaling(2.0)

      call d2%set_name("Data Set 2")
      call d2%set_line_style(LINE_DASHED)
      call d2%set_line_width(2.0)

      ! Add the data sets to the plot
      call plt%push(d1)
      call plt%push(d2)

      ! Let GNUPLOT draw the plot
      call plt%draw()
   end subroutine plot1

   subroutine plot2
      ! Local Variables
      type(plot_2d) :: plt
      type(vector_field_plot_data) :: ds1
      class(plot_axis), pointer :: xAxis, yAxis
      type(rainbow_colormap) :: cmap
      real(real64), allocatable, dimension(:,:,:) :: pts
      real(real64), allocatable, dimension(:,:) :: dx, dy
      real(real64) :: dxdt(2)
      integer(int32) :: i, j

      ! Create a grid of points defining the vector locations
      pts = meshgrid( &
         linspace(-2.0d0, 2.0d0, 20), &
         linspace(-5.0d0, 5.0d0, 20))

      ! Compute the values of each derivative
      allocate(dx(size(pts, 1), size(pts, 2)))
      allocate(dy(size(pts, 1), size(pts, 2)))
      do j = 1, size(pts, 2)
         do i = 1, size(pts, 1)
            call eqn([pts(i,j,1), pts(i,j,2)], dxdt)
            dx(i,j) = dxdt(1)
            dy(i,j) = dxdt(2)
         end do
      end do

      ! Define arrow properties
      call ds1%set_arrow_size(0.0001d0)  ! 1.0 by default
      call ds1%set_fill_arrow(.true.) ! .false. by default

      ! Create the plot
      call plt%initialize(term=GNUPLOT_TERMINAL_QT)
      call plt%set_font_size(14)
      xAxis => plt%get_x_axis()
      yAxis => plt%get_y_axis()

      ! Define axis labels
      call xAxis%set_title("x(t)")
      call yAxis%set_title("dx/dt")

      ! Set plot style information
      call xAxis%set_zero_axis(.true.)
      call yAxis%set_zero_axis(.true.)
      call plt%set_draw_border(.false.)
      call plt%set_show_gridlines(.false.)

      ! Define the colormap
    !   call plt%set_colormap(cmap)

      ! Add the data to the plot - color by the magnitude of gradient
      call ds1%define_data(pts(:,:,1), pts(:,:,2), dx, dy, sqrt(dx**2 + dy**2))
      call plt%push(ds1)

      call plt%draw()
   contains
      ! Van der Pol Equation
      ! x" - mu * (1 - x^2) * x' + x = 0
      subroutine eqn(x, dxdt)
         real(real64), intent(in) :: x(2)
         real(real64), intent(out) :: dxdt(2)

         real(real64), parameter :: mu = 2.0d0

         dxdt(1) = x(2)
         dxdt(2) = mu * (1.0d0 - x(1)**2) * x(2) - x(1)
      end subroutine

   end subroutine plot2
end module example_plot_1
