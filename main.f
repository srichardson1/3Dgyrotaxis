function u_1(x,y,z) !set value for the flow velocity in the x direction
      real :: u_1
      real :: x
      real :: y
      real :: z
      real :: a = -0.5
      u_1= a*cos(x)*sin(z)
      return
      end function

      function u_2(x,y,z)  !set value for the flow velocity in the y direction
      real :: u_2
      real :: x
      real :: y
      real :: z

      u_2=0
      return
      end function

      function u_3(x,y,z)  !set value for the flow velocity in the z direction
      real :: u_3
      real :: x
      real :: y
      real :: z
      real :: a=0.5
      u_3=a*sin(x)*cos(z)
      return
      end function


      function w_1(x,y,z) !set value for the vorticity in the x direction
      real :: w_1
      real :: x
      real ::y
      real :: z
      w_1=0
      return
      end function

      function w_2(x,y,z)  !set value for the vorticity in the y direction
      real :: w_2
      real :: x
      real :: y
      real :: z
      w_2=-cos(x)*cos(z)
      return
      end function

      function w_3(x,y,z)  !set value for the vorticity in the z direction
      real :: w_3
      real :: x
      real :: y
      real :: z
      w_3=0
      return
      end function



      function p_1(phi,theta)
      real :: p_1
      real :: theta
      real :: phi
      p_1=sin(theta)*cos(phi)
      return
      end function

      function p_2(phi,theta)
      real :: p_2
      real :: theta
      real :: phi
      p_2=sin(theta)*sin(phi)
      return
      end function

      function p_3(phi,theta)
      real :: p_3
      real :: theta
      real :: phi
      p_3=cos(theta)
      return
      end function

      function e_11(x,y,z) !set values for the rate of strain tensor
      real :: e_11
      real :: x
      real :: y
      real :: z
      e_11=0    !0.5*sin(x)*sin(z)
      return
      end function

      function e_12(x,y,z)  !set values for the rate of strain tensor
      real :: e_12
      real :: x
      real :: y
      real :: z
      e_12=0.5*(cos(x)-sin(y)) !0
      return
      end function

      function e_13(x,y,z)  !set values for the rate of strain tensor
      real :: e_13
      real :: x
      real ::y
      real :: z
      e_13=0.5*(cos(z)-sin(x)) !0
      return
      end function


      function e_22(x,y,z)  !set values for the rate of strain tensor
      real :: e_22
      real :: x
      real ::y
      real :: z
      e_22=0
      return
      end function

      function e_23(x,y,z)  !set values for the rate of strain tensor
      real :: e_23
      real :: x
      real :: y
      real :: z
      e_23=0.5*(cos(y)-sin(z))
      return
      end function

      function e_33(x,y,z)  !set values for the rate of strain tensor
      real :: e_33
      real :: x
      real :: y
      real :: z
      e_33=0    !-0.5*sin(x)*sin(z)
      return
      end function



!*******************************************************************************


      module cdata
      implicit none
      integer, parameter :: particles = 1000 !set the number of particles
      integer, parameter :: shots = 10 !set how often to be written to file
      type main
      real :: x(3) !3 dimensional flow
      real :: x_old(3)
      real :: phi !define angle phi
      real :: theta
      real :: velocity
      end type
      real :: t
      real,parameter :: h=0.1 !timestep used in the RK4 method
      real :: dt = 0.1 !timestep used
      type(main), dimension(particles) :: y
      end module
!******************************************************************************



      !>draw random variables from a \f$U(\alpha,\beta)\f$ distribution
      module stats
      contains
      real function runif(alpha,beta)
      implicit none
      real, intent(IN) :: alpha, beta
      real :: u
      call random_number(u)
      runif=alpha+u*(beta-alpha)
      end function
      end module

!************************************************************************

      program scott
      use omp_lib
      use cdata
      use stats
      implicit none
      integer :: i, itime
      real :: pi = 3.141592653 !set value for pi
      real :: u_1, u_2, u_3 !values for velocity
      real :: w_1, w_2, w_3 !values for vorticity
      real :: p_1, p_2, p_3 !values for swimming direction
      real :: e_11, e_12, e_13, e_22, e_23, e_33 !rate of strain tensor values
      real :: k_1, k_2, k_3, k_4 !RK4 for theta equation
      real :: l_1, l_2, l_3, l_4 !RK4 for phi equation
      real :: j_1, j_2, j_3, j_4 !RK4 for x equation
      real :: h_1, h_2, h_3, h_4 !RK4 for y equation
      real :: g_1, g_2, g_3, g_4 !RK4 for z equation
      real :: G = 2 !measure or orientational stability
      real :: V = 0.7 !swimming speed relative to the flow speed
      integer :: alpha = 0. !determines the shape of the cells
      integer :: n
      real :: velocity
!*************************************************************************

      !initialise particle positions and swimming direction------
      do i=1, particles
      y(i)%x(1)=runif(-pi,pi)
      y(i)%x(2)=runif(-pi,pi)
      y(i)%x(3)=runif(-pi,pi)
      y(i)%phi=pi/2
      y(i)%theta=pi/2
      end do

!**************************************************************************



!$OMP PARALLEL
!$OMP DO

      do itime=1,1000
      t=itime*dt

!$OMP DO

      do i=1,particles






      k_1=(-1/(sin(y(i)%theta)))&
      *((((1/(2.0*G))*(1-(p_3(y(i)%phi,y(i)%theta))**2)))&
      +(0.5*((w_1(y(i)%x(1),y(i)%x(2),y(i)%x(3))*p_2(y(i)%phi,&
      y(i)%theta))-((w_2(y(i)%x(1),y(i)%x(2),y(i)%x(3))*p_1&
      (y(i)%phi,y(i)%theta)))))+(alpha*(((-e_11(y(i)%x(1),&
      y(i)%x(2),y(i)%x(3)))*((p_1(y(i)%phi,y(i)%theta))**2)*&
      (p_3(y(i)%phi,y(i)%theta))-(2*(e_12(y(i)%x(1),y(i)%x(2),&
      y(i)%x(3))*(p_1(y(i)%phi,y(i)%theta))*(p_2(y(i)%phi,y(i)&
      %theta))*(p_3(y(i)%phi,y(i)%theta))))+((e_13(y(i)%x(1),&
      y(i)%x(2),y(i)%x(3)))*((p_1(y(i)%phi,y(i)%theta))*&
      (1-(2*(p_3(y(i)%phi,y(i)%theta))**2))))-((e_22(y(i)%x(1)&
      ,y(i)%x(2),y(i)%x(3)))*(p_3(y(i)%phi,y(i)%theta))*&
      ((p_2(y(i)%phi,y(i)%theta))**2))+((e_23(y(i)%x(1),y(i)%x(2),&
      y(i)%x(3)))*((p_2(y(i)%phi,y(i)%theta))*(1-(2*&
      ((p_3(y(i)%phi,y(i)%theta))**2)))))+((e_33(y(i)%x(1),&
      y(i)%x(2),y(i)%x(3)))*((p_3(y(i)%phi,y(i)%theta))*&
      (1-((p_3(y(i)%phi,y(i)%theta))**2))))))))




      l_1=(1/sin(y(i)%theta))*(0.5*(((w_3(y(i)%x(1),y(i)%x(2),&
      y(i)%x(3)))*((p_1(y(i)%phi,y(i)%theta))**2+((p_2(y(i)%phi,&
      y(i)%theta))**2)))-(w_1(y(i)%x(1),y(i)%x(2),y(i)%x(3))*p_3&
      (y(i)%phi,y(i)%theta)*p_1(y(i)%phi,y(i)%theta))-&
      (w_2(y(i)%x(1),y(i)%x(2),y(i)%x(3))*p_2(y(i)%phi,y(i)%theta)&
      *p_3(y(i)%phi,y(i)%theta)))+(alpha*((((-e_11(y(i)%x(1),y(i)%x(2),&
      y(i)%x(3)))*(p_1(y(i)%phi,y(i)%theta))*(p_2(y(i)%phi,&
      y(i)%theta)))+((e_12(y(i)%x(1),y(i)%x(2),y(i)%x(3)))*&
      (((p_1(y(i)%phi,y(i)%theta))**2)-((p_2(y(i)%phi,y(i)%theta))&
      **2)))-(e_13(y(i)%x(1),y(i)%x(2),y(i)%x(3))*(p_3(y(i)%phi,&
      y(i)%theta))*(p_2(y(i)%phi,y(i)%theta)))+((e_22(y(i)%x(1)&
      ,y(i)%x(2),y(i)%x(3)))*(p_2(y(i)%phi,y(i)%theta))*&
      (p_1(y(i)%phi,y(i)%theta)))+((e_23(y(i)%x(1),y(i)%x(2),&
      y(i)%x(3)))*(p_1(y(i)%phi,y(i)%theta))*(p_3(y(i)%phi,&
      y(i)%theta)))))))

      j_1=(V*(p_1(y(i)%phi,y(i)%theta)))+&
      u_1(y(i)%x(1),y(i)%x(2),y(i)%x(3))

      h_1=(V*(p_2(y(i)%phi,y(i)%theta)))+&
      u_2(y(i)%x(1),y(i)%x(2),y(i)%x(3))

      g_1=(V*(p_3(y(i)%phi,y(i)%theta)))+&
      u_3(y(i)%x(1),y(i)%x(2),y(i)%x(3))




      k_2=(-1/(sin(y(i)%theta+(0.5*h*k_1))))&
      *((((1/(2.0*G))*(1-(p_3(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))**2)))+(0.5*((w_1(y(i)%x(1)+&
      (0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1))&
      *p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))-&
      ((w_2(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1))*p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))))&
      +(alpha*(((-e_11(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),&
      y(i)%x(3)+(0.5*h*g_1)))*((p_1(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))**2)*(p_3(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))-(2*(e_12(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)&
      +(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1))*(p_1(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))*(p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+&
      (0.5*h*k_1)))*(p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))))&
      +((e_13(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1)))*((p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))&
      *(1-(2*(p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2))))&
      -((e_22(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)&
      +(0.5*h*g_1)))*(p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))&
      *((p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2))+&
      ((e_23(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1)))*((p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))&
      *(1-(2*((p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2)))))&
      +((e_33(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1)))*((p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))*&
      (1-((p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2))))))))


      l_2=(1/sin(y(i)%theta)+(0.5*h*k_1))*(0.5*(((w_3(y(i)%x(1)+&
      (0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1)))*&
      ((p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2+&
      ((p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))**2)))-&
      (w_1(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1))*p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))&
      *p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))-&
      (w_2(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+&
      (0.5*h*g_1))*p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))&
      *p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))+(alpha*&
      ((((-e_11(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1),&
      y(i)%x(3)+(0.5*h*g_1)))*(p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+&
      (0.5*h*k_1)))*(p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)&
      )))+((e_12(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1)&
      ,y(i)%x(3)+(0.5*h*g_1)))*(((p_1(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))**2)-((p_2(y(i)%phi+(0.5*h*l_1),&
      y(i)%theta+(0.5*h*k_1)))**2)))-(e_13(y(i)%x(1)+(0.5*h*j_1),&
      y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1))*(p_3(y(i)%phi+&
      (0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))*(p_2(y(i)%phi+(0.5*h*l_1)&
      ,y(i)%theta+(0.5*h*k_1))))+((e_22(y(i)%x(1)+(0.5*h*j_1),&
      y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1)))*(p_2(y(i)%phi&
      +(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))*(p_1(y(i)%phi+(0.5*h*l_1)&
      ,y(i)%theta+(0.5*h*k_1))))+((e_23(y(i)%x(1)+(0.5*h*j_1)&
      ,y(i)%x(2)+(0.5*h*h_1),y(i)%x(3)+(0.5*h*g_1)))*(p_1(y(i)%phi&
      +(0.5*h*l_1),y(i)%theta+(0.5*h*k_1)))*(p_3(y(i)%phi+&
      (0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))))))


      j_2=(V*(p_1(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))+&
      u_1(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1)&
      ,y(i)%x(3)+(0.5*h*g_1))



      h_2=(V*(p_2(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))+&
      u_2(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1)&
      ,y(i)%x(3)+(0.5*h*g_1))



      g_2=(V*(p_3(y(i)%phi+(0.5*h*l_1),y(i)%theta+(0.5*h*k_1))))+&
      u_3(y(i)%x(1)+(0.5*h*j_1),y(i)%x(2)+(0.5*h*h_1)&
      ,y(i)%x(3)+(0.5*h*g_1))


      k_3=(-1/(sin(y(i)%theta+(0.5*h*k_2))))&
      *((((1/(2.0*G))*(1-(p_3(y(i)%phi+(0.5*h*l_2),&
      y(i)%theta+(0.5*h*k_2)))**2)))+(0.5*((w_1(y(i)%x(1)+&
      (0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2))&
      *p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))-&
      ((w_2(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+&
      (0.5*h*g_2))*p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))))))&
      +(alpha*(((-e_11(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),&
      y(i)%x(3)+(0.5*h*g_2)))*((p_1(y(i)%phi+(0.5*h*l_2),&
      y(i)%theta+(0.5*h*k_2)))**2)*(p_3(y(i)%phi+(0.5*h*l_2),&
      y(i)%theta+(0.5*h*k_2)))-(2*(e_12(y(i)%x(1)+(0.5*h*j_2),&
      y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2))*(p_1(y(i)%phi+&
      (0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))*(p_2(y(i)%phi+(0.5*h*l_2),&
      y(i)%theta+(0.5*h*k_2)))*(p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+&
      (0.5*h*k_2)))))+((e_13(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+&
      (0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*((p_1(y(i)%phi+(0.5*h*l_2)&
      ,y(i)%theta+(0.5*h*k_2)))*(1-(2*(p_3(y(i)%phi+(0.5*h*l_2),y(i)&
      %theta+(0.5*h*k_2)))**2))))-((e_22(y(i)%x(1)+(0.5*h*j_2),&
      y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*(p_3(y(i)%phi+&
      (0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))*((p_2(y(i)%phi+(0.5*h*l_2)&
      ,y(i)%theta+(0.5*h*k_2)))**2))+((e_23(y(i)%x(1)+(0.5*h*j_2),&
      y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*((p_2(y(i)%phi+&
      (0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))*(1-(2*((p_3(y(i)%phi+&
      (0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))**2)))))+((e_33(y(i)%x(1)+&
      (0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*&
      ((p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))*&
      (1-((p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))**2))))))))




      l_3=(1/sin(y(i)%theta)+(0.5*h*k_2))*(0.5*(((w_3(y(i)%x(1)+&
      (0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*&
      ((p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))**2+&
      ((p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))**2)))-&
      (w_1(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+&
      (0.5*h*g_2))*p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))&
      *p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2)))-&
      (w_2(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),y(i)%x(3)+&
      (0.5*h*g_2))*p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))&
      *p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))))+(alpha*&
      ((((-e_11(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),&
      y(i)%x(3)+(0.5*h*g_2)))*(p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+&
      (0.5*h*k_2)))*(p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2&
      ))))+((e_12(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2),&
      y(i)%x(3)+(0.5*h*g_2)))*(((p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+&
      (0.5*h*k_2)))**2)-((p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+&
      (0.5*h*k_2)))**2)))-(e_13(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+&
      (0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2))*(p_3(y(i)%phi+(0.5*h*l_2)&
      ,y(i)%theta+(0.5*h*k_2)))*(p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta&
      +(0.5*h*k_2))))+((e_22(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+&
      (0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*(p_2(y(i)%phi+(0.5*h*l_2)&
      ,y(i)%theta+(0.5*h*k_2)))*(p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta&
      +(0.5*h*k_2))))+((e_23(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+&
      (0.5*h*h_2),y(i)%x(3)+(0.5*h*g_2)))*(p_1(y(i)%phi+(0.5*h*l_2)&
      ,y(i)%theta+(0.5*h*k_2)))*(p_3(y(i)%phi+(0.5*h*l_2),&
      y(i)%theta+(0.5*h*k_2))))))))




      j_3=(V*(p_1(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))))+&
      u_1(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2)&
      ,y(i)%x(3)+(0.5*h*g_2))




      h_3=(V*(p_2(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))))+&
      u_2(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2)&
      ,y(i)%x(3)+(0.5*h*g_2))




      g_3=(V*(p_3(y(i)%phi+(0.5*h*l_2),y(i)%theta+(0.5*h*k_2))))+&
      u_3(y(i)%x(1)+(0.5*h*j_2),y(i)%x(2)+(0.5*h*h_2)&
      ,y(i)%x(3)+(0.5*h*g_2))


      k_4=(-1/(sin(y(i)%theta+(h*k_3))))&
      *((((1/(2.0*G))*(1-(p_3(y(i)%phi+(h*l_3),&
      y(i)%theta+(h*k_3)))**2)))+(0.5*((w_1(y(i)%x(1)+&
      (h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3))&
      *p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))-&
      ((w_2(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+&
      (h*g_3))*p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))))+&
      (alpha*(((-e_11(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),&
      y(i)%x(3)+(h*g_3)))*((p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))&
      **2)*(p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))-(2*(e_12&
      (y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3))*&
      (p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))*(p_2(y(i)%phi+&
      (h*l_3),y(i)%theta+(h*k_3)))*(p_3(y(i)%phi+(h*l_3),y(i)%theta+&
      (h*k_3)))))+((e_13(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)&
      +(h*g_3)))*((p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))*&
      (1-(2*(p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))**2))))-&
      ((e_22(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3)))&
      *(p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))*((p_2(y(i)%phi+&
      (h*l_3),y(i)%theta+(h*k_3)))**2))+((e_23(y(i)%x(1)+(h*j_3),&
      y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3)))*((p_2(y(i)%phi+(h*l_3)&
      ,y(i)%theta+(h*k_3)))*(1-(2*((p_3(y(i)%phi+(h*l_3),y(i)%theta&
      +(h*k_3)))**2)))))+((e_33(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3)&
      ,y(i)%x(3)+(h*g_3)))*((p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))&
      *(1-((p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))**2))))))))



      l_4=(1/sin(y(i)%theta)+(h*k_3))*(0.5*(((w_3(y(i)%x(1)+&
      (h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3)))*&
      ((p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))**2+&
      ((p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))**2)))-&
      (w_1(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+&
      (h*g_3))*p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))&
      *p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))-&
      (w_2(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+&
      (h*g_3))*p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))&
      *p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+&
      (alpha*((((-e_11(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),&
      y(i)%x(3)+(h*g_3)))*(p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))&
      *(p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+((e_12(y(i)%x(1)+&
      (h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3)))*(((p_1(y(i)%phi+&
      (h*l_3),y(i)%theta+(h*k_3)))**2)-((p_2(y(i)%phi+(h*l_3),&
      y(i)%theta+(h*k_3)))**2)))-(e_13(y(i)%x(1)+(h*j_3),y(i)%x(2)+&
      (h*h_3),y(i)%x(3)+(h*g_3))*(p_3(y(i)%phi+(h*l_3),y(i)%theta&
      +(h*k_3)))*(p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+&
      ((e_22(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3),y(i)%x(3)+(h*g_3)))&
      *(p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))*(p_1(y(i)%phi+(h*l_3)&
      ,y(i)%theta+(h*k_3))))+((e_23(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3)&
      ,y(i)%x(3)+(h*g_3)))*(p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3)))&
      *(p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))))))


      j_4=(V*(p_1(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+&
      u_1(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3)&
      ,y(i)%x(3)+(h*g_3))


      h_4=(V*(p_2(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+&
      u_2(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3)&
      ,y(i)%x(3)+(h*g_3))






      g_4=(V*(p_3(y(i)%phi+(h*l_3),y(i)%theta+(h*k_3))))+&
      u_3(y(i)%x(1)+(h*j_3),y(i)%x(2)+(h*h_3)&
      ,y(i)%x(3)+(h*g_3))



      y(i)%x_old(2)=y(i)%x(2)


      y(i)%phi=y(i)%phi+((1.0/6.0)*(l_1+(2*l_2)+(2*l_3)+l_4))
      y(i)%theta=y(i)%theta+((1.0/6.0)*(k_1+(2*k_2)+(2*k_3)+k_4))
      y(i)%x(1)=y(i)%x(1)+((1.0/6.0)*(j_1+(2*j_2)+(2*j_3)+j_4))
      y(i)%x(2)=y(i)%x(2)+((1.0/6.0)*(h_1+(2*h_2)+(2*h_3)+h_4))
      y(i)%x(3)=y(i)%x(3)+((1.0/6.0)*(g_1+(2*g_2)+(2*g_3)+g_4))

      y(i)%velocity=(y(i)%x(2)-y(i)%x_old(2))/dt

!*************************************************************************
!$OMP DO
      do n=1,5000
      IF (y(i)%x(1)>pi) THEN         !ensure all solutions are kept
      y(i)%x(1)=y(i)%x(1)-(2*pi)     !in a box with dimensions pi
      ELSE IF (y(i)%x(1)<-pi) THEN   !in the x,y and z direction
      y(i)%x(1)=y(i)%x(1)+(2*pi)
      END IF

      IF (y(i)%x(2)>pi) THEN
      y(i)%x(2)=y(i)%x(2)-(2*pi)
      ELSE IF (y(i)%x(2)<-pi) THEN
      y(i)%x(2)=y(i)%x(2)+(2*pi)
      END IF

      IF (y(i)%x(3)>pi) THEN
      y(i)%x(3)=y(i)%x(3)-(2*pi)
      ELSE IF (y(i)%x(3)<-pi) THEN
      y(i)%x(3)=y(i)%x(3)+(2*pi)
      END IF
      end do
!$OMP END DO
!**************************************************************************



      end do

!$OMP END DO


      if (mod(itime,shots)==0) then
      call output(itime/shots)
      end if

      end do
!$OMP END DO
!$OMP END PARALLEL


      end program


!***************************************************************************

      subroutine output(filenumber)
      use cdata
      implicit none
      integer, intent(IN) :: filenumber

      character (len=40) :: print_file
      write(unit=print_file,fmt="(a,i4.4,a)")"par",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted'&
      ,access='stream')
      write(98) t
      write(98) particles
      write(98) y(:)%x(1)
      write(98) y(:)%x(2)
      write(98) y(:)%x(3)
      write(98) y(:)%velocity
      close(98)
      end subroutine





