      module cdata
      implicit none
      integer, parameter :: particles = 400 !set the number of particles
      integer, parameter :: shots = 10 !set how often to be written to file
      type main
      real :: x(3) !3 dimensional flow
      real :: p, b !angles phi and rho
      end type
      real :: t
      real,parameter :: h=0.1 !timestep
      type(main), dimension(particles) :: y
      end module
!*************************************************************************

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
!***************************************************************************
      program test
      use cdata
      use stats
      implicit none
      real :: q_11  !set values for using RK4 method
      real :: q_12  !*******************************
      real :: q_21
      real :: q_31
      real :: q_41
      real :: q_22
      real :: q_32
      real :: q_42
      real :: w_11
      real :: w_12
      real :: w_13
      real :: w_21
      real :: w_22
      real :: w_23
      real :: w_31
      real :: w_32
      real :: w_33
      real :: w_41
      real :: w_42
      real :: w_43  !set values for using RK4 method
!******************************************************************************
      real :: psi = 2.6 !measure of orientational stability
      real :: phi = 0.6 !swimming speed relative to the flow speed
      real :: a = 0.1 !set value of alpha
      real :: pi=3.141592653 !set the value pi for the program
      integer:: i, itime
      integer :: n

!initialise particle positions and swimming direction------
      do i=1, particles
      y(i)%x(1)=runif(-pi,pi)
      y(i)%x(2)=0
      y(i)%x(3)=runif(-pi,pi)
      y(i)%p=pi/2 !angle theta
      y(i)%b=pi/4 !angle phi
      end do

!-----------------------------------

      do itime=1,2000


      do i=1,particles


      q_11=(0.5)*(-(cos(y(i)%p))/(sin(y(i)%b))-
     &cos(y(i)%x(1))*cos(y(i)%x(3))+(a*sin(2*y(i)%p)
     &*sin(y(i)%x(1))*sin(y(i)%x(3))))
      q_12=(0.5)*(((1/psi)*(sin(y(i)%p))*cos(y(i)%b))
     &+(a*cos(y(i)%b)*sin(y(i)%x(1))*sin(y(i)%x(3))*
     &sin(y(i)%b)*cos(2*(y(i)%p))))
      w_11=(phi*cos(y(i)%p)*sin(y(i)%b))-
     &(0.5*cos(y(i)%x(1))*sin(y(i)%x(3)))
      w_12=phi*cos(y(i)%b)
      w_13=(phi*sin(y(i)%p)*sin(y(i)%b))+
     &(0.5*sin(y(i)%x(1))*cos(y(i)%x(3)))

      q_21=(0.5)*(-(cos(y(i)%p+(0.5*h*q_11)))/(sin(y(i)%b+
     &(0.5*h*q_12)))-cos(y(i)%x(1)+(0.5*h*w_11))*cos(y(i)%x(3)+
     &(0.5*h*w_13))+(a*sin(2*(y(i)%p+(0.5*h*q_11)))
     &*sin(y(i)%x(1)+(0.5*h*w_11))*sin(y(i)%x(3)+(0.5*h*w_13))))


      q_22=(0.5)*(((1/psi)*(sin(y(i)%p+(0.5*h*q_11)))*
     &cos(y(i)%b+(0.5*h*q_12)))+(a*cos(y(i)%b+(0.5*h*q_12))
     &*sin(y(i)%x(1)+(0.5*h*w_11))*sin(y(i)%x(3)+(0.5*h*w_13))
     &*sin(y(i)%b+(0.5*h*q_12))*cos(2*(y(i)%p+(0.5*h*q_11)))))


      w_21=(phi*cos(y(i)%p+(0.5*h*q_11))*sin(y(i)%b+
     &(0.5*h*q_12)))-(0.5*cos(y(i)%x(1)+(0.5*h*w_11))*
     &sin(y(i)%x(3)+(0.5*h*w_13)))

      w_22=phi*cos(y(i)%b+(0.5*h*q_12))

      w_23=(phi*sin(y(i)%p+(0.5*h*q_11))*sin(y(i)%b+
     &(0.5*h*q_12)))+(0.5*sin(y(i)%x(1)+(0.5*h*w_11))*
     &cos(y(i)%x(3)+(0.5*h*w_13)))


      q_31=(0.5)*(-(cos(y(i)%p+(0.5*h*q_21)))/(sin(y(i)%b+
     &(0.5*h*q_22)))-cos(y(i)%x(1)+(0.5*h*w_21))*cos(y(i)%x(3)+
     &(0.5*h*w_23))+(a*sin(2*(y(i)%p+(0.5*h*q_21)))
     &*sin(y(i)%x(1)+(0.5*h*w_21))*sin(y(i)%x(3)+(0.5*h*w_23))))

      q_32=(0.5)*(((1/psi)*(sin(y(i)%p+(0.5*h*q_21)))*
     &cos(y(i)%b+(0.5*h*q_22)))+(a*cos(y(i)%b+(0.5*h*q_22))
     &*sin(y(i)%x(1)+(0.5*h*w_21))*sin(y(i)%x(3)+(0.5*h*w_23))
     &*sin(y(i)%b+(0.5*h*q_22))*cos(2*(y(i)%p+(0.5*h*q_21)))))

      w_31=(phi*cos(y(i)%p+(0.5*h*q_21))*sin(y(i)%b+
     &(0.5*h*q_22)))-(0.5*cos(y(i)%x(1)+(0.5*h*w_21))*
     &sin(y(i)%x(3)+(0.5*h*w_23)))

      w_32=phi*cos(y(i)%b+(0.5*h*q_22))

      w_33=(phi*sin(y(i)%p+(0.5*h*q_21))*sin(y(i)%b+
     &(0.5*h*q_22)))+(0.5*sin(y(i)%x(1)+(0.5*h*w_21))*
     &cos(y(i)%x(3)+(0.5*h*w_23)))


      q_41=(0.5)*(-(cos(y(i)%p+(h*q_31)))/(sin(y(i)%b+
     &(h*q_32)))-cos(y(i)%x(1)+(h*w_31))*cos(y(i)%x(3)+
     &(h*w_33))+(a*sin(2*(y(i)%p+(h*q_31)))
     &*sin(y(i)%x(1)+(h*w_31))*sin(y(i)%x(3)+(h*w_33))))

      q_42=(0.5)*(((1/psi)*(sin(y(i)%p+(h*q_31)))*
     &cos(y(i)%b+(h*q_32)))+(a*cos(y(i)%b+(h*q_32))
     &*sin(y(i)%x(1)+(h*w_31))*sin(y(i)%x(3)+(h*w_33))
     &*sin(y(i)%b+(h*q_32))*cos(2*(y(i)%p+(h*q_31)))))

      w_41=(phi*cos(y(i)%p+(h*q_31))*sin(y(i)%b+
     &(h*q_32)))-(0.5*cos(y(i)%x(1)+(h*w_31))*
     &sin(y(i)%x(3)+(h*w_33)))

      w_42=phi*cos(y(i)%b+(h*q_32))

      w_43=(phi*sin(y(i)%p+(h*q_31))*sin(y(i)%b+
     &(h*q_32)))+(0.5*sin(y(i)%x(1)+(h*w_31))*
     &cos(y(i)%x(3)+(h*w_33)))





      y(i)%p=y(i)%p+((1.0/6.0)*(q_11+(2*q_21)+(2*q_31)+q_41))
      y(i)%b=y(i)%b+((1.0/6.0)*(q_12+(2*q_22)+(2*q_32)+q_42))
      y(i)%x(1)=y(i)%x(1)+((1.0/6.0)*(w_11+(2*w_21)+(2*w_31)+w_41))
      y(i)%x(2)=y(i)%x(2)+((1.0/6.0)*(w_12+(2*w_22)+(2*w_32)+w_42))
      y(i)%x(3)=y(i)%x(3)+((1.0/6.0)*(w_13+(2*w_23)+(2*w_33)+w_43))





















      do n=1,5000
      IF (y(i)%x(1)>pi) THEN         !ensure all solutions are kept
      y(i)%x(1)=y(i)%x(1)-(2*pi)     !in a box with dimensions pi
      ELSE IF (y(i)%x(1)<-pi) THEN   !in the x,y and z direction
      y(i)%x(1)=y(i)%x(1)+(2*pi)     !***************
      END IF

      IF (y(i)%x(2)>pi) THEN         !*********************
      y(i)%x(2)=y(i)%x(2)-(2*pi)
      ELSE IF (y(i)%x(2)<-pi) THEN
      y(i)%x(2)=y(i)%x(2)+(2*pi)
      END IF                         !**********************

      IF (y(i)%x(3)>pi) THEN         !**********************
      y(i)%x(3)=y(i)%x(3)-(2*pi)
      ELSE IF (y(i)%x(3)<-pi) THEN
      y(i)%x(3)=y(i)%x(3)+(2*pi)
      END IF                         !**********************

      end do




      end do

      if (mod(itime,shots)==0) then
      call output(itime/shots)
      end if

      end do



      end program test

      subroutine output(filenumber)
      use cdata
      implicit none
      integer, intent(IN) :: filenumber

      character (len=40) :: print_file
      write(unit=print_file,fmt="(a,i4.4,a)")"par",filenumber,".log"
      open(unit=98,file=print_file,status='replace',form='unformatted'
     &,access='stream')
      write(98) t
      write(98) particles
      write(98) y(:)%x(1)
      write(98) y(:)%x(2)
      write(98) y(:)%x(3)
      close(98)
      end subroutine