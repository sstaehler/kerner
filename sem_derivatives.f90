!=========================================================================================
module sem_derivatives
  use global_parameters,      only : dp
  use commpi,                 only : pabort
  use finite_elem_mapping,    only : inv_jacobian

  implicit none
  private

  public :: axisym_gradient
  public :: dsdf_axis

  interface dsdf_axis
    module procedure  :: dsdf_axis
    module procedure  :: dsdf_axis_td
  end interface

  interface axisym_gradient
    module procedure  :: axisym_gradient
    module procedure  :: axisym_gradient_td
  end interface

  interface mxm
    module procedure  :: mxm
    module procedure  :: mxm_atd
    module procedure  :: mxm_btd
  end interface

  interface mxm_ipol0
    module procedure  :: mxm_ipol0
    module procedure  :: mxm_ipol0_atd
    module procedure  :: mxm_ipol0_btd
  end interface

contains

!-----------------------------------------------------------------------------------------
function strain_monopole(u, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the strain tensor for displacement u excited bz a monopole source
  ! in Voigt notation: [dsus, dpup, dzuz, dzup, dsuz, dsup]
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: u(0:npol,0:npol, 3)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: strain_monopole(0:npol,0:npol,6)
  
  real(kind=dp)                 :: grad_buff1(0:npol,0:npol,2)
  real(kind=dp)                 :: grad_buff2(0:npol,0:npol,2)
  
  ! 1: dsus, 2: dzus
  grad_buff1 = axisym_gradient(u(:,:,1), G, GT, xi, eta, npol, nodes, element_type)
  
  ! 1: dsuz, 2: dzuz
  grad_buff2 = axisym_gradient(u(:,:,3), G, GT, xi, eta, npol, nodes, element_type)


  strain_monopole(:,:,1) = grad_buff1(:,:,1)
  strain_monopole(:,:,2) = 0 !@TODO implement f/s
  strain_monopole(:,:,3) = grad_buff2(:,:,2)
  strain_monopole(:,:,4) = 0
  strain_monopole(:,:,5) = grad_buff1(:,:,2) + grad_buff2(:,:,1)
  strain_monopole(:,:,6) = 0

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function dsdf_axis_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(1:nsamp,0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for 
                                                     ! axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial 
                                               ! elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: dsdf_axis_td(1:nsamp,0:npol)

  real(kind=dp)                 :: inv_j_npol(0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm_ipol0_1(1:nsamp,0:npol)
  real(kind=dp)                 :: mxm_ipol0_2(1:nsamp,0:npol)

  ipol = 0
  do jpol = 0, npol
     inv_j_npol(jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
  enddo

  mxm_ipol0_1 = mxm_ipol0(GT,f)
  mxm_ipol0_2 = mxm_ipol0(f,G)

  do jpol = 0, npol
     dsdf_axis_td(:,jpol) =   inv_j_npol(jpol,1,1) * mxm_ipol0_1(:,jpol) &
                            + inv_j_npol(jpol,2,1) * mxm_ipol0_2(:,jpol)
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function dsdf_axis(f, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the partial derivative of scalar field f for ipol = 0
  ! needed for l'hospitals rule to compute f/s = df/ds at the axis s = 0
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: dsdf_axis(0:npol)

  real(kind=dp)                 :: inv_j_npol(0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm_ipol0_1(0:npol)
  real(kind=dp)                 :: mxm_ipol0_2(0:npol)

  ipol = 0
  do jpol = 0, npol
     inv_j_npol(jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
  enddo

  mxm_ipol0_1 = mxm_ipol0(GT,f)
  mxm_ipol0_2 = mxm_ipol0(f,G)
  dsdf_axis(:) =   inv_j_npol(:,1,1) * mxm_ipol0_1(:) &
                 + inv_j_npol(:,2,1) * mxm_ipol0_2(:)

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function axisym_gradient_td(f, G, GT, xi, eta, npol, nsamp, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol, nsamp
  real(kind=dp), intent(in)     :: f(1:nsamp,0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient_td(1:nsamp,0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(1:nsamp,0:npol,0:npol)
  real(kind=dp)                 :: mxm2(1:nsamp,0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)

  do jpol = 0, npol
     do ipol = 0, npol
        axisym_gradient_td(:,ipol,jpol,1) =   &
                inv_j_npol(ipol,jpol,1,1) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,1) * mxm2(:,ipol,jpol)
        axisym_gradient_td(:,ipol,jpol,2) =   &
                inv_j_npol(ipol,jpol,1,2) * mxm1(:,ipol,jpol) &
              + inv_j_npol(ipol,jpol,2,2) * mxm2(:,ipol,jpol)
     enddo
  enddo

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function axisym_gradient(f, G, GT, xi, eta, npol, nodes, element_type)
  ! Computes the axisymmetric gradient of scalar field f
  ! grad = \nabla {f} = \partial_s(f) \hat{s} + \partial_z(f) \hat{z}
  
  integer, intent(in)           :: npol
  real(kind=dp), intent(in)     :: f(0:npol,0:npol)
  real(kind=dp), intent(in)     :: G(0:npol,0:npol)  ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: GT(0:npol,0:npol) ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: xi(0:npol)  ! GLL for non-axial and GLJ for axial elements
  real(kind=dp), intent(in)     :: eta(0:npol) ! same for all elements (GLL)
  real(kind=dp), intent(in)     :: nodes(4,2)
  integer, intent(in)           :: element_type
  real(kind=dp)                 :: axisym_gradient(0:npol,0:npol,1:2)

  real(kind=dp)                 :: inv_j_npol(0:npol,0:npol,2,2)
  integer                       :: ipol, jpol
  real(kind=dp)                 :: mxm1(0:npol,0:npol)
  real(kind=dp)                 :: mxm2(0:npol,0:npol)

  do ipol = 0, npol
     do jpol = 0, npol
        inv_j_npol(ipol,jpol,:,:) = inv_jacobian(xi(ipol), eta(jpol), nodes, element_type)
     enddo
  enddo

!        | dxi  / ds  dxi  / dz |
! J^-1 = |                      |
!        | deta / ds  deta / dz |

  mxm1 = mxm(GT,f)
  mxm2 = mxm(f,G)
  axisym_gradient(:,:,1) =   inv_j_npol(:,:,1,1) * mxm1 &
                           + inv_j_npol(:,:,2,1) * mxm2
  axisym_gradient(:,:,2) =   inv_j_npol(:,:,1,2) * mxm1 &
                           + inv_j_npol(:,:,2,2) * mxm2

end function
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! a is time dependent
pure function mxm_atd(a, b)

  real(kind=dp), intent(in)  :: a(1:,0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_atd(1:size(a,1), 0:size(a,2)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_atd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,2) -1
        do k = 0, size(a,3) -1
           mxm_atd(:,i,j) = mxm_atd(:,i,j) + a(:,i,k) * b(k,j)
        enddo
     end do
  end do 

end function mxm_atd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
!! b is time dependent
pure function mxm_btd(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(1:,0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_btd(1:size(b,1),0:size(a,1)-1,0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_btd = 0

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        do k = 0, size(a,2) -1
           mxm_btd(:,i,j) = mxm_btd(:,i,j) + a(i,k) * b(:,k,j)
        enddo
     end do
  end do 

end function mxm_btd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c.
pure function mxm(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm(0:size(a,1)-1,0:size(b,2)-1)    !< Result
  integer                    :: i, j

  do j = 0, size(b,2) -1
     do i = 0, size(a,1) -1
        mxm(i,j) = sum(a(i,:) * b(:,j))
     end do
  end do 

end function mxm
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
!! a is time dependent
pure function mxm_ipol0_atd(a, b)

  real(kind=dp), intent(in)  :: a(1:,0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0_atd(1:size(a,1), 0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_ipol0_atd = 0
  i = 0

  do j = 0, size(b,2) -1
     do k = 0, size(a,3) -1
        mxm_ipol0_atd(:,j) = mxm_ipol0_atd(:,j) + a(:,i,k) * b(k,j)
     enddo
  end do 

end function mxm_ipol0_atd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
!! b is time dependent
pure function mxm_ipol0_btd(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(1:,0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0_btd(1:size(b,1),0:size(b,2)-1)  !< Result
  integer                    :: i, j, k

  mxm_ipol0_btd = 0

  i = 0
  do j = 0, size(b,2) -1
     do k = 0, size(a,2) -1
        mxm_ipol0_btd(:,j) = mxm_ipol0_btd(:,j) + a(i,k) * b(:,k,j)
     enddo
  end do 

end function mxm_ipol0_btd
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!> Multiplies matrizes a and b to have c, but only computes the component (0,:) of c
pure function mxm_ipol0(a, b)

  real(kind=dp), intent(in)  :: a(0:,0:), b(0:,0:)                  !< Input matrices
  real(kind=dp)              :: mxm_ipol0(0:size(b,2)-1)    !< Result
  integer                    :: i, j

  i = 0
  do j = 0, size(b,2) -1
     mxm_ipol0(j) = sum(a(i,:) * b(:,j))
  end do 

end function mxm_ipol0
!-----------------------------------------------------------------------------------------

end module
!=========================================================================================
