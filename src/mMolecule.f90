!   This file is part of Playmol.
!
!    Playmol is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Playmol is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Playmol. If not, see <http://www.gnu.org/licenses/>.
!
!    Author: Charlles R. A. Abreu (abreu at eq.ufrj.br)
!            Applied Thermodynamics and Molecular Simulation
!            Federal University of Rio de Janeiro, Brazil

module mMolecule

use mGlobal
use mStruc

implicit none

type tMolecule
  integer :: N = 0
  type(StrucList) :: list = StrucList( "molecule" )
  type(StrucList) :: bonds = StrucList( "bond", 2, .true. )
  type(StrucList) :: xyz = StrucList( "coordinate" )
  contains
    procedure :: add_atom => tMolecule_add_atom
    procedure :: index => tMolecule_index
    procedure :: fuse => tMolecule_fuse
    procedure :: split => tMolecule_split
    procedure :: number_of_atoms => tMolecule_number_of_atoms
    procedure :: count => tMolecule_count
    procedure :: per_molecule => tMolecule_per_molecule
    procedure :: set_geometry => tMolecule_set_geometry
    procedure :: coordinates => tMolecule_coordinates
    procedure, nopass :: align => tMolecule_align
end type tMolecule

contains

  !=================================================================================================

  subroutine tMolecule_add_atom( me, atom )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(in)    :: atom
    character(sl) :: arg(2)
    me%N = me%N + 1
    arg = [atom,int2str(me%N)]
    call me % list % add( 2, arg )
  end subroutine tMolecule_add_atom

  !=================================================================================================

  function tMolecule_index( me, string ) result( imol )
    class(tMolecule), intent(in) :: me
    character(sl),    intent(in) :: string
    integer                      :: imol
    integer       :: last
    character(sl) :: atom
    type(Struc), pointer :: ptr
    last = len_trim(string)
    if ((last > 5).and.(string(1:4) == "mol(").and.(string(last:last) == ")")) then
      atom = string(5:last-1)
      if (has_macros(atom)) call error( "invalid atom name", atom )
      call me % list % search( [atom], ptr )
      if (.not.associated(ptr)) call error( "atom", atom, "does not exist" )
      imol = str2int( ptr % params )
    else if (is_int(string)) then
      imol = str2int(string)
      if (imol < 1) call error( "molecule index cannot be zero or negative" )
      if (imol > me%N) call error( "molecule index must lie between 1 and ", int2str(me%N) )
    else
      call error( "molecule definition as", string, "is invalid" )
    end if
  end function tMolecule_index

  !=================================================================================================

  subroutine tMolecule_fuse( me, atoms )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(inout) :: atoms(2)
    integer :: i, mol(2), imin, imax
    type(Struc), pointer :: ptr
    if (atoms(1) == atoms(2)) call error( "atom", atoms(1), "cannot bind to itself" )
    do i = 1, 2
      call me % list % search( atoms(i:i), ptr )
      mol(i) = str2int( ptr % params )
    end do
    imin = minval(mol)
    imax = maxval(mol)
    if (imin == imax) then
      call writeln( "Closing cycle(s) in molecule", int2str(imin) )
    else
      call writeln( "Fusing molecule", int2str(imax), "to molecule", int2str(imin) )
      call rename_molecule( trim(int2str(imax)), trim(int2str(imin)) )
      if (imax < me%N) then
        call writeln( "Lowering indices of molecules", int2str(imax+1), "to", int2str(me%N) )
        do i = imax+1, me%N
          call rename_molecule( trim(int2str(i)), trim(int2str(i-1)) )
        end do
      end if
      me%N = me%N - 1
    end if
    call me % bonds % add( 2, atoms, silent = .true. )
    contains
      !---------------------------------------------------------------------------------------------
      subroutine rename_molecule( old, new )
        character(*), intent(in) :: old, new
        ptr => me % list % first
        do while (associated(ptr))
          if (ptr % params == old) ptr % params = new
          ptr => ptr % next
        end do
      end subroutine rename_molecule
      !---------------------------------------------------------------------------------------------
  end subroutine tMolecule_fuse

  !=================================================================================================

  subroutine tMolecule_split( me, arg )
    class(tMolecule), intent(inout) :: me
    character(sl),    intent(in)    :: arg(2)
    integer :: i, j, k, imin, imax, natoms
    character(sl) :: cmol, newmol
    type(Struc), pointer :: ptr, pcoords
    integer, allocatable :: mol(:), sort(:)
    character(sl), allocatable :: atom(:), id(:), params(:)
    ! Remove the bond between these atoms (if any) from the list:
    call me % bonds % remove( arg, silent = .true. )
    ! Determine the molecule two which these atoms belong:
    call me % list % search( arg(1:1), ptr )
    cmol = me % list % parameters( arg(1:1) )
    ! Determine the number of atoms of this molecule:
    natoms = 0
    ptr => me % list % first
    do while (associated(ptr))
      if (ptr%params == cmol) natoms = natoms + 1
      ptr => ptr % next
    end do
    ! Save the atoms as if they were monoatomic molecules:
    allocate( mol(natoms), atom(natoms) )
    ptr => me % list % first
    do i = 1, natoms
      mol(i) = i
      do while (ptr%params /= cmol)
        ptr => ptr % next
      end do
      atom(i) = ptr%id(1)
      ptr => ptr % next
    end do
    ! Search for non-deleted bonds and reunite the atoms as molecules:
    do i = 1, natoms-1
      do j = i+1, natoms
        if (me % bonds % find([atom(i),atom(j)])) then
          if (mol(i) /= mol(j)) then
            imin = min(mol(i),mol(j))
            imax = max(mol(i),mol(j))
            forall (k=1:natoms,mol(k) == imax) mol(k) = imin
            forall (k=1:natoms,mol(k) > imax) mol(k) = mol(k) - 1
          end if
        end if
      end do
    end do
    ! If the removed bond used to define two separate regions, there will now be two molecules:
    if (any(mol == 2)) then
      me%N = me%N + 1
      newmol = int2str(me%N)
      call writeln( "Splitting molecule", cmol, "into molecules", cmol, "and", newmol )
      ptr => me % list % first
      do while (associated(ptr))
        if (any((atom == ptr%id(1)).and.(mol == 2))) ptr%params = newmol
        ptr => ptr % next
      end do
      ! Sort the arrays in accordance with the order in mol:
      call writeln( "Reorganizing coordinates...", advance = .false. )
      allocate( sort(natoms) )
      sort = sorted( mol )
      atom = atom(sort)
      allocate( id(natoms), params(natoms) )
      ptr => me % xyz % first
      do while (associated(ptr))
        if (any(atom == ptr%id(1))) then
          pcoords => ptr
          ! Get the current order of atoms:
          do i = 1, natoms
            id(i) = ptr % id(1)
            params(i) = ptr % params
            ptr => ptr % next
          end do
          ! Reset the atoms in the new order:
          do i = 1, natoms
            pcoords % id(1) = atom(i)
            j = 1
            do while (id(j) /= atom(i))
              j = j + 1
            end do
            pcoords % params = params(j)
            pcoords => pcoords % next
          end do
        else
          ptr => ptr % next
        end if
      end do
      call writeln( " done." )
    else
      call writeln( "Opening cycle(s) in molecule", cmol )
    end if
  end subroutine tMolecule_split

  !=================================================================================================

  function tMolecule_number_of_atoms( me ) result( natoms )
    class(tMolecule), intent(in) :: me
    integer                      :: natoms(me%N)
    type(Struc), pointer :: atom
    integer :: i
    natoms = 0
    atom => me % list % first
    do while (associated(atom))
      i = str2int(atom%params)
      natoms(i) = natoms(i) + 1
      atom => atom % next
    end do
  end function tMolecule_number_of_atoms

  !=================================================================================================

  function tMolecule_count( me ) result( nmols )
    class(tMolecule), intent(in) :: me
    integer                      :: nmols(me%N)
    integer  :: imol, iatom, natoms(me%N)
    type(Struc), pointer :: ptr
    natoms = me % number_of_atoms()
    nmols = 0
    ptr => me % xyz % first
    do while (associated(ptr))
      imol = str2int(me % list % parameters( ptr%id ) )
      nmols(imol) = nmols(imol) + 1
      do iatom = 1, natoms(imol)
        ptr => ptr % next
      end do
    end do
  end function tMolecule_count

  !=================================================================================================

  function tMolecule_per_molecule( me, val_list ) result( total )
    class(tMolecule), intent(in) :: me
    type(StrucList),  intent(in) :: val_list
    real(rb)                     :: total(me%N)
    integer :: imol
    type(Struc), pointer :: ptr, pval
    total = 0.0_rb
    pval => val_list % first
    do while (associated(pval))
      ptr => me % list % first
      do while (associated(ptr))
        if (ptr % match_id(pval % id)) then
          imol = str2int( ptr % params )
          total(imol) = total(imol) + str2real( pval % params )
        end if
        ptr => ptr % next
      end do
      pval => pval % next
    end do
  end function tMolecule_per_molecule

  !=================================================================================================

  subroutine tMolecule_set_geometry( me, data, ndata, silent )
    class(tMolecule), intent(inout)        :: me
    character(sl),    intent(in)           :: data(:,:)
    integer,          intent(in)           :: ndata(size(data,1))
    logical,          intent(in), optional :: silent
    integer       :: N, i, j, k, narg, imol, imolprev, iatom, ind(3), atomsInMolecule
    character(sl) :: arg(size(data,2)), catom
    integer :: natoms(me%N)
    logical :: new_molecule
    type(Struc), pointer :: atom, ptr
    character(sl), allocatable :: prev(:), name(:), molAtoms(:)
    real(rb), allocatable :: R(:,:)
    real(rb) :: L, theta, phi, R1(3), R2(3), R3(3), x(3), y(3), z(3)
    logical :: print
    print = .not.present(silent)
    if (.not.print) print = .not.silent
    natoms = me % number_of_atoms()
    allocate( prev(maxval(natoms)) )
    N = size(ndata)
    if (print) call writeln( "Number of provided geometric data: ", int2str(N) )
    allocate( name(N), R(3,N) )
    new_molecule = .true.
    atomsInMolecule = 0
    imolprev = 0
    allocate( molAtoms(0) )
    do i = 1, N
      atomsInMolecule = atomsInMolecule + 1
      narg = ndata(i)
      arg = data(i,:)
      if ((narg < 1).or.(narg == 2).or.(narg == 6).or.(narg > 7)) &
        call error( "invalid geometric info format" )
      catom = arg(1)
      call me % list % search( [catom], atom )
      if (.not.associated(atom)) call error( "invalid atom", catom )
      if (new_molecule) then
        imol = str2int(atom%params)
        iatom = 1
        if (print) call writeln( "Processing", int2str(natoms(imol)), &
                      "geometric data for molecule", trim(int2str(imol))//":" )
        if (imol /= imolprev) then
          deallocate( molAtoms )
          allocate( molAtoms(natoms(imol)) )
          ptr => me % list % first
          j = 0
          do while (associated(ptr))
            if (str2int(ptr%params) == imol) then
              j = j + 1
              molAtoms(j) = ptr % id(1)
            end if
            ptr => ptr % next
          end do
        end if
      else if (trim(atom%params) /= trim(int2str(imol))) then
        call error( "atom", catom, "does not belong to molecule", int2str(imol) )
      else if (any(str_find([catom],prev(1:iatom)) > 0)) then
        call error( "repeated geometric info for atom", catom )
      else
        iatom = iatom + 1
      end if
      if (print) call writeln( "Data provided for atom", catom, ":", join(arg(2:narg)) )
      name(i) = catom
      select case (narg)
        case (1) ! Atom at origin
          R(:,i) = 0.0_rb
        case (3) ! Bond
          call check_atoms( arg(2:2), ind(1:1) )
          L = str2real(arg(3))
          x = real([1,0,0],rb)
          R(:,i) = R(:,ind(1)) + L*x
        case (4) ! Coordinates
          do j = 1, 3
            R(j,i) = str2real(arg(j+1))
          end do
        case (5) ! Bond and angle
          call check_atoms( arg([2,4]), ind(1:2) )
          L = str2real(arg(3))
          theta = str2real(arg(5))
          R1 = R(:,ind(1))
          R2 = R(:,ind(2))
          x = (R1 - R2)/norm(R1 - R2)
          y = real([0,1,0],rb)
          if (abs(x(2)-1.0_rb) < 0.01_rb) y = real([1,0,0],rb)
          y = y - scalar(y,x)*x
          y = y / norm(y)
          R(:,i) = R1 + L*(cosine(180-theta)*x + sine(180-theta)*y)
        case (7) ! Bond, angle, and dihedral
          call check_atoms( arg([2,4,6]), ind(1:3) )
          L = str2real(arg(3))
          theta = str2real(arg(5))
          phi = str2real(arg(7))
          R1 = R(:,ind(1))
          R2 = R(:,ind(2))
          R3 = R(:,ind(3))
          x = (R1 - R2)/norm(R1 - R2)
          y = R3 - R2 - scalar(R3 - R2,x)*x
          y = y / norm(y)
          z = cross(x,y)
          R(:,i) = R1 + L*(cosine(180-theta)*x + sine(180-theta)*(cosine(phi)*y + sine(phi)*z))
        case default
          call error( "invalid geometric info" )
      end select
      prev(iatom) = catom
      new_molecule = iatom == natoms(imol)
      if (new_molecule) then
        do j = 1, atomsInMolecule
          arg(1) = molAtoms(j)
          k = i
          do while (name(k) /= arg(1))
            k = k - 1
          end do
          arg(2:4) = real2str(R(:,k))
          call me % xyz % add( 4, arg(1:4), me % list, repeatable = .true., silent = .not.print )
        end do
        atomsInMolecule = 0
      end if
    end do
    if (.not.new_molecule) then
      call error( "geometric info for molecule", int2str(imol), "is incomplete" )
    end if
    contains
      subroutine check_atoms( atom, ind )
        character(sl), intent(in) :: atom(:)
        integer,       intent(inout) :: ind(size(atom))
        integer :: j, k
        do j = 1, size(atom)
          if (.not. me % list % find( [atom(j)])) call error( "invalid atom", atom(j) )
          if (any(name(1:i-1) == atom(j))) then
            k = 1
            do while (name(k) /= atom(j))
              k = k + 1
            end do
            ind(j) = k
          else
            call error( "no coordinates have been defined for atom", atom(j) )
          end if
        end do
      end subroutine
  end subroutine tMolecule_set_geometry

  !=================================================================================================

  subroutine tMolecule_coordinates( me, imol, N, Coord, atom, option )
    class(tMolecule), intent(inout) :: me
    integer,          intent(in)    :: imol, N, option
    real(rb),         intent(inout) :: Coord(3,N)
    character(sl),    intent(inout) :: atom(N)
    type(Struc), pointer :: current
    logical :: found
    integer :: i, narg
    character(sl) :: arg(3)
    current => me % xyz % first
    found = .false.
    do while (associated(current).and.(.not.found))
      found = str2int(me % list % parameters( current%id )) == imol
      if (.not.found) current => current % next
    end do
    if (.not.found) call error( "no coordinates for molecule", int2str(imol) )
    do i = 1, N
      if (option == 1) then ! Retrieve coordinates:
        atom(i) = current % id(1)
        call split( current % params, narg, arg )
        Coord(:,i) = [str2real(arg(1)), str2real(arg(2)), str2real(arg(3))]
      else ! Set coordinates:
        current % id(1) = atom(i)
        current % params = join(real2str(Coord(:,i)))
      end if
      current => current % next
    end do
  end subroutine tMolecule_coordinates

  !=================================================================================================

  subroutine tMolecule_align(  N, Mass, Coord, axis )
    integer,  intent(in)    :: N, axis(3)
    real(rb), intent(in)    :: Mass(N)
    real(rb), intent(inout) :: Coord(3,N)
    integer :: i
    real(rb) :: MolMass, Rcm(3), inertia(3,3), delta(3,N), MoI(3), A(3,3)
    if (N == 1) then
      Coord = 0.0_rb
    else
      ! Compute molecular mass and the center-of-mass position:
      MolMass = sum(Mass)
      forall (i=1:3) Rcm(i) = sum(Mass*Coord(i,:))/MolMass
      ! Compute inertia tensor:
      inertia = 0.0_rb
      do i = 1, N
        delta(:,i) = Coord(:,i) - Rcm
        call add_inertia( inertia, Mass(i), delta(:,i) )
      end do
      ! Diagonalize the inertia tensor and compute rotation matrix:
      call diagonalization( inertia, A, MoI )
      ! Recalculate positions in the body-fixed frame:
      Coord(axis,:) = matmul( transpose(A), delta )
      ! Invert one axis if necessary:
      if (all(axis == [1,3,2]).or.all(axis == [2,1,3]).or.all(axis == [3,2,1])) then
        Coord(1,:) = -Coord(1,:)
      end if
    end if
    contains
    !-----------------------------------------------------------------------------------------------
    subroutine add_inertia( inertia, mass, delta )
      real(rb), intent(inout) :: inertia(3,3)
      real(rb), intent(in)    :: mass, delta(3)
      inertia(1,1) = inertia(1,1) + mass*(delta(2)**2 + delta(3)**2)
      inertia(2,2) = inertia(2,2) + mass*(delta(1)**2 + delta(3)**2)
      inertia(3,3) = inertia(3,3) + mass*(delta(1)**2 + delta(2)**2)
      inertia(1,2) = inertia(1,2) - mass*delta(1)*delta(2)
      inertia(1,3) = inertia(1,3) - mass*delta(1)*delta(3)
      inertia(2,3) = inertia(2,3) - mass*delta(2)*delta(3)
      inertia(2,1) = inertia(1,2)
      inertia(3,1) = inertia(1,3)
      inertia(3,2) = inertia(2,3)
    end subroutine add_inertia
!---------------------------------------------------------------------------------------------------
  end subroutine  tMolecule_align

  !=================================================================================================
  ! Numerical diagonalization of 3x3 matrcies
  ! Copyright (C) 2006  Joachim Kopp
  !=================================================================================================

  pure subroutine diagonalization( matrix, q, w )
    real(rb), intent(in)  :: matrix(3,3)
    real(rb), intent(out) :: q(3,3), w(3)

    ! ----------------------------------------------------------------------------
    ! Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
    ! matrix A using Cardano's method for the eigenvalues and an analytical
    ! method based on vector cross products for the eigenvectors.
    ! Only the diagonal and upper triangular parts of A need to contain
    ! meaningful values. However, all of A may be used as temporary storage
    ! and may hence be destroyed.
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   matrix: The symmetric input matrix
    !   q: Storage buffer for eigenvectors
    !   w: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------

    real(rb), parameter :: zero = 0.0_rb, one = 1.0_rb, two = 2.0_rb, three = 3.0_rb
    real(rb), parameter :: half = one/two, third = one/three
    real(rb), parameter :: eps = epsilon(one)
    real(rb), parameter :: sqrt3 = sqrt(3.0_rb)

    integer  :: i, j
    real(rb) :: a(3,3), norm, n1, n2, w1, w2, w3, thresh, t, wmax8eps
    real(rb) :: m, c1, c0, de, dd, ee, ff, p, sqrtp, r, c, s, phi
    logical  :: success

    a = matrix

    ! Calculate the eigenvalues of a symmetric 3x3 matrix a using Cardano's
    ! analytical algorithm. Only the diagonal and upper triangular parts of A are
    ! accessed. The access is read-only.
    de = a(1,2)*a(2,3)
    dd = a(1,2)**2
    ee = a(2,3)**2
    ff = a(1,3)**2
    m  = a(1,1) + a(2,2) + a(3,3)
    c1 = (a(1,1)*a(2,2) + a(1,1)*a(3,3) + a(2,2)*a(3,3)) - (dd + ee + ff)
    c0 = 27.0_rb*(a(3,3)*dd + a(1,1)*ee + a(2,2)*ff - a(1,1)*a(2,2)*a(3,3) - two*a(1,3)*de)

    p = m*m - 3.0_rb*c1
    r = m*(p - 1.5_rb*c1) - half*c0
    sqrtp = sqrt(abs(p))
    phi = third*atan2(sqrt(abs(6.75_rb*c1*c1*(p - c1) + c0*(r + 0.25_rb*c0))),r)

    c = sqrtp*cos(phi)
    s = (one/sqrt3)*sqrtp*sin(phi)

    p = third*(m - c)
    w1 = p + c
    w2 = p - s
    w3 = p + s

    ! Sort eigenvalues:
    if (abs(w1) > abs(w3)) call swap( w1, w3 )
    if (abs(w1) > abs(w2)) call swap( w1, w2 )
    if (abs(w2) > abs(w3)) call swap( w2, w3 )
    w = [w1, w2, w3]

    wmax8eps = 8.0_rb*eps*abs(w1)
    thresh = wmax8eps**2

    ! Prepare calculation of eigenvectors
    n1 = a(1,2)**2 + a(1,3)**2
    n2 = a(1,2)**2 + a(2,3)**2
    q(1,1) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    q(1,2) = q(1,1)
    q(2,1) = a(1,3)*a(1,2) - a(2,3)*a(1,1)
    q(2,2) = q(2,1)
    q(3,2) = a(1,2)**2

    ! Calculate first eigenvector by the formula v(1) = (A - lambda(1)).e1 x (A - lambda(1)).e2
    a(1,1) = a(1,1) - w1
    a(2,2) = a(2,2) - w1
    q(:,1) = [q(1,2) + a(1,3)*w1, q(2,2) + a(2,3)*w1, a(1,1)*a(2,2) - q(3,2)]
    call compute_eigenvector( q(:,1), a, n1, n2 )

    ! Prepare calculation of second eigenvector     
    t = w1 - w2

    ! Is this eigenvalue degenerate?
    if (abs(t) > wmax8eps) then

      ! For non-degenerate eigenvalue, calculate second eigenvector by the formula
      !         v[1] = (A - lambda[1]).e1 x (A - lambda[1]).e2
      a(1,1) = a(1,1) + t
      a(2,2) = a(2,2) + t
      q(:,2) = [q(1,2) + a(1,3)*w2, q(2,2) + a(2,3)*w2, a(1,1)*a(2,2) - q(3,2)]
      call compute_eigenvector( q(:,2), a, n1, n2 )

    else

      ! For degenerate eigenvalue, calculate second eigenvector according to
      !         v[1] = v(1) x (A - lambda[1]).e[i]

      ! This would really get too complicated if we could not assume all of A to
      !       contain meaningful values.
      a(2,1) = a(1,2)
      a(3,1) = a(1,3)
      a(3,2) = a(2,3)
      a(1,1) = a(1,1) + w1
      a(2,2) = a(2,2) + w1
      i = 0
      success = .false.
      do while ((i < 3).and.(.not.success))
        i = i + 1
        a(i,i) = a(i,i) - w2
        n1 = sum(a(:,i)**2)
        success = n1 > thresh
        if (success) then
          q(:,2) = cross_product( q(:,1), a(:,i) )
          norm = sum(q(:,2)**2)
          success = norm > (256.0_rb*eps)**2*n1
          if (success) q(:,2) = q(:,2)*sqrt(one/norm)
        end if
      end do

      ! This means that any vector orthogonal to v(1) is an EV.
      if (.not.success) then
        i = 1
        do while (q(i,1) == zero)
          i = i + 1
        end do
        j = 1 + mod(i,3)
        norm = one/sqrt(q(i,1)**2 + q(j,1)**2)
        q(i,2) =  q(j,1)*norm
        q(j,2) = -q(i,1)*norm
        q(1+mod(i+1,3),2) = zero
      end if
    end if

    ! Calculate third eigenvector according to v[2] = v(1) x v[1]
    q(:,3) = cross_product( q(:,1), q(:,2) )

    contains
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine compute_eigenvector( q, a, n1tmp, n2tmp )
        real(rb), intent(inout) :: q(3)
        real(rb), intent(in)    :: a(3,3), n1tmp, n2tmp

        real(rb) :: norm, n1, n2, error, t, f

        norm = sum(q**2)
        n1 = n1tmp + a(1,1)**2
        n2 = n2tmp + a(2,2)**2
        error = n1*n2

        ! If the first column is zero, then (1, 0, 0) is an eigenvector
        if (n1 <= thresh) then
          q = [one, zero, zero]

        ! If the second column is zero, then (0, 1, 0) is an eigenvector
        else if (n2 <= thresh) then
          q = [zero, one, zero]

        ! If angle between A(*,1) and A(*,2) is too small, don't use
        !  cross product, but calculate v ~ (1, -A0/A1, 0)
        else if (norm < (64.0_rb*eps)**2*error) then
          t = abs(a(1,2))
          f = -a(1,1)/a(1,2)
          if (abs(a(2,2)) > t) then
            t = abs(a(2,2))
            f = -a(1,2)/a(2,2)
          end if
          if (abs(a(2,3)) > t) f = -a(1,3)/a(2,3)
          norm = one/sqrt(one + f**2)
          q = [norm, f*norm, zero]

        ! This is the standard branch
        else
          q = q*sqrt(one/norm)
        end if

      end subroutine compute_eigenvector
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pure subroutine swap( a, b )
        real(rb), intent(inout) :: a, b
        real(rb) :: c
        c = a; a = b; b = c
      end subroutine swap
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine diagonalization

  !=================================================================================================

  pure function cross_product(a, b) result( c )
    real(rb), intent(in) :: a(3), b(3)
    real(rb)             :: c(3)
    c = [ a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1) ]
  end function cross_product

  !=================================================================================================

end module mMolecule
