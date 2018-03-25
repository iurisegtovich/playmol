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

! TODO: Implement improper dihedrals (check order conformity with LAMMPS)
! TODO: Implement other openmm torsion models

  subroutine tPlaymol_write_openmm( me, unit, keywords )
    class(tPlaymol), intent(inout)        :: me
    integer,         intent(in)           :: unit
    character(*),    intent(in)           :: keywords

    integer :: ntotal
    real(rb) :: length, energy, angle
    logical :: guess
    character(sl) :: lj14, coul14

    integer, allocatable :: natoms(:)
    character(sl), allocatable :: atom(:), atom_type(:), raw_atom(:), charge(:), element(:), mass(:)

    call process( keywords )

    natoms = me % molecules % number_of_atoms()
    ntotal = sum(natoms)
    allocate( atom(ntotal), &
              atom_type(ntotal), &
              raw_atom(ntotal), &
              charge(ntotal), &
              element(ntotal), &
              mass(ntotal) )

    block
      integer :: i, imol, n, mol(ntotal)
      type(Struc), pointer :: current
      n = 0
      current => me % molecules % list % first
      do while (associated(current))
        imol = str2int(current % params)
        if (natoms(imol) > 0) then
          n = n + 1
          mol(n) = imol
          atom(n) = current % id(1)
        end if
        current => current % next
      end do
      atom = atom(sorted(mol))
      do i = 1, ntotal
        atom_type(i) = me % atom_list % parameters( atom(i:i) )
        raw_atom(i) = me % raw_atom_list % parameters( atom(i:i) )
        charge(i) = me % charge_list % parameters( atom(i:i), default = "0" )
        call me % element_and_mass( atom_type(i), element(i), mass(i) )
        if (guess.and.(element(i) == "UA")) element(i) = element_guess( mass(i) )
      end do
    end block

    write(unit,'("<ForceField>")')
    call atom_types()

    call residues()
    call bond_types()

    call angle_types()

    call dihedral_types()

    call non_bonded_model()

    write(unit,'("</ForceField>")')

    contains

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine process( keywords )
        character(*), intent(in) :: keywords
        integer :: i, narg
        real(rb) :: rval
        character(sl) :: keyword, value, arg(40)
        call split( keywords, narg, arg)
        if (mod(narg, 2) == 1) call error( "invalid write openmm command" )
        do i = 1, narg/2
          keyword = arg(2*i-1)
          value = arg(2*i)
          select case (keyword)
            case ("length", "energy", "angle", "lj14", "coul14")
              rval = str2real(value)
              if (rval < 0.0_rb) call error( "invalid", keyword, "parameter value" )
              select case (keyword)
                case ("length")
                  length = rval
                case ("energy")
                  energy = rval
                case ("angle")
                  angle = rval
                case ("lj14")
                  lj14 = value
                case ("coul14")
                  coul14 = value
              end select
            case ("elements")
              if (.not.any(value == ["yes", "no "])) call error( "invalid write openmm command" )
              guess = (value == "yes")
            case default
              call error( "invalid write openmm command" )
          end select
        end do
      end subroutine process

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine atom_types()
        integer :: i, n
        type(StrucList) :: local_list = StrucList( "atom type", 1 )
        character(sl), parameter :: p3(3) = [character(sl) :: "name", "class", "mass"], &
                                    p4(4) = [character(sl) :: "name", "class", "element", "mass"]

        write(unit,'(2X,"<AtomTypes>")')
        n = 0
        do i = 1, ntotal 
          if (.not. local_list % find(atom_type(i:i))) then
            call local_list % add(1, atom_type(i:i), silent = .true.)
            if ((element(i) == "EP").or.(element(i) == "UA")) then
              call items(4, "Type", p3, [atom_type(i), atom_type(i), mass(i)])
            else
              call items(4, "Type", p4, [atom_type(i), atom_type(i), element(i), mass(i)])
            end if
          end if
        end do
        call local_list % destroy(silent = .true.)
        write(unit,'(2X,"</AtomTypes>")')
      end subroutine atom_types

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine residues()
        integer :: imol
        character(sl) :: name(me%molecules%N)

        write(unit,'(2X,"<Residues>")')
        call me % get_molecule_names( name )
        do imol = 1, me%molecules%N
          call residue( imol, name(imol) )
        end do
        write(unit,'(2X,"</Residues>")')

      end subroutine residues

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine residue( imol, resname )
        integer,       intent(in) :: imol
        character(sl), intent(in) :: resname

        integer :: i, j, first, last
        integer, allocatable :: indx(:), pos(:)
        type(Struc), pointer :: current

        character(sl), parameter :: pa(3) = [character(sl) :: "name", "type", "charge"], &
                                    pb(2) = [character(sl) :: "atomName1", "atomName2"]

        last = sum(natoms(1:imol))
        first = last - natoms(imol) + 1

        write(unit,'(4X,"<Residue ",A,">")') trim(item("name", resname))

        do i = first, last
          call items(6, "Atom", pa, [raw_atom(i), atom_type(i), charge(i)])
        end do

        ! Virtual sites:
        do i = first, last
          if (element(i) == "EP") call virtual_site( i )
        end do

        ! Bonds:
        indx = [(i,i=first,last)]
        current => me % bond_list % first
        do while (associated(current))
          pos = [(pack(indx, atom(indx) == current%id(j)), j=1, 2)]
          if (size(pos) == 2) call items(6, "Bond", pb, raw_atom(pos))
          current => current % next
        end do

        write(unit,'(4X,"</Residue>")')
      end subroutine residue

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine virtual_site( i )
        integer, intent(in) :: i

        real(rb), parameter :: tol = 1.0E-4_rb

        integer :: k, n
        character(sl) :: string
        integer :: partner(3)
        real(rb) :: a(3,3), b(3), w(3), axis(3,3)
        character(sl) :: xyz(3), average(3)
        integer, allocatable :: pos(:)
        character(sl), allocatable :: properties(:), values(:)
        type(Struc), pointer :: current

        string = me % molecules % xyz % parameters( [atom(i)] )
        call split(string, k, xyz)
        b = [(str2real(xyz(k)),k=1,3)]
        n = 0
        current => me % link_list % first
        do while (associated(current).and.(n < 4))
          pos = pack([2,1], current%id == atom(i))
          if (size(pos) > 0) then
            n = n + 1
            string = me % molecules % xyz % parameters( [current%id(pos(1))] )
            call split(string, k, xyz)
            a(:,n) = [(str2real(xyz(k)),k=1,3)]
            partner(n:n) = pack([(k,k=1,size(atom))], atom == current%id(pos(1)))
          end if
          current => current % next
        end do
        if (n == 2) then
          ! Test for colinearity:
          if (colinear(a(:,1), a(:,2), b)) then
            w(1:2) = gaussian_elimination( a(1:2,1:2), b(1:2) )
            average(1:2) = [character(sl) :: "1", "2"]
          else
            call error( "VirtualSite type average2 requires colinearity")
          end if
        else if (n == 3) then
          axis(:,1) = unit_vector(a(:,1), a(:,2))
          axis(:,2) = unit_vector(a(:,1), a(:,3))
          axis(:,3) = cross_product(axis(:,1), axis(:,2))
          w = gaussian_elimination( axis, b - a(:,1) )
          if (abs(w(3)) < tol) then
            w = gaussian_elimination( a, b )
            average = [character(sl) :: "1", "2", "3"]
          else
            average = [character(sl) :: "12", "13", "Cross"]
          end if
        else
          call error( "Extra particle must be linked to 2 or 3 atoms" )
        end if
        properties = [character(sl) :: "type", "siteName", &
                      ("atomName"//int2str(k), k=1, n), ("weight"//average(k), k=1, n)]
        values = [character(sl) :: "average"//average(n), raw_atom(i), &
                  raw_atom(partner), (float2str(w(k)),k=1,3)]
        call items(6, "VirtualSite", properties, values)
      end subroutine virtual_site

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine bond_types()

        integer :: narg
        real(rb) :: K, r0
        character(sl) :: arg(20)
        type(StrucList) :: list
        type(Struc), pointer :: current
        character(sl), parameter :: p(4) = [character(sl) :: "type1", "type2", "length", "k"]

        write(unit,'(2X,"<HarmonicBondForce>")')
        list = local_list( me % bond_list, me % bond_type_list, .false. )
        current => list % first
        do while (associated(current))
          call split( current%params, narg, arg )
          if ((narg == 2).and.all(is_real(arg(1:2)))) then
            K = 2.0_rb*str2real(arg(1)) * energy/length**2
            r0 = str2real(arg(2)) * length
          else if (arg(1) == "harmonic") then
            K = 2.0_rb*str2real(arg(2)) * energy/length**2
            r0 = str2real(arg(3)) * length
          else if (arg(1) == "zero") then
            current => current % next
            cycle
          else
            call error( "harmonic bond model required" )
          end if
          call items(4, "Bond", p, [current%id, real2str([r0, K])])
        current => current % next
        end do
        call list % destroy(silent = .true.)
        write(unit,'(2X,"</HarmonicBondForce>")')
      end subroutine bond_types

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine angle_types()

        integer :: narg
        real(rb) :: K, theta0
        character(sl) :: arg(20)
        type(StrucList) :: list
        type(Struc), pointer :: current
        character(sl), parameter :: p(5) = ["type1", "type2", "type3", "angle", "k    "]

        write(unit,'(2X,"<HarmonicAngleForce>")')
        list = local_list( me % angle_list, me % angle_type_list, .false. )
        current => list % first
        do while (associated(current))
          call split( current%params, narg, arg )
          if ((narg == 2).and.all(is_real(arg(1:2)))) then
            K = 2.0_rb*str2real(arg(1)) * energy
            theta0 = str2real(arg(2)) * angle
          else if (arg(1) == "harmonic") then
            K = 2.0_rb*str2real(arg(2)) * energy
            theta0 = str2real(arg(3)) * angle
          else
            call error( "harmonic angle model required" )
          end if
          call items(4, "Angle", p, [current%id, real2str([theta0, K])])
          current => current % next
        end do
        call list % destroy(silent = .true.)
        write(unit,'(2X,"</HarmonicAngleForce>")')
      end subroutine angle_types

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine dihedral_types()

        integer :: narg, n, i, j
        real(rb) :: K, phase
        logical :: same
        character(sl) :: arg(20), type_id(4), props(9), values(9)
        type(StrucList) :: list
        type(Struc), pointer :: current
        character(sl), parameter :: empty(4) = "", &
                                    pt(4) = ["type1", "type2", "type3", "type4"], &
                                    pp(3) = [character(sl) :: "periodicity", "phase", "k"]
        
        write(unit,'(2X,"<PeriodicTorsionForce>")')
        list = local_list( me % dihedral_list, me % dihedral_type_list, .true. )
        current => list % first
        do while (associated(current))
          type_id = current%id
          i = 0
          same = .true.
          do while (associated(current) .and. same)
            i = i + 1
            call split( current%params, narg, arg )
            if ((arg(1) == "harmonic").or.((narg == 3).and.all(is_real(arg(1:narg))))) then
              K = str2real(arg(2)) * energy
              n = str2int(arg(4))
              phase = 90*(1 - str2int(arg(3))) * angle
            else if ((arg(1) == "charmm").or.((narg == 4).and.all(is_real(arg(1:narg))))) then
              K = str2real(arg(2)) * energy
              n = str2int(arg(3))
              phase = str2int(arg(4)) * angle
            else
              call error( "harmonic or charmm dihedral model required" )
            end if
            forall (j=1:3) props(3*(i-1)+j) = trim(pp(j))//int2str(i)
            values(3*i-2:3*i) = [int2str(n), real2str(phase), real2str(K)]
            current => current % next
            if (associated(current)) same = all(current%id == type_id)
          end do
          call items(4, "Proper", [pt, props(1:3*i)], &
                                  [merge(empty, type_id, type_id == "*"), values(1:3*i)])
        end do
        call list % destroy(silent = .true.)
        write(unit,'(2X,"</PeriodicTorsionForce>")')
      end subroutine dihedral_types

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine non_bonded_model()

        integer :: narg
        character(sl) :: arg(20), eps, sig
        type(StrucList) :: list
        type(Struc), pointer :: current
        real(rb), parameter :: tol = 1.0E-8_rb
        character(sl), parameter :: zero = "0", &
                                    p(4) = [character(sl) :: "type", "charge", "sigma", "epsilon"]

        write(unit,'(2X,"<NonbondedForce ",A,X,A,">")') trim(item("coulomb14scale", coul14)), &
                                                        trim(item("lj14scale", lj14))
        list = local_list( me % atom_list, me % atom_type_list, .false. )
        current => list % first
        do while (associated(current))
          call split( current%params, narg, arg )
          if ((narg == 2).and.all(is_real(arg(1:2)))) then
            eps = real2str(str2real(arg(1)) * energy)
            sig = real2str(str2real(arg(2)) * length)
          else if (arg(1)(1:2) == "lj") then
            eps = real2str(str2real(arg(2)) * energy)
            sig = real2str(str2real(arg(3)) * length)
          else if ((arg(1) == "zero").or.(arg(1)(1:4) == "coul")) then
            eps = zero
            sig = zero
          else
            call error( "Lennard-Jones potential model required" )
          end if
          call items(4, "Atom", p, [current%id(1), zero, sig, eps])
          current => current % next
        end do
        call list % destroy(silent = .true.)
        write(unit,'(2X,"</NonbondedForce>")')
      end subroutine non_bonded_model

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      subroutine items( ident, title, property, value )
        integer,      intent(in) :: ident
        character(*), intent(in) :: title, property(:), value(:)
        integer :: i
        write(unit,'("'//repeat(" ",ident)//'","<",A,X,A,"/>")') trim(title), &
          trim(join([(item(property(i), value(i)), i=1, size(property))]))
      end subroutine items

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      elemental character(sl) function item( property, value )
        character(*), intent(in) :: property, value
        item = trim(property)//"="""//trim(value)//""""
      end function item

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      function element_guess( mass ) result( element )
        character(sl), intent(in) :: mass
        character(sl)             :: element
        element = me%elements(minloc(abs(me%masses - str2real(mass)), dim = 1))
      end function element_guess

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      function find( a, b )
        character(sl), intent(in) :: a(:), b(:)
        logical                   :: find(size(b))
        integer :: i
        find = .false.
        do i = 1, size(a)
          find = find .or. (b == a(i))
        end do
      end function find

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      function local_list( struc_list, type_list, macros ) result( list )
        type(StrucList), intent(in) :: struc_list, type_list
        logical,         intent(in) :: macros
        type(StrucList)             :: list

        integer :: i, indx(ntotal)
        integer, allocatable :: pos(:)
        character(sl) :: type_id(struc_list%number), arg(struc_list%number+1)
        type(Struc), pointer :: current, ptr

        list%number = struc_list%number

        indx = [(i,i=1,ntotal)]
        current => struc_list % first
        do while (associated(current))
          pos = pack(indx, find(current%id, atom))
          if (size(pos) == list%number) then
            call me % get_types( current%id, type_id )
            if (.not. list % find( type_id )) then
              ptr => type_list % first
              do while (associated(ptr))
                if (ptr % match_id( type_id )) then
                  if (macros) then
                    arg = [merge(ptr%id, type_id, ptr%id == "*"), ptr%params]
                  else
                    arg = [type_id, ptr%params]
                  end if
                  call list % add(list%number+1, arg, repeatable = .true., silent = .true.)
                end if
                ptr => ptr % next
              end do
            end if
          end if
          current => current % next
        end do
      end function local_list

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  end subroutine tPlaymol_write_openmm
