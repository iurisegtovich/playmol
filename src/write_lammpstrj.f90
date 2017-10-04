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

  subroutine tPlaymol_write_lammpstrj( me, unit ) !box snapshot in lammpstrj dump format
    class(tPlaymol), intent(inout) :: me
    integer,         intent(in)    :: unit
    integer :: molcount(me%molecules%N), iatom, itype, i, imol, jmol, narg, natoms(me%molecules%N)
    real(rb) :: mass(me%molecules%N)
    character(sl) :: limits, arg(1)
    type(Struc), pointer :: current, atom_type
    logical :: found
    if (.not.me % box % exists()) call error( "simulation box has not been defined" )
    write(unit,'("ITEM: TIMESTEP",/,"0")')
    write(unit,'("ITEM: NUMBER OF ATOMS")')
    write(unit,'(A)') trim(int2str(me % molecules % xyz % count))
    write(unit,'("ITEM: BOX BOUNDS pp pp pp")')
    molcount = me % molecules % count()
    mass = me % molecules % per_molecule( me % atom_masses )
    call me % box % compute( sum(molcount*mass) )
    do i = 1, 3
      limits = join(real2str( me%box%length(i)*[-0.5_rb,+0.5_rb] ))
      write(unit,'(A)') trim(limits)
    end do
!    write(unit,'("ITEM: ATOMS id mol type x y z ix iy iz ")')
    write(unit,'("ITEM: ATOMS id mol type x y z")')
    natoms = me % molecules % number_of_atoms()
    current => me % molecules % xyz % first !first atom
    iatom = 0
    jmol = 0
    do while (associated(current))
      jmol = jmol + 1 !increasing counter for each molecule in the box (according to bond or link info - i guess)
      imol = str2int(me % molecules % list % parameters( current % id ) ) !which molecule is being point to now
      do i = 1, natoms(imol) !each atom in this molecule
        iatom = iatom + 1 !increasing counter for each atom in the box
        call split( me % atom_list % parameters( current % id ), narg, arg )
        !arg is now the string identification of atom type being looped in the block below; seeking match between atom types and current iatom
        atom_type => me % atom_type_list % first
        itype = 0
        found = .false.
        do while (associated(atom_type).and.(.not.found))
          if (atom_type % usable) then
            itype = itype + 1
            found = atom_type % match_id( arg )
!            print*, "FOUND: atom with in-box-number", iatom, "is of itype/type", itype, "/", arg
          end if
          if (.not.found) atom_type => atom_type % next
        end do
!        write(unit,'(3(A,X),"0 0 0")') trim(join(int2str([iatom,jmol,itype]))), trim(current%params) !integer identification
        write(unit,'(3(A,X))') trim(join([int2str([iatom,jmol]),atom_type%id(1)])), trim(current%params) !string identification
        current => current % next
      end do
    end do
  end subroutine tPlaymol_write_lammpstrj

