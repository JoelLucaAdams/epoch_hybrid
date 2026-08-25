! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
!
! hy_resistivity.F90
!
! This module holds all the subroutines related to resistivity, and also the
! subroutines for updating the ionisation state and Coulomb logarithm (as these
! are both used in the reduced Lee-More resistivity model)

MODULE hy_resistivity
#ifdef HYBRID

  USE hy_shared

  IMPLICIT NONE

  REAL(num), PRIVATE, PARAMETER :: rlm_c_hot = &
      128.0_num*epsilon0**2/q0**4 * SQRT(0.5_num*pi*m0)
  REAL(num), PRIVATE, PARAMETER :: rlm_c_ve = 3.0_num*kb/m0
  REAL(num), PRIVATE, PARAMETER :: rlm_c_sphere = 4.0_num*pi/3.0_num

CONTAINS

  SUBROUTINE update_ionisation

    ! Loops over all cells and updates the ionisation charge state of the
    ! background solid. In the case of compound solids, the ion properties are
    ! averaged over species

    INTEGER :: ix, iy, iz

    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng
          ion_charge(ix,iy,iz) = thomas_fermi_ionisation(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE update_ionisation



  SUBROUTINE update_coulomb_logarithm

    ! Loops over all cells and updates the Coulomb logarithm of the background
    ! solid. In the case of compound solids, the ion properties are averaged
    ! over species

    INTEGER :: ix, iy, iz

    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng
          ion_cou_log(ix,iy,iz) = lee_more_coulomb_log(ix,iy,iz)
        END DO
      END DO
    END DO

  END SUBROUTINE update_coulomb_logarithm



  SUBROUTINE update_resistivity

    ! Loops over the resistivity grid and updates each point based on the
    ! resistivity model in that cell

    INTEGER :: ix, iy, iz

    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng
          SELECT CASE (resistivity_model(ix,iy,iz))
          CASE(c_resist_vacuum)
            resistivity(ix,iy,iz) = 0.0_num
          CASE(c_resist_milchberg)
            resistivity(ix,iy,iz) = calc_resistivity_milchberg(hy_te(ix,iy,iz))
          CASE(c_resist_plastic)
            resistivity(ix,iy,iz) = calc_resistivity_plastic(hy_te(ix,iy,iz))
          CASE(c_resist_rlm)
            resistivity(ix,iy,iz) = calc_resistivity_rlm(ix,iy,iz)
          CASE(c_resist_table)
            resistivity(ix,iy,iz) = calc_resistivity_table(ix,iy,iz)
          END SELECT
        END DO
      END DO
    END DO

  END SUBROUTINE update_resistivity



  FUNCTION calc_resistivity_milchberg(te)

    ! Calculates the resistivity for the background temperature te, using the
    ! data from H. M. Milchberg, et al, 1988. Phys. Rev. Lett., 61(20), p.2364.
    !
    ! This applies to Al only. The fit to this data comes from Davies, et al,
    ! (2002). Phys. Rev. E, 65(2), 026407

    REAL(num), INTENT(IN) :: te
    REAL(num) :: te_ev
    REAL(num) :: calc_resistivity_milchberg

    ! Convert te from Kelvin to eV
    te_ev = kelvin_to_ev * te

    ! Calculate resistivity
    calc_resistivity_milchberg = te_ev &
        / (5.0e6_num + 170.0_num*te_ev**2.5_num + 3.0e5_num*te_ev)

  END FUNCTION calc_resistivity_milchberg



  FUNCTION calc_resistivity_plastic(te)

    ! Calculates the resistivity for the background temperature te, using the
    ! heuristic model for plastic from Davies, et al, (1999). Phys. Rev. E,
    ! 59(5), 6032.

    REAL(num), INTENT(IN) :: te
    REAL(num) :: te_ev
    REAL(num) :: calc_resistivity_plastic

    ! Convert te from Kelvin to eV
    te_ev = kelvin_to_ev * te

    ! Calculate resistivity
    calc_resistivity_plastic = 1.0_num / (4.3e5_num + 1.3e3_num*te_ev**1.5_num)

  END FUNCTION calc_resistivity_plastic



  FUNCTION calc_resistivity_rlm(ix, iy, iz)

    ! Calculates the resistivity for the background temperature te, using a
    ! reduced form of the Lee-More resistivity. The approximations used here are
    ! similar to those used in the hybrid PIC code Zephyros.
    !
    ! The full Lee-More model involves calculating the chemical potential, and
    ! can be found in Lee, More (1984). Phys. Fluids. 27(5). Note that this
    ! paper is written in Gaussian CGS units, whereas EPOCH is in SI

    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: rlm_z_star, rlm_cou_log, rlm_te, rlm_ni
    REAL(num) :: ve
    REAL(num) :: ion_sphere_rad
    REAL(num) :: tau_cold, tau_hot, tau
    REAL(num) :: calc_resistivity_rlm

    ! Copy out array variables
    rlm_z_star = ion_charge(ix,iy,iz)
    rlm_cou_log = ion_cou_log(ix,iy,iz)
    rlm_te = hy_te(ix,iy,iz)
    rlm_ni = ion_ni(ix,iy,iz)

    ! Ion sphere radius (interatomic distance)
    ion_sphere_rad = (rlm_c_sphere*rlm_ni)**(-1.0_num/3.0_num)

    ! Thermal electron velocity (see between Lee-More eq. 20 & 21)
    ve = SQRT(rlm_c_ve * rlm_te)

    ! Collision time just above melting point (Lee-More section 3 end), with fit
    ! parameter rlm_1
    tau_cold = ion_sphere_rad / ve * rlm_1

    ! Collision time in the high temperature limit (Lee-More eq. 27 & 28a)
    ! This is tau*A^alpha, as used in Lee-More 23a
    tau_hot = rlm_c_hot * (kb*rlm_te)**1.5_num &
        / (rlm_z_star**2 * rlm_ni * rlm_cou_log)

    ! Approximation: as temperature falls, tau_hot becomes unphysically small,
    ! so switch to tau_cold
    tau = MAX(tau_cold, tau_hot)

    ! Inverted conductivity eq. 23a. Added a global scaling fit parameter
    calc_resistivity_rlm = m0/(rlm_z_star * rlm_ni * tau * q0**2) * rlm_2

  END FUNCTION calc_resistivity_rlm


  SUBROUTINE setup_resistivity_tables

    ! Reads each solid's resistivity table from disk (rank 0 only) then
    ! broadcasts the data to all MPI ranks.
    !
    ! File format (plain ASCII):
    !   Line 1: n_rho n_te n_ti
    !   Line 2: density values [kg/m^3], space-separated
    !   Line 3: electron temperature values [eV], space-separated
    !   Line 4: ion temperature values [eV], space-separated
    !   Lines 5+: resistivity values [Ohm.m], one row per (i_rho, i_te) pair,
    !             with n_ti values per row, varying i_rho slowest, i_te next,
    !             i_ti fastest.

    USE mpi

    INTEGER :: isolid, i_rho, i_te, n_rho, n_te, n_ti

    DO isolid = 1, solid_count
      IF (solid_array(isolid)%res_model /= c_resist_table) CYCLE

      IF (rank == 0) THEN
        OPEN(unit=lu, &
            file=TRIM(solid_array(isolid)%resistivity_table_location), &
            status='OLD')
        READ(lu,*) n_rho, n_te, n_ti
        ALLOCATE(solid_array(isolid)%rho_table(n_rho))
        ALLOCATE(solid_array(isolid)%te_table(n_te))
        ALLOCATE(solid_array(isolid)%ti_table(n_ti))
        ALLOCATE(solid_array(isolid)%eta_table(n_rho, n_te, n_ti))

        READ(lu,*) solid_array(isolid)%rho_table
        READ(lu,*) solid_array(isolid)%te_table
        READ(lu,*) solid_array(isolid)%ti_table
        DO i_te = 1, n_te
          DO i_rho = 1, n_rho
            READ(lu,*) solid_array(isolid)%eta_table(i_rho, i_te, :)
          END DO
        END DO
        CLOSE(unit=lu)
      END IF

      ! Array sizes are cast to all other ranks
      CALL MPI_BCAST(n_rho, 1, MPI_INTEGER, 0, comm, errcode)
      CALL MPI_BCAST(n_te,  1, MPI_INTEGER, 0, comm, errcode)
      CALL MPI_BCAST(n_ti,  1, MPI_INTEGER, 0, comm, errcode)

      ! All other ranks know the size from the MPI_BCAST and allocate arrays
      IF (rank /= 0) THEN
        ALLOCATE(solid_array(isolid)%rho_table(n_rho))
        ALLOCATE(solid_array(isolid)%te_table(n_te))
        ALLOCATE(solid_array(isolid)%ti_table(n_ti))
        ALLOCATE(solid_array(isolid)%eta_table(n_rho, n_te, n_ti))
      END IF

      ! Arrays are cast to all other ranks
      CALL MPI_BCAST(solid_array(isolid)%rho_table, n_rho, mpireal, 0, &
          comm, errcode)
      CALL MPI_BCAST(solid_array(isolid)%te_table, n_te, mpireal, 0, &
          comm, errcode)
      CALL MPI_BCAST(solid_array(isolid)%ti_table, n_ti, mpireal, 0, &
          comm, errcode)
      CALL MPI_BCAST(solid_array(isolid)%eta_table, n_rho * n_te * n_ti, &
          mpireal, 0, comm, errcode)
    END DO

  END SUBROUTINE setup_resistivity_tables



  FUNCTION calc_resistivity_table(ix, iy, iz)

    ! Calculates resistivity using trilinear interpolation on the table loaded
    ! for the dominant solid at this cell. The dominant solid index is stored in
    ! solid_index_model, set once at initialisation alongside resistivity_model,
    ! using the same highest-electron-density criterion.

    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: calc_resistivity_table

    INTEGER :: i_sol
    REAL(num) :: rho_mass, te, ti

    i_sol = solid_index_model(ix,iy,iz)
    rho_mass = solid_array(i_sol)%ion_density(ix,iy,iz) &
        * solid_array(i_sol)%mass_no * amu
    te = hy_te(ix,iy,iz)
    IF (use_ion_temp) THEN
      ti = hy_ti(ix,iy,iz)
    ELSE
      ti = te
    END IF

    calc_resistivity_table = interp_trilinear(rho_mass, te, ti, i_sol)

  END FUNCTION calc_resistivity_table



  FUNCTION interp_trilinear(rho_in, te_in, ti_in, i_sol)

    ! Trilinear (3-parameter linear Lagrange) interpolation on the resistivity
    ! table for solid i_sol.

    REAL(num), INTENT(IN) :: rho_in, te_in, ti_in
    INTEGER, INTENT(IN) :: i_sol
    REAL(num) :: interp_trilinear

    INTEGER :: i1r, i2r, i1e, i2e, i1i, i2i
    REAL(num) :: fr, fe, fi
    REAL(num) :: w000, w100, w010, w110, w001, w101, w011, w111

    CALL find_bracket_indices(rho_in, solid_array(i_sol)%rho_table, &
        i1r, i2r, fr)
    CALL find_bracket_indices(te_in,  solid_array(i_sol)%te_table,  &
        i1e, i2e, fe)
    CALL find_bracket_indices(ti_in,  solid_array(i_sol)%ti_table,  &
        i1i, i2i, fi)

    w000 = (1.0_num-fr) * (1.0_num-fe) * (1.0_num-fi)
    w100 = fr           * (1.0_num-fe) * (1.0_num-fi)
    w010 = (1.0_num-fr) * fe           * (1.0_num-fi)
    w110 = fr           * fe           * (1.0_num-fi)
    w001 = (1.0_num-fr) * (1.0_num-fe) * fi
    w101 = fr           * (1.0_num-fe) * fi
    w011 = (1.0_num-fr) * fe           * fi
    w111 = fr           * fe           * fi

    interp_trilinear = &
        w000 * solid_array(i_sol)%eta_table(i1r, i1e, i1i) + &
        w100 * solid_array(i_sol)%eta_table(i2r, i1e, i1i) + &
        w010 * solid_array(i_sol)%eta_table(i1r, i2e, i1i) + &
        w110 * solid_array(i_sol)%eta_table(i2r, i2e, i1i) + &
        w001 * solid_array(i_sol)%eta_table(i1r, i1e, i2i) + &
        w101 * solid_array(i_sol)%eta_table(i2r, i1e, i2i) + &
        w011 * solid_array(i_sol)%eta_table(i1r, i2e, i2i) + &
        w111 * solid_array(i_sol)%eta_table(i2r, i2e, i2i)

  END FUNCTION interp_trilinear



  SUBROUTINE find_bracket_indices(x_in, x, i1, i2, fx)

    ! Bisection search to find adjacent bracket indices i1, i2 in sorted
    ! array x such that x(i1) <= x_in <= x(i2).  Returns the linear fraction
    ! fx = (x_in - x(i1)) / (x(i2) - x(i1)).  Values outside the array range
    ! are clamped to the nearest boundary with a one-time rank-0 warning.

    REAL(num), INTENT(IN) :: x_in
    REAL(num), INTENT(IN) :: x(:)
    INTEGER, INTENT(OUT) :: i1, i2
    REAL(num), INTENT(OUT) :: fx

    INTEGER :: nx
    REAL(num) :: xdif1, xdif2, xdifm
    INTEGER :: im
    LOGICAL, SAVE :: warning = .TRUE.

    nx = SIZE(x)
    xdif1 = x(1) - x_in
    xdif2 = x(nx) - x_in

    IF (xdif1 * xdif2 < 0.0_num) THEN
      i1 = 1
      i2 = nx
      DO
        im = (i1 + i2) / 2
        xdifm = x(im) - x_in
        IF (xdif1 * xdifm < 0.0_num) THEN
          i2 = im
        ELSE
          i1 = im
          xdif1 = xdifm
        END IF
        IF (i2 - i1 == 1) EXIT
      END DO
      fx = (x_in - x(i1)) / (x(i2) - x(i1))
    ELSE
      IF (warning .AND. rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Resistivity table lookup out of range. Clamping to boundary.'
          WRITE(io,*) 'No further warnings will be issued.'
        END DO
        CALL abort_code(c_err_io_error)
        warning = .FALSE.
      END IF
      IF (xdif1 >= 0.0_num) THEN
        i1 = 1
        i2 = MIN(2, nx)
        fx = 0.0_num
      ELSE
        i1 = MAX(nx - 1, 1)
        i2 = nx
        fx = 1.0_num
      END IF
    END IF

  END SUBROUTINE find_bracket_indices

  FUNCTION thomas_fermi_ionisation(ix, iy, iz)

    ! Calculates the average ionisation state of the background ions using the
    ! algorithm outlined in table IV, More, R. M. (1985). "Pressure ionization,
    ! resonances, and the continuity of bound and free states". In Advances in
    ! atomic and molecular physics (Vol. 21, pp. 305-356)
    !
    ! For ease of comparison, we use their variable names with the prefix tf_
    !
    ! This subroutine uses an average global background, where solid_array
    ! properties have been averaged

    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: tf_t0, tf_tf
    REAL(num) :: tf_a, tf_b, tf_c
    REAL(num) :: tf_q1, tf_q
    REAL(num) :: tf_x
    REAL(num) :: thomas_fermi_ionisation

    tf_t0 = hy_te(ix,iy,iz) * kelvin_to_ev * &
        ion_z_avg(ix,iy,iz)**(-4.0_num/3.0_num)
    tf_tf = tf_t0/(1.0_num + tf_t0)

    tf_a = 0.003323_num*tf_t0**0.9718_num + 9.26148e-5_num*tf_t0**3.10165_num
    tf_b = -EXP(-1.7630_num + 1.43175_num*tf_tf + 0.31546_num*tf_tf**7)
    tf_c = -0.366667_num*tf_tf + 0.983333_num

    tf_q1 = tf_a * ion_reduced_density(ix,iy,iz)**tf_b
    tf_q = (ion_reduced_density(ix,iy,iz)**tf_c + tf_q1**tf_c)**(1.0_num/tf_c)

    tf_x = 14.3139_num * tf_q**0.6624_num
    thomas_fermi_ionisation = ion_z_avg(ix,iy,iz)*tf_x &
        / (1.0_num + tf_x + SQRT(1.0_num + 2.0_num*tf_x))

  END FUNCTION thomas_fermi_ionisation



  FUNCTION lee_more_coulomb_log(ix, iy, iz)

    ! Calculates the coulomb logarithm using the method outlined in Lee, More
    ! (1984). Phys. Fluids. 27(5). Note that this paper is written in Gaussian
    ! CGS units, whereas EPOCH is in SI

    INTEGER, INTENT(IN) :: ix, iy, iz
    REAL(num) :: rlm_z_star, rlm_ti, rlm_te, rlm_ni
    REAL(num) :: fermi_temp, ve, ion_term, idebye_huckel2
    REAL(num) :: ion_sphere_rad, debye_huckel, class_impact, uncert_lim
    REAL(num) :: b_min, b_max
    REAL(num) :: lee_more_coulomb_log

    ! Copy out array variables
    rlm_z_star = ion_charge(ix,iy,iz)
    rlm_te = hy_te(ix,iy,iz)
    rlm_ni = ion_ni(ix,iy,iz)

    IF (use_ion_temp) THEN
      rlm_ti = hy_ti(ix,iy,iz)
    ELSE
      rlm_ti = rlm_te
    END IF

    ! Ion sphere radius (interatomic distance)
    ion_sphere_rad = (4.0_num*pi*rlm_ni/3.0_num)**(-1.0_num/3.0_num)

    ! Fermi temperature
    fermi_temp = h_bar**2/(2.0_num*m0*q0*kb) &
        * (3.0_num*pi*rlm_z_star*rlm_ni)**(2.0_num/3.0_num)

    ! Debye-Huckel screening length (Lee-More eq. 19)
    IF (rlm_ti < c_tiny) THEN
      ion_term = 0.0_num
    ELSE
      ion_term = rlm_z_star/rlm_ti
    END IF
    idebye_huckel2 = rlm_z_star*rlm_ni*q0**2/(epsilon0*kb) &
        * (1.0_num/(kb*SQRT(rlm_te**2 + fermi_temp**2)) + ion_term)
    debye_huckel = 1.0_num/SQRT(idebye_huckel2)

    ! Thermal electron velocity (see between Lee-More eq. 20 & 21)
    ve = SQRT(3.0_num*kb*rlm_te/m0)

    ! Classical distance of closest approach (Lee-More eq. 20)
    class_impact = rlm_z_star*q0**2/(4.0_num*pi*epsilon0*m0*ve**2)

    ! High energy b_min limit is uncertainty principle (Lee-More eq. 21)
    uncert_lim = h_planck/(2.0_num*m0*ve)

    ! Impact parameters
    ! Minimum (Lee-More eq. 22)
    b_min = MAX(class_impact, uncert_lim)
    ! Maximum (see between Lee-More eq. 19 and 20 - Debye-Huckel screening
    ! breaks down when less than the interatomic distance)
    b_max = MAX(debye_huckel, ion_sphere_rad)

    ! Coulomb logarithm (Lee-More eq. 17 and following discussion)
    lee_more_coulomb_log = MAX(0.5_num*LOG(1.0_num + (b_max/b_min)**2), 2.0_num)

  END FUNCTION lee_more_coulomb_log

#endif
END MODULE hy_resistivity
