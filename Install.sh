## Install/unInstall package classes in LAMMPS
#  Copies files from here to src/ or deletes those files in src/
#  Recommended:
#   - Install this package "make yes-USER-GFMD" after all other packages
#   - Dont change the install flags while USER-GFMD package is installed


# mode = 0/1/2 for uninstall/install/update
mode=$1


echo "#define REV \"\$Id: REV `svnversion` \$\"" > REV


## Hard-coded install flags (0/1) specify which optional aspects of USER-GFMD to install.
# These should NOT BE CHANGED while USER-GFMD is installed.  Type make no-USER-GFMD first.
# (The concern: "make yes", then user turns off a flag, then later "make no" does not
#               revert the patch for that component.)
usercuda=0
fftw3=0
ncvar=0
manybody=0


## Four small functions are defined in this script, and then called below 
# add_or_rm()
# dependency_check()
# toggle_patch()
# set_or_unset()


# function: add_or_rm arg1 [arg2]
#    a wrapper for either 'rm' (if mode is 0) or 'cp' (if mode is 1 or 2)
# arg1 = file, arg2 = path to LAMMPS src/ dir
add_or_rm () {

  file=`basename $1`
  srcdir=..
  if (test -n "$2") then 
    srcdir=$2 
  fi

  # Mode 0
  if (test $mode = 0) then
    if ((test -e $srcdir/$file) && (! cmp -s $1 $srcdir/$file)) then # These four lines help developers
      echo " Almost removed a file that has changed. Moved instead to $file.del"
      mv $srcdir/$file $srcdir/$file.del
    fi
    rm -f $srcdir/$file
  fi

  # Mode 1/2: Copy into src dir if identical file not already there 
  if (test $mode = 1 || test $mode = 2) then
    if (! cmp -s $1 $srcdir/$file) then
      cp $1 $srcdir -p
      if (test $mode = 2) then
        echo "  Updating src/$file"
      fi
    fi
  fi
}


# function: dependency_check arg1 arg2 arg3
#   change mode from 1 to -1 (or 2 to 0) if the file arg1 does not exist
#   arg1 = name of file that should exist to proceed, arg2 = arg3 = strings for error message
dependency_check () {

  if (! test -e $1) then
    echo " Note: $2=1 in USER-GFMD/Install.sh, but dependency '"$3"' not installed ($1 not found)"

    if (test $mode = 1) then
      echo " so GFMD+"$2" will not be installed."
      mode=-1

    elif (test $mode = 2) then
      echo " so any GFMD+"$2" files will be uninstalled."
      mode=0
    fi
  fi
}

# function: toggle_patch arg1 arg2
#    calls either patch -R (if mode is 0) or patch (if mode is 1 or 2)
# arg1 = filetopatch, arg2 = patchfile
toggle_patch () {

  if ((test ! $mode = -1) && (test ! -e $1)) then
    echo "  Warning: USER-GFMD/Install.sh flag: Expected to un/patch file $1 which does not exist."

  else 
    case $mode in 
    0)
      patch -R $1 $2
      ;;
    1)
      patch $1 $2
      ;;
    2) 
      echo "(If a warning about previously-applied patch prints next, this file was"
      echo " already updated.  Do [n]ot assume reversed, do [n]ot apply anyway.)"
      patch $1 $2
      ;;
    esac
  fi
}


# function: set_or_unset arg1
#   adds (or removes) "-Darg1" to the PKG_INC variable in Makefile.package if mode = 1 (or 0)
#   This is standard in LAMMPS (eg USER-OMP/Install.sh) to control #ifdef statements in the .cpp files
# arg1 = macro name to define in Makefile.package for the compiler preprocessor
set_or_unset () {

  if (test -e ../Makefile.package) then

    if (test $mode = 1) then
      sed -i -e 's/[^ \t]*'$1'[^ \t]* //' ../Makefile.package
      sed -i -e 's|^PKG_INC =[ \t]*|&-D'$1' |' ../Makefile.package

    elif (test $mode = 0) then
      sed -i -e 's/[^ \t]*'$1'[^ \t]* //' ../Makefile.package
    fi
  fi
}



# Add/Remove .cpp files to/from src directory
add_or_rm mathutils/linearalgebra.cpp
add_or_rm stiffnesskernels/geometry.cpp
add_or_rm stiffnesskernels/li_berger.cpp
add_or_rm mathutils/table2d.cpp
add_or_rm atom_vec_gfmd.cpp
add_or_rm gfmd_solver.cpp
add_or_rm solvers/gfmd_solver_fft.cpp
add_or_rm solvers/gfmd_solver_dynamic.cpp
add_or_rm solvers/gfmd_solver_static.cpp  
add_or_rm surface_stiffness.cpp
add_or_rm stiffnesskernels/debug_stiffness.cpp
add_or_rm stiffnesskernels/isotropic_stiffness.cpp
add_or_rm stiffnesskernels/sc100_stiffness.cpp
add_or_rm stiffnesskernels/fcc100_stiffness.cpp
add_or_rm stiffnesskernels/fcc100ft_stiffness.cpp
add_or_rm stiffnesskernels/ft_stiffness.cpp
add_or_rm crystal_surface.cpp
add_or_rm surfaces/dia100_surface.cpp
add_or_rm surfaces/dia111_surface.cpp
add_or_rm surfaces/fcc100_surface.cpp
add_or_rm surfaces/fcc111_surface.cpp
add_or_rm surfaces/sc100_surface.cpp
add_or_rm force_constants.cpp
add_or_rm forceconstants/fc_lj_cut.cpp
add_or_rm forceconstants/fc_lj_cut_fd.cpp
add_or_rm forceconstants/fc_lj_smooth.cpp
add_or_rm forceconstants/fc_pair_potential.cpp
add_or_rm forceconstants/fc_finite_differences.cpp
add_or_rm forceconstants/pair_lj_cut_gf.cpp
add_or_rm gfmd_misc.cpp
add_or_rm fix_gfmd.cpp
add_or_rm extras/fix_contact_rigid.cpp
add_or_rm extras/fix_contact_sphere.cpp
add_or_rm extras/fix_wall_map.cpp
# add_or_rm extras/min_cg_gfmd.cpp These dont compile for me.  Error is
# add_or_rm extras/min_tr_gfmd.cpp  line (99): error: class "LAMMPS_NS::FixGFMD" has no member "prec_gradient"
add_or_rm extras/displace_xeq.cpp
add_or_rm extras/fix_addforce_dynamic.cpp

# Add/Remove .h files from src/ directory
add_or_rm mathutils/complexcomp.h
add_or_rm mathutils/vec.h
add_or_rm mathutils/mat.h
add_or_rm mathutils/linearalgebra.h
add_or_rm stiffnesskernels/geometry.h
add_or_rm stiffnesskernels/li_berger.h
add_or_rm mathutils/table2d.h
add_or_rm atom_vec_gfmd.h
add_or_rm gfmd_solver.h
add_or_rm solvers/gfmd_solver_fft.h
add_or_rm solvers/gfmd_solver_dynamic.h
add_or_rm solvers/gfmd_solver_static.h
add_or_rm surface_stiffness.h
add_or_rm stiffnesskernels/debug_stiffness.h
add_or_rm stiffnesskernels/isotropic_stiffness.h
add_or_rm stiffnesskernels/sc100_stiffness.h
add_or_rm stiffnesskernels/fcc100_stiffness.h
add_or_rm stiffnesskernels/fcc100ft_stiffness.h
add_or_rm stiffnesskernels/ft_stiffness.h
add_or_rm crystal_surface.h
add_or_rm surfaces/dia100_surface.h
add_or_rm surfaces/dia111_surface.h
add_or_rm surfaces/fcc100_surface.h
add_or_rm surfaces/fcc111_surface.h
add_or_rm surfaces/sc100_surface.h
add_or_rm force_constants.h
add_or_rm forceconstants/fc_lj_cut.h
add_or_rm forceconstants/fc_lj_cut_fd.h
add_or_rm forceconstants/fc_lj_smooth.h
add_or_rm forceconstants/fc_pair_potential.h
add_or_rm forceconstants/fc_finite_differences.h
add_or_rm forceconstants/pair_lj_cut_gf.h
add_or_rm gfmd_misc.h
add_or_rm fix_gfmd.h
add_or_rm gfmd_grid.h
add_or_rm extras/fix_contact_rigid.h
add_or_rm extras/fix_contact_sphere.h
add_or_rm extras/fix_wall_map.h
# add_or_rm folder/min_cg_gfmd.h These dont compile for me.  Error is 
# add_or_rm folder/min_tr_gfmd.h line (99): error: class "LAMMPS_NS::FixGFMD" has no member "prec_gradient"
add_or_rm extras/displace_xeq.h
add_or_rm extras/fix_addforce_dynamic.h

# Apply the tiny patches that LAMMPS requires for new styles
toggle_patch ../atom.cpp atom.cpp.patch
toggle_patch ../atom.h atom.h.patch
toggle_patch ../pair_lj_smooth.h forceconstants/pair_lj_smooth.h.patch

if (test -e ../comm_brick.cpp) then
  # We have a new LAMMPS version. Patch comm_brick.cpp
  toggle_patch ../comm_brick.cpp comm_brick.cpp.patch
else
  toggle_patch ../comm.cpp comm.cpp.patch
fi

if ((test $mode -ge 1) && (test ! -e ../fft3d_wrap.h)) then
  echo " Warning: fft3d_wrap.h not found in src/, but GFMD requires 'make yes-KSPACE' before compile"
fi

#######################################################
# Optional aspects of USER-GFMD to install/ uninstall #
#######################################################

if (test $ncvar = 1) then 
  add_or_rm extras/gfmd_analyzer.cpp
  add_or_rm extras/gfmd_analyzer.h
fi


if (test $manybody = 1) then
  echo "Many body functionality not included"
fi


# 
if (test $usercuda = 1) then
  echo "Many body functionality not included"
fi

#
if (test $fftw3 = 1) then 

  set_or_unset GFMD_FFTW3

  add_or_rm stiffnesskernels/nonperiodic_stiffness.cpp
  add_or_rm stiffnesskernels/nonperiodic_stiffness.h
fi


