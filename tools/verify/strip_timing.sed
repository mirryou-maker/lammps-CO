# Drop everything that legitimately varies between two runs of the same physics:
# wall-clock numbers, the per-section timing table, and host/build banners.
/CPU = /d
/^[A-Za-z][A-Za-z ]*| *[0-9]/d
/^Other  */d
/^LAMMPS /d
/^OMP_NUM/d
/^  using /d
/^Total wall/d
/^Loop time/d
/^Performance/d
/^Section /d
/^Nlocal:/d
/^Nghost:/d
/^Neighs:/d
/^Ave /d
/^Dangerous/d
/^WARNING/d
/^Per MPI rank memory/d
/[0-9]\.[0-9]* *% /d
