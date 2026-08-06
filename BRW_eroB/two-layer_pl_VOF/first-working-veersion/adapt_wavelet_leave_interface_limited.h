/**
# Wavelet adaptation with protected VOF and embedded interfaces

This is a current-Basilisk reimplementation of the
`adapt_wavelet_leave_interface_limited.h` strategy used in
MGYBY/Basilisk_practice and originally attributed there to César
Pairetti.  It follows the current `adapt_wavelet()` implementation in
`grid/tree-common.h`, including block-field handling and dependency-aware
field refinement, and adds:

* a spatially varying maximum level supplied by `MLFun(x,y,z)`;
* forced refinement of mixed cells of every field in `vol_frac`;
* optional same-level padding around protected interfaces.

The header must be included after `grid/quadtree.h`.
*/

#ifndef BASILISK_ADAPT_WAVELET_LEAVE_INTERFACE_LIMITED_H
#define BASILISK_ADAPT_WAVELET_LEAVE_INTERFACE_LIMITED_H

#if !TREE
# error "adapt_wavelet_leave_interface_limited.h requires a tree grid"
#endif

typedef struct Adapt_limited {
  scalar * slist;       // fields used for the wavelet error estimate
  scalar * vol_frac;    // VOF/embedded fractions whose mixed cells are protected
  double * max;         // one tolerance per entry of slist
  int (*MLFun)(double, double, double); // legacy centre-based local maximum
  int (*MLFunDelta)(double, double, double, double); // optional cell-size-aware maximum
  int (*MinFun)(double, double, double); // optional local minimum level
  int minlevel;         // global minimum level
  int padding;          // same-level neighbour padding around mixed cells
  double interface_eps; // mixed-cell test: eps < c < 1-eps
  scalar * list;        // fields to update; defaults to all
} Adapt_limited;

static inline int awl_local_minlevel (const Adapt_limited * p,
                                      double xx, double yy, double zz)
{
  int lev = p->MinFun ? p->MinFun (xx, yy, zz) : p->minlevel;
  return max(1, lev);
}

static inline int awl_local_maxlevel (const Adapt_limited * p,
                                      double xx, double yy, double zz,
                                      double delta)
{
  int lev = p->MLFunDelta ? p->MLFunDelta (xx, yy, zz, delta) :
    (p->MLFun ? p->MLFun (xx, yy, zz) : depth());
  const int minlev = awl_local_minlevel (p, xx, yy, zz);
  if (lev < minlev)
    lev = minlev;
  return lev;
}

trace
astats adapt_wavelet_limited (struct Adapt_limited p)
{
  scalar * ilist = p.list;
  scalar * list = p.list;

  if (p.minlevel < 1)
    p.minlevel = 1;
  if (p.padding < 0)
    p.padding = 0;
  if (p.interface_eps <= 0.)
    p.interface_eps = 1.e-6;

  if (is_constant(cm)) {
    if (list == NULL || list == all)
      list = list_copy (all);
    boundary (list);
    restriction (p.slist);
  }
  else {
    if (list == NULL || list == all) {
      list = list_copy ({cm, fm});
      for (scalar s in all)
        list = list_add (list, s);
    }
    boundary (list);
    scalar * listr = list_concat (p.slist, {cm});
    restriction (listr);
    free (listr);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  for (scalar s in list)
    listc = list_add_depend (listc, s);

  tree->refined.n = 0;
  static const int refined = 1 << user;
  static const int too_fine = 1 << (user + 1);
  static const int too_coarse = 1 << (user + 2);
  static const int just_fine = 1 << (user + 3);

  foreach_cell() {
    if (!is_active(cell))
      continue;

    if (is_leaf(cell)) {
      if (cell.flags & too_coarse) {
        cell.flags &= ~too_coarse;
        refine_cell (point, listc, refined, &tree->refined);
        st.nf++;
      }
      continue;
    }

    if (cell.flags & refined) {
      cell.flags &= ~too_coarse;
      continue;
    }

    bool local = is_local(cell);
    if (!local)
      foreach_child()
        if (is_local(cell)) {
          local = true;
          break;
        }
    if (!local)
      continue;

    int i = 0;
    for (scalar s in p.slist) {
      const double emax = p.max[i++];
      double sc[(1 << dimension)*s.block];
      double * b = sc;
      foreach_child()
        foreach_blockf(s)
          *b++ = s[];

      s.prolongation (point, s);
      b = sc;
      foreach_child()
        foreach_blockf(s) {
          const int localmin = awl_local_minlevel (&p, x, y, z);
          const int localmax = awl_local_maxlevel (&p, x, y, z, Delta);
          const double e = fabs(*b - s[]);
          if (level < localmin) {
            cell.flags &= ~too_fine;
            cell.flags |= too_coarse;
          }
          else if (e > emax && level < localmax) {
            cell.flags &= ~too_fine;
            cell.flags |= too_coarse;
          }
          else if ((e <= emax/1.5 || level > localmax) &&
                   !(cell.flags & (too_coarse | just_fine))) {
            if (level >= p.minlevel)
              cell.flags |= too_fine;
          }
          else if (!(cell.flags & too_coarse)) {
            cell.flags &= ~too_fine;
            cell.flags |= just_fine;
          }
          s[] = *b++;
        }
    }

    /* Protect mixed VOF/embedded-fraction cells independently of the
       wavelet error.  A protected cell at its allowed maximum level is
       kept but not refined further. */
    foreach_child() {
      bool mixed = false;
      for (scalar vf in p.vol_frac) {
        const double c = vf[];
        if (c > p.interface_eps && c < 1. - p.interface_eps) {
          mixed = true;
          break;
        }
      }

      if (mixed) {
        const int localmax = awl_local_maxlevel (&p, x, y, z, Delta);
        cell.flags &= ~too_fine;
        if (level < localmax) {
          cell.flags |= too_coarse;
          cell.flags &= ~just_fine;
        }

        if (p.padding > 0)
          foreach_neighbor (p.padding) {
            const int neighbour_max = awl_local_maxlevel (&p, x, y, z, Delta);
            cell.flags &= ~too_fine;
            if (level < neighbour_max) {
              cell.flags |= too_coarse;
              cell.flags &= ~just_fine;
            }
          }
      }
    }

    foreach_child() {
      cell.flags &= ~just_fine;
      if (!is_leaf(cell)) {
        cell.flags &= ~too_coarse;
        const int localmin = awl_local_minlevel (&p, x, y, z);
        const int localmax = awl_local_maxlevel (&p, x, y, z, Delta);
        if (level < localmin) {
          cell.flags &= ~too_fine;
          cell.flags |= too_coarse;
        }
        else if (level >= localmax)
          cell.flags |= too_fine;
      }
      else if (!is_active(cell))
        cell.flags &= ~too_coarse;
    }
  }

  mpi_boundary_refine (listc);

  /* Coarsening pass, including the standard 2:1-symmetry cleanup. */
  for (int l = depth(); l >= 0; l--) {
    foreach_cell()
      if (!is_boundary(cell)) {
        if (level == l) {
          if (!is_leaf(cell)) {
            if (cell.flags & refined)
              cell.flags &= ~(refined | too_fine);
            else if (cell.flags & too_fine) {
              if (is_local(cell) && coarsen_cell (point, listc))
                st.nc++;
              cell.flags &= ~too_fine;
            }
          }
          if (cell.flags & too_fine)
            cell.flags &= ~too_fine;
          else if (level > 0 && (aparent(0).flags & too_fine))
            aparent(0).flags &= ~too_fine;
          continue;
        }
        else if (is_leaf(cell))
          continue;
      }
    mpi_boundary_coarsen (l, too_fine);
  }

  free (listc);
  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (list);

  if (list != ilist)
    free (list);

  return st;
}

#endif
