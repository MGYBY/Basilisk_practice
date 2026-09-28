/* Native C only. Do NOT #include this file from a Basilisk C source.
 * The GFS wrapper has already verified a zero padding record and synchronized
 * the MPI writers. File contents and the in-place repair algorithm are
 * unchanged from v2p6/v2p7; only the compilation boundary has changed.
 * qcc sees the small public prototype, never these POSIX headers.
 */
#ifndef _POSIX_C_SOURCE
# define _POSIX_C_SOURCE 200809L
#endif
#ifndef _FILE_OFFSET_BITS
# define _FILE_OFFSET_BITS 64
#endif
#include "gfs_tail_truncate.h"
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

int gfs_tail_truncate (const char * filename, long record_start)
{
  if (!filename || record_start < 0) {
    errno = EINVAL;
    return -1;
  }
  FILE * fp = fopen (filename, "r+b");
  if (!fp)
    return -1;
  int failed = 0;
  if (fseek (fp, record_start, SEEK_SET) != 0 ||
      fwrite ("}\n", 1, 2, fp) != 2 || fflush (fp) != 0)
    failed = 1;
  if (!failed && ftruncate (fileno(fp), (off_t)record_start + 2) != 0)
    failed = 1;
  const int saved_errno = errno;
  if (fclose (fp) != 0)
    return -1;
  if (failed) {
    errno = saved_errno;
    return -1;
  }
  return 0;
}
