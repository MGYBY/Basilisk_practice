/* Public interface only: this header is parsed by qcc.
 * Keep POSIX/system headers and the implementation in gfs_tail_truncate.c.
 * Build that ordinary-C file with the native compiler and link its object.
 * Shared build fix for top-solid v2p6a and no-top-solid v2p7a.
 */
#ifndef FRONT_RUNNER_GFS_TAIL_TRUNCATE_H
#define FRONT_RUNNER_GFS_TAIL_TRUNCATE_H
extern int gfs_tail_truncate (const char * filename, long record_start);
#endif
