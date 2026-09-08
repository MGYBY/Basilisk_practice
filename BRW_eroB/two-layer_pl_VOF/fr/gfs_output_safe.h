/**
# Safe Gerris/GfsView output for the uploaded Basilisk TREE+MPI snapshot

The uploaded Basilisk snapshot contains an off-by-one inconsistency between
`grid/tree-mpi.h::z_indexing()` and `output.h::output_gfs()`: the former
returns the number of indexed cells, while the latter treats the return
value as the largest zero-based index and adds one.  The resulting MPI GFS
file contains one zero-filled cell record immediately before the final
closing brace.  Gerris/GfsView then reports e.g.

    expecting a closing brace

This wrapper leaves the upstream writer unchanged and removes *only* a
trailing, entirely zero-filled record.  A valid final tree record cannot be
all zero because it must carry the Gerris leaf flag.  The repair is done by
stream-copying to a temporary file, so no POSIX `ftruncate()` dependency is
introduced into Basilisk C.
*/
#ifndef BASILISK_GFS_OUTPUT_SAFE_H
#define BASILISK_GFS_OUTPUT_SAFE_H

static size_t gfs_binary_cell_size (scalar * fields)
{
  size_t nfields = 0;
  for (scalar s in fields)
    if (s.name)
      nfields++;

  /* flags + solid marker + one double per named variable */
  return sizeof(unsigned) + sizeof(double)*(1u + nfields);
}

static void gfs_copy_prefix_and_close (const char * filename,
                                       long prefix_length)
{
  const size_t chunk_size = 1024u*1024u;
  unsigned char * buffer = (unsigned char *) malloc (chunk_size);
  char * temporary = (char *) malloc (strlen(filename) + 16u);
  if (!buffer || !temporary) {
    fprintf (ferr, "gfs-output-safe: memory allocation failed.\n");
    free (buffer);
    free (temporary);
    exit (1);
  }

  sprintf (temporary, "%s.tmp-gfs", filename);
  FILE * input = fopen (filename, "rb");
  FILE * output = fopen (temporary, "wb");
  if (!input || !output) {
    fprintf (ferr, "gfs-output-safe: cannot create repaired '%s'.\n",
             filename);
    if (input) fclose (input);
    if (output) fclose (output);
    free (buffer);
    free (temporary);
    exit (1);
  }

  long remaining = prefix_length;
  while (remaining > 0) {
    const size_t wanted = remaining < (long)chunk_size ?
      (size_t)remaining : chunk_size;
    const size_t got = fread (buffer, 1, wanted, input);
    if (got != wanted || fwrite (buffer, 1, got, output) != got) {
      fprintf (ferr, "gfs-output-safe: failed while copying '%s'.\n",
               filename);
      fclose (input);
      fclose (output);
      remove (temporary);
      free (buffer);
      free (temporary);
      exit (1);
    }
    remaining -= (long)got;
  }

  if (fwrite ("}\n", 1, 2, output) != 2 ||
      fclose (input) != 0 || fclose (output) != 0 ||
      rename (temporary, filename) != 0) {
    fprintf (ferr, "gfs-output-safe: cannot replace '%s'.\n", filename);
    remove (temporary);
    free (buffer);
    free (temporary);
    exit (1);
  }

  free (buffer);
  free (temporary);
}

static bool gfs_remove_trailing_zero_record (const char * filename,
                                             scalar * fields)
{
#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  bool repaired = false;
  if (pid() == 0) {
    FILE * fp = fopen (filename, "rb");
    if (!fp) {
      fprintf (ferr, "gfs-output-safe: cannot open '%s'.\n", filename);
      exit (1);
    }

    if (fseek (fp, 0, SEEK_END) != 0) {
      fprintf (ferr, "gfs-output-safe: cannot seek in '%s'.\n", filename);
      fclose (fp);
      exit (1);
    }

    const long end = ftell (fp);
    const size_t record_size = gfs_binary_cell_size (fields);
    if (end < (long)(record_size + 2u)) {
      fprintf (ferr, "gfs-output-safe: '%s' is unexpectedly short.\n",
               filename);
      fclose (fp);
      exit (1);
    }

    unsigned char closing[2];
    if (fseek (fp, end - 2L, SEEK_SET) != 0 ||
        fread (closing, 1, 2, fp) != 2 ||
        closing[0] != '}' || closing[1] != '\n') {
      fprintf (ferr,
               "gfs-output-safe: '%s' does not end with the expected '}\\n'.\n",
               filename);
      fclose (fp);
      exit (1);
    }

    unsigned char * record = (unsigned char *) malloc (record_size);
    if (!record) {
      fprintf (ferr, "gfs-output-safe: memory allocation failed.\n");
      fclose (fp);
      exit (1);
    }

    const long record_start = end - 2L - (long)record_size;
    if (fseek (fp, record_start, SEEK_SET) != 0 ||
        fread (record, 1, record_size, fp) != record_size) {
      fprintf (ferr, "gfs-output-safe: cannot inspect '%s'.\n", filename);
      free (record);
      fclose (fp);
      exit (1);
    }

    bool all_zero = true;
    for (size_t j = 0; j < record_size; j++)
      if (record[j] != 0u) {
        all_zero = false;
        break;
      }
    free (record);
    fclose (fp);

    if (all_zero) {
      gfs_copy_prefix_and_close (filename, record_start);
      repaired = true;
      fprintf (ferr,
               "# gfs-output-safe: removed one %zu-byte MPI padding record "
               "from %s\n",
               record_size, filename);
    }
  }

#if _MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
  return repaired;
}

static void output_gfs_safe (const char * filename,
                             scalar * fields,
                             bool translate)
{
  output_gfs (file = (char *) filename,
              list = fields,
              translate = translate);
  gfs_remove_trailing_zero_record (filename, fields);
}

#endif
