/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reax_c.h"
#include "error.h"
#include "reaxc_ffield.h"
#include "reaxc_tool_box.h"

typedef union ff_entry_t {
    int i;
    real r;
    char *s;
} ff_entry;

typedef struct ff_reader_t {
    unsigned int (*read_record)(struct ff_reader_t *ctxt, ff_entry_t *record);
    int (*get_int)(struct ff_reader_t *ctxt, ff_entry_t *entry);
    real (*get_real)(struct ff_reader_t *ctxt, ff_entry_t *entry);
    char *(*get_string)(struct ff_reader_t *ctxt, ff_entry_t *entry);
    void (*destroy)(struct ff_reader_t *ctxt);
    FILE *file;
    void *private_data;
} ff_reader;

typedef struct ff_text_reader_data_t {
    char    *s;
    char   **tmp;
} ff_text_reader_data;

static unsigned int text_reader_read_record(struct ff_reader_t *ctxt, ff_entry_t *record)
{
    ff_text_reader_data *d = (ff_text_reader_data *) ctxt->private_data;
    unsigned int n, i;

    fgets(d->s,MAX_LINE,ctxt->file);
    n = Tokenize(d->s,&(d->tmp));

    for (i = 0; i < n; i++) {
        record[i].s = d->tmp[i];
    }

    return n;
}

static int text_reader_get_int(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return atoi(entry->s);
}

static real text_reader_get_real(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return atof(entry->s);
}

static char *text_reader_get_string(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return entry->s;
}

static void text_reader_destroy(struct ff_reader_t *ctxt)
{
    ff_text_reader_data *d = (ff_text_reader_data *) ctxt->private_data;
    unsigned int i;

    for( i = 0; i < MAX_TOKENS; i++ )
        free( d->tmp[i] );
    free( d->tmp );
    free( d->s );
}

static void init_text_reader(ff_reader *ctxt, FILE *fp)
{
    ff_text_reader_data *d;
    unsigned int i;

    d = (ff_text_reader_data *) malloc(sizeof(ff_text_reader_data));

    d->s = (char*) malloc(sizeof(char)*MAX_LINE);
    d->tmp = (char**) malloc(sizeof(char*)*MAX_TOKENS);
    for (i=0; i < MAX_TOKENS; i++) {
        d->tmp[i] = (char*) malloc(sizeof(char)*MAX_TOKEN_LEN);
    }

    ctxt->read_record = text_reader_read_record;
    ctxt->get_int = text_reader_get_int;
    ctxt->get_real = text_reader_get_real;
    ctxt->get_string = text_reader_get_string;
    ctxt->destroy = text_reader_destroy;

    ctxt->file = fp;
    ctxt->private_data = (void *) d;
}

typedef struct ff_binary_reader_data_t {
    char s[MAX_LINE];
    unsigned int spos;
} ff_binary_reader_data;

static unsigned int binary_reader_read_record(struct ff_reader_t *ctxt, ff_entry_t *record)
{
    ff_binary_reader_data *d = (ff_binary_reader_data *) ctxt->private_data;
    unsigned int n, i, len;
    unsigned char type;
    int ival;
    float fval;
    double dval;

    d->spos = 0;
    d->s[0] = '\0';
    fread(&n, sizeof(n), 1, ctxt->file);

    for (i = 0; i < n; i++) {
        fread(&type, sizeof(type), 1, ctxt->file);
        switch (type) {
        case 'I':
            fread(&ival, sizeof(ival), 1, ctxt->file);
            record[i].i = ival;
            break;
        case 'F':
            fread(&fval, sizeof(fval), 1, ctxt->file);
            record[i].r = fval;
            break;
        case 'D':
            fread(&dval, sizeof(dval), 1, ctxt->file);
            record[i].r = dval;
            break;
        case 'S':
            fread(&len, sizeof(len), 1, ctxt->file);
            if (len > 0) {
                d->spos++;
                record[i].s = d->s + d->spos;
                fread(d->s + d->spos, len, 1, ctxt->file);
                d->spos += len;
                d->s[d->spos] = '\0';
            } else {
                record[i].s = d->s + d->spos;
            }
            break;
        }
    }

    return n;
}

static int binary_reader_get_int(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return entry->i;
}

static real binary_reader_get_real(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return entry->r;
}

static char *binary_reader_get_string(struct ff_reader_t *ctxt, ff_entry_t *entry)
{
    return entry->s;
}

static void binary_reader_destroy(struct ff_reader_t *ctxt)
{
}


static void init_binary_reader(ff_reader *ctxt, FILE *fp)
{
    ff_binary_reader_data *d;
    unsigned int i;

    d = (ff_binary_reader_data *) malloc(sizeof(ff_binary_reader_data));

    ctxt->read_record = binary_reader_read_record;
    ctxt->get_int = binary_reader_get_int;
    ctxt->get_real = binary_reader_get_real;
    ctxt->get_string = binary_reader_get_string;
    ctxt->destroy = binary_reader_destroy;

    ctxt->file = fp;
    ctxt->private_data = (void *) d;
}

static int ff_read_header(FILE *fp)
{
    char header[7];
    int c;
    fgets(header, sizeof(header), fp);
    if (strcmp(header, "BINFF\n") == 0) {
        return 1;
    }

    if (strchr(header, '\n') == NULL) {
        do {
            c = getc(fp);
        } while ((c != '\n') && (c != EOF));
    }
    return 0;
}


char Read_Force_Field( char *ffield_file, reax_interaction *reax,
                       control_params *control )
{
  FILE    *fp;
  ff_reader rd;
  ff_entry *rec;
  char *tmp;
  char ****tor_flag;
  int      c, i, j, k, l, m, n, o, p, cnt;
  int lgflag = control->lgflag;
  int errorflag = 1;
  real     val;
  MPI_Comm comm;
  int error_code = 0;

  comm = MPI_COMM_WORLD;

  /* open force field file */
  if ( (fp = fopen( ffield_file, "r" ) ) == NULL ) {
    fprintf( stderr, "error opening the force field file! terminating...\n" );
    error_code = FILE_NOT_FOUND;
    goto end;
  }

  /* reading first header comment */
  if (ff_read_header(fp) == 1) {
    init_binary_reader(&rd, fp);
  } else {
    init_text_reader(&rd, fp);
  }

  rec = (ff_entry *) malloc(MAX_TOKENS * sizeof(ff_entry));

  /* line 2 is number of global parameters */
  c = (*rd.read_record)(&rd, rec);

  /* reading the number of global parameters */
  n = (*rd.get_int)(&rd, &(rec[0]));
  if (n < 1) {
    fprintf( stderr, "WARNING: number of globals in ffield file is 0!\n" );
    goto cleanup_reader;
  }

  reax->gp.n_global = n;
  if (!reax->gp.l) {
    reax->gp.l = (real*) malloc(sizeof(real)*n);
  }
  else {
    memset(reax->gp.l, 0, sizeof(real)*n);
  }

  /* see reax_types.h for mapping between l[i] and the lambdas used in ff */
  for (i=0; i < n; i++) {
    c = (*rd.read_record)(&rd, rec);

    reax->gp.l[i] = (*rd.get_real)(&rd, &(rec[0]));
  }

  control->bo_cut    = 0.01 * reax->gp.l[29];
  control->nonb_low  = reax->gp.l[11];
  control->nonb_cut  = reax->gp.l[12];

  /* next line is number of atom types and some comments */
  c = (*rd.read_record)(&rd, rec);
  reax->num_atom_types = (*rd.get_int)(&rd, &(rec[0]));

  /* 3 lines of comments */
  c = (*rd.read_record)(&rd, rec);
  c = (*rd.read_record)(&rd, rec);
  c = (*rd.read_record)(&rd, rec);

  tor_flag  = (char****)
    scalloc( reax->num_atom_types, sizeof(char***), "tor_flag", comm );
  for( i = 0; i < reax->num_atom_types; i++ ) {
    tor_flag[i]  = (char***)
      scalloc( reax->num_atom_types, sizeof(char**), "tor_flag[i]", comm );

    for( j = 0; j < reax->num_atom_types; j++ ) {
      tor_flag[i][j]  = (char**)
        scalloc( reax->num_atom_types, sizeof(char*), "tor_flag[i,j]", comm );

      for (k=0; k < reax->num_atom_types; k++) {
        tor_flag[i][j][k]  = (char*)
          scalloc( reax->num_atom_types, sizeof(char), "tor_flag[i,j,k]",
                   comm );
      }
    }
  }

  if (!reax->sbp) {
    /* Allocating structures in reax_interaction */
    reax->sbp = (single_body_parameters*)
      scalloc( reax->num_atom_types, sizeof(single_body_parameters), "sbp",
	       comm );
    reax->tbp = (two_body_parameters**)
      scalloc( reax->num_atom_types, sizeof(two_body_parameters*), "tbp", comm );
    reax->thbp= (three_body_header***)
      scalloc( reax->num_atom_types, sizeof(three_body_header**), "thbp", comm );
    reax->hbp = (hbond_parameters***)
      scalloc( reax->num_atom_types, sizeof(hbond_parameters**), "hbp", comm );
    reax->fbp = (four_body_header****)
      scalloc( reax->num_atom_types, sizeof(four_body_header***), "fbp", comm );

    for( i = 0; i < reax->num_atom_types; i++ ) {
      reax->tbp[i] = (two_body_parameters*)
        scalloc( reax->num_atom_types, sizeof(two_body_parameters), "tbp[i]",
	         comm );
      reax->thbp[i]= (three_body_header**)
        scalloc( reax->num_atom_types, sizeof(three_body_header*), "thbp[i]",
	         comm );
      reax->hbp[i] = (hbond_parameters**)
        scalloc( reax->num_atom_types, sizeof(hbond_parameters*), "hbp[i]",
	         comm );
      reax->fbp[i] = (four_body_header***)
        scalloc( reax->num_atom_types, sizeof(four_body_header**), "fbp[i]",
	         comm );

      for( j = 0; j < reax->num_atom_types; j++ ) {
        reax->thbp[i][j]= (three_body_header*)
	  scalloc( reax->num_atom_types, sizeof(three_body_header), "thbp[i,j]",
                   comm );
        reax->hbp[i][j] = (hbond_parameters*)
	  scalloc( reax->num_atom_types, sizeof(hbond_parameters), "hbp[i,j]",
                   comm );
        reax->fbp[i][j] = (four_body_header**)
	  scalloc( reax->num_atom_types, sizeof(four_body_header*), "fbp[i,j]",
		   comm );

        for (k=0; k < reax->num_atom_types; k++) {
	  reax->fbp[i][j][k] = (four_body_header*)
	    scalloc( reax->num_atom_types, sizeof(four_body_header), "fbp[i,j,k]",
		     comm );
        }
      }
    }
  }
  else {
    /* Arrays already allocated, clear them */
    memset( reax->sbp, 0, reax->num_atom_types * sizeof(single_body_parameters) );

    for( i = 0; i < reax->num_atom_types; i++ ) {
      memset( reax->tbp[i], 0, reax->num_atom_types * sizeof(two_body_parameters) );

      for( j = 0; j < reax->num_atom_types; j++ ) {
        memset( reax->thbp[i][j], 0, reax->num_atom_types * sizeof(three_body_header) );
        memset( reax->hbp[i][j], 0, reax->num_atom_types * sizeof(hbond_parameters) );

        for (k=0; k < reax->num_atom_types; k++) {
          memset( reax->fbp[i][j][k], 0, reax->num_atom_types * sizeof(four_body_header) );
        }
      }
    }
  }

  reax->gp.vdw_type = 0;


  for( i = 0; i < reax->num_atom_types; i++ ) {
    /* line one */
    c = (*rd.read_record)(&rd, rec);

    tmp = (*rd.get_string)(&rd, &(rec[0]));
    for( j = 0; j < (int)(strlen(tmp)); ++j )
      reax->sbp[i].name[j] = toupper( tmp[j] );

    val = (*rd.get_real)(&rd, &(rec[1])); reax->sbp[i].r_s        = val;
    val = (*rd.get_real)(&rd, &(rec[2])); reax->sbp[i].valency    = val;
    val = (*rd.get_real)(&rd, &(rec[3])); reax->sbp[i].mass       = val;
    val = (*rd.get_real)(&rd, &(rec[4])); reax->sbp[i].r_vdw      = val;
    val = (*rd.get_real)(&rd, &(rec[5])); reax->sbp[i].epsilon    = val;
    val = (*rd.get_real)(&rd, &(rec[6])); reax->sbp[i].gamma      = val;
    val = (*rd.get_real)(&rd, &(rec[7])); reax->sbp[i].r_pi       = val;
    val = (*rd.get_real)(&rd, &(rec[8])); reax->sbp[i].valency_e  = val;
    reax->sbp[i].nlp_opt = 0.5 * (reax->sbp[i].valency_e-reax->sbp[i].valency);

    /* line two */
    c = (*rd.read_record)(&rd, rec);

    val = (*rd.get_real)(&rd, &(rec[0])); reax->sbp[i].alpha      = val;
    val = (*rd.get_real)(&rd, &(rec[1])); reax->sbp[i].gamma_w    = val;
    val = (*rd.get_real)(&rd, &(rec[2])); reax->sbp[i].valency_boc= val;
    val = (*rd.get_real)(&rd, &(rec[3])); reax->sbp[i].p_ovun5    = val;
    val = (*rd.get_real)(&rd, &(rec[4]));
    val = (*rd.get_real)(&rd, &(rec[5])); reax->sbp[i].chi        = val;
    val = (*rd.get_real)(&rd, &(rec[6])); reax->sbp[i].eta        = 2.0 * val;
    val = (*rd.get_real)(&rd, &(rec[7])); reax->sbp[i].p_hbond = (int) val;

    /* line 3 */
    c = (*rd.read_record)(&rd, rec);

    val = (*rd.get_real)(&rd, &(rec[0])); reax->sbp[i].r_pi_pi    = val;
    val = (*rd.get_real)(&rd, &(rec[1])); reax->sbp[i].p_lp2      = val;
    val = (*rd.get_real)(&rd, &(rec[2]));
    val = (*rd.get_real)(&rd, &(rec[3])); reax->sbp[i].b_o_131    = val;
    val = (*rd.get_real)(&rd, &(rec[4])); reax->sbp[i].b_o_132    = val;
    val = (*rd.get_real)(&rd, &(rec[5])); reax->sbp[i].b_o_133    = val;
    val = (*rd.get_real)(&rd, &(rec[6]));
    val = (*rd.get_real)(&rd, &(rec[7]));

    /* line 4  */
    c = (*rd.read_record)(&rd, rec);

    /* Sanity check */
    if (c < 3) {
      fprintf(stderr, "Inconsistent ffield file (reaxc_ffield.cpp) \n");
      error_code = FILE_NOT_FOUND;
      goto cleanup_tor_flag;
    }

    val = (*rd.get_real)(&rd, &(rec[0])); reax->sbp[i].p_ovun2    = val;
    val = (*rd.get_real)(&rd, &(rec[1])); reax->sbp[i].p_val3     = val;
    val = (*rd.get_real)(&rd, &(rec[2]));
    val = (*rd.get_real)(&rd, &(rec[3])); reax->sbp[i].valency_val= val;
    val = (*rd.get_real)(&rd, &(rec[4])); reax->sbp[i].p_val5     = val;
    val = (*rd.get_real)(&rd, &(rec[5])); reax->sbp[i].rcore2     = val;
    val = (*rd.get_real)(&rd, &(rec[6])); reax->sbp[i].ecore2     = val;
    val = (*rd.get_real)(&rd, &(rec[7])); reax->sbp[i].acore2     = val;

    /* line 5, only if lgvdw is yes */
    if (lgflag) {
      c = (*rd.read_record)(&rd, rec);

      /* Sanity check */
      if (c > 3) {
        fprintf(stderr, "Inconsistent ffield file (reaxc_ffield.cpp) \n");
        error_code = FILE_NOT_FOUND;
        goto cleanup_tor_flag;
      }

      val = (*rd.get_real)(&rd, &(rec[0])); reax->sbp[i].lgcij           = val;
      val = (*rd.get_real)(&rd, &(rec[1])); reax->sbp[i].lgre           = val;
    }

    if( reax->sbp[i].rcore2>0.01 && reax->sbp[i].acore2>0.01 ){ // Inner-wall
      if( reax->sbp[i].gamma_w>0.5 ){ // Shielding vdWaals
        if( reax->gp.vdw_type != 0 && reax->gp.vdw_type != 3 ) {
          if (errorflag)
            fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n"        \
                   "Force field parameters for element %s\n"                \
                   "indicate inner wall+shielding, but earlier\n"        \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   reax->sbp[i].name );
          errorflag = 0;
        }
        else{
          reax->gp.vdw_type = 3;
        }
      }
      else {  // No shielding vdWaals parameters present
        if( reax->gp.vdw_type != 0 && reax->gp.vdw_type != 2 )
          fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n"        \
                   "Force field parameters for element %s\n"                \
                   "indicate inner wall without shielding, but earlier\n" \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   reax->sbp[i].name );
        else{
          reax->gp.vdw_type = 2;
        }
      }
    }
    else{ // No Inner wall parameters present
      if( reax->sbp[i].gamma_w>0.5 ){ // Shielding vdWaals
        if( reax->gp.vdw_type != 0 && reax->gp.vdw_type != 1 )
          fprintf( stderr, "Warning: inconsistent vdWaals-parameters\n"        \
                   "Force field parameters for element %s\n"                \
                   "indicate  shielding without inner wall, but earlier\n" \
                   "atoms indicate different vdWaals-method.\n"                \
                   "This may cause division-by-zero errors.\n"                \
                   "Keeping vdWaals-setting for earlier atoms.\n",
                   reax->sbp[i].name );
        else{
          reax->gp.vdw_type = 1;
        }
      }
      else{
        fprintf( stderr, "Error: inconsistent vdWaals-parameters\n"\
                 "No shielding or inner-wall set for element %s\n",
                 reax->sbp[i].name );
        error_code = INVALID_INPUT;
        goto cleanup_tor_flag;
      }
    }
  }

  /* Equate vval3 to valf for first-row elements (25/10/2004) */
  for( i = 0; i < reax->num_atom_types; i++ )
    if( reax->sbp[i].mass < 21 ){
      reax->sbp[i].valency_val = reax->sbp[i].valency_boc;
    }

  /* next line is number of two body combination and some comments */
  c = (*rd.read_record)(&rd, rec);
  l = (*rd.get_int)(&rd, &(rec[0]));

  /* a line of comments */
  c = (*rd.read_record)(&rd, rec);

  for (i=0; i < l; i++) {
    /* line 1 */
    c = (*rd.read_record)(&rd, rec);

    j = (*rd.get_int)(&rd, &(rec[0])) - 1;
    k = (*rd.get_int)(&rd, &(rec[1])) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types) {

      val = (*rd.get_real)(&rd, &(rec[2])); reax->tbp[j][k].De_s      = val;
      reax->tbp[k][j].De_s      = val;
      val = (*rd.get_real)(&rd, &(rec[3])); reax->tbp[j][k].De_p      = val;
      reax->tbp[k][j].De_p      = val;
      val = (*rd.get_real)(&rd, &(rec[4])); reax->tbp[j][k].De_pp     = val;
      reax->tbp[k][j].De_pp     = val;
      val = (*rd.get_real)(&rd, &(rec[5])); reax->tbp[j][k].p_be1     = val;
      reax->tbp[k][j].p_be1     = val;
      val = (*rd.get_real)(&rd, &(rec[6])); reax->tbp[j][k].p_bo5     = val;
      reax->tbp[k][j].p_bo5     = val;
      val = (*rd.get_real)(&rd, &(rec[7])); reax->tbp[j][k].v13cor    = val;
      reax->tbp[k][j].v13cor    = val;

      val = (*rd.get_real)(&rd, &(rec[8])); reax->tbp[j][k].p_bo6     = val;
      reax->tbp[k][j].p_bo6     = val;
      val = (*rd.get_real)(&rd, &(rec[9])); reax->tbp[j][k].p_ovun1 = val;
      reax->tbp[k][j].p_ovun1 = val;

      /* line 2 */
      c = (*rd.read_record)(&rd, rec);

      val = (*rd.get_real)(&rd, &(rec[0])); reax->tbp[j][k].p_be2     = val;
      reax->tbp[k][j].p_be2     = val;
      val = (*rd.get_real)(&rd, &(rec[1])); reax->tbp[j][k].p_bo3     = val;
      reax->tbp[k][j].p_bo3     = val;
      val = (*rd.get_real)(&rd, &(rec[2])); reax->tbp[j][k].p_bo4     = val;
      reax->tbp[k][j].p_bo4     = val;
      val = (*rd.get_real)(&rd, &(rec[3]));

      val = (*rd.get_real)(&rd, &(rec[4])); reax->tbp[j][k].p_bo1     = val;
      reax->tbp[k][j].p_bo1     = val;
      val = (*rd.get_real)(&rd, &(rec[5])); reax->tbp[j][k].p_bo2     = val;
      reax->tbp[k][j].p_bo2     = val;
      val = (*rd.get_real)(&rd, &(rec[6])); reax->tbp[j][k].ovc       = val;
      reax->tbp[k][j].ovc       = val;

      val = (*rd.get_real)(&rd, &(rec[7]));
    }
  }

  for (i=0; i < reax->num_atom_types; i++)
    for (j=i; j < reax->num_atom_types; j++) {
      reax->tbp[i][j].r_s = 0.5 *
        (reax->sbp[i].r_s + reax->sbp[j].r_s);
      reax->tbp[j][i].r_s = 0.5 *
        (reax->sbp[j].r_s + reax->sbp[i].r_s);

      reax->tbp[i][j].r_p = 0.5 *
        (reax->sbp[i].r_pi + reax->sbp[j].r_pi);
      reax->tbp[j][i].r_p = 0.5 *
        (reax->sbp[j].r_pi + reax->sbp[i].r_pi);

      reax->tbp[i][j].r_pp = 0.5 *
        (reax->sbp[i].r_pi_pi + reax->sbp[j].r_pi_pi);
      reax->tbp[j][i].r_pp = 0.5 *
        (reax->sbp[j].r_pi_pi + reax->sbp[i].r_pi_pi);


      reax->tbp[i][j].p_boc3 =
        sqrt(reax->sbp[i].b_o_132 *
             reax->sbp[j].b_o_132);
      reax->tbp[j][i].p_boc3 =
        sqrt(reax->sbp[j].b_o_132 *
             reax->sbp[i].b_o_132);

      reax->tbp[i][j].p_boc4 =
        sqrt(reax->sbp[i].b_o_131 *
             reax->sbp[j].b_o_131);
      reax->tbp[j][i].p_boc4 =
        sqrt(reax->sbp[j].b_o_131 *
             reax->sbp[i].b_o_131);

      reax->tbp[i][j].p_boc5 =
        sqrt(reax->sbp[i].b_o_133 *
             reax->sbp[j].b_o_133);
      reax->tbp[j][i].p_boc5 =
        sqrt(reax->sbp[j].b_o_133 *
             reax->sbp[i].b_o_133);


      reax->tbp[i][j].D =
        sqrt(reax->sbp[i].epsilon *
             reax->sbp[j].epsilon);

      reax->tbp[j][i].D =
        sqrt(reax->sbp[j].epsilon *
             reax->sbp[i].epsilon);

      reax->tbp[i][j].alpha =
        sqrt(reax->sbp[i].alpha *
             reax->sbp[j].alpha);

      reax->tbp[j][i].alpha =
        sqrt(reax->sbp[j].alpha *
             reax->sbp[i].alpha);

      reax->tbp[i][j].r_vdW =
        2.0 * sqrt(reax->sbp[i].r_vdw * reax->sbp[j].r_vdw);

      reax->tbp[j][i].r_vdW =
        2.0 * sqrt(reax->sbp[j].r_vdw * reax->sbp[i].r_vdw);

      reax->tbp[i][j].gamma_w =
        sqrt(reax->sbp[i].gamma_w *
             reax->sbp[j].gamma_w);

      reax->tbp[j][i].gamma_w =
        sqrt(reax->sbp[j].gamma_w *
             reax->sbp[i].gamma_w);

      reax->tbp[i][j].gamma =
        pow(reax->sbp[i].gamma *
            reax->sbp[j].gamma,-1.5);

      reax->tbp[j][i].gamma =
        pow(reax->sbp[j].gamma *
            reax->sbp[i].gamma,-1.5);

      // additions for additional vdWaals interaction types - inner core

      reax->tbp[i][j].rcore = reax->tbp[j][i].rcore =
        sqrt( reax->sbp[i].rcore2 * reax->sbp[j].rcore2 );

      reax->tbp[i][j].ecore = reax->tbp[j][i].ecore =
        sqrt( reax->sbp[i].ecore2 * reax->sbp[j].ecore2 );

      reax->tbp[i][j].acore = reax->tbp[j][i].acore =
        sqrt( reax->sbp[i].acore2 * reax->sbp[j].acore2 );

      // additions for additional vdWalls interaction types lg correction

      reax->tbp[i][j].lgcij = reax->tbp[j][i].lgcij =
        sqrt( reax->sbp[i].lgcij * reax->sbp[j].lgcij );

      reax->tbp[i][j].lgre = reax->tbp[j][i].lgre = 2.0 *
        sqrt( reax->sbp[i].lgre*reax->sbp[j].lgre );

    }

  c = (*rd.read_record)(&rd, rec);
  l = (*rd.get_int)(&rd, &(rec[0]));

  for (i=0; i < l; i++) {
    c = (*rd.read_record)(&rd, rec);

    j = (*rd.get_int)(&rd, &(rec[0])) - 1;
    k = (*rd.get_int)(&rd, &(rec[1])) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types)        {
      val = (*rd.get_real)(&rd, &(rec[2]));
      if (val > 0.0) {
        reax->tbp[j][k].D = val;
        reax->tbp[k][j].D = val;
      }

      val = (*rd.get_real)(&rd, &(rec[3]));
      if (val > 0.0) {
        reax->tbp[j][k].r_vdW = 2 * val;
        reax->tbp[k][j].r_vdW = 2 * val;
      }

      val = (*rd.get_real)(&rd, &(rec[4]));
      if (val > 0.0) {
        reax->tbp[j][k].alpha = val;
        reax->tbp[k][j].alpha = val;
      }

      val = (*rd.get_real)(&rd, &(rec[5]));
      if (val > 0.0) {
        reax->tbp[j][k].r_s = val;
        reax->tbp[k][j].r_s = val;
      }

      val = (*rd.get_real)(&rd, &(rec[6]));
      if (val > 0.0) {
        reax->tbp[j][k].r_p = val;
        reax->tbp[k][j].r_p = val;
      }

      val = (*rd.get_real)(&rd, &(rec[7]));
      if (val > 0.0) {
        reax->tbp[j][k].r_pp = val;
        reax->tbp[k][j].r_pp = val;
      }

      val = (*rd.get_real)(&rd, &(rec[8]));
      if (val >= 0.0) {
        reax->tbp[j][k].lgcij = val;
        reax->tbp[k][j].lgcij = val;
      }
    }
  }

  for( i = 0; i < reax->num_atom_types; ++i )
    for( j = 0; j < reax->num_atom_types; ++j )
      for( k = 0; k < reax->num_atom_types; ++k )
        reax->thbp[i][j][k].cnt = 0;

  c = (*rd.read_record)(&rd, rec);
  l = (*rd.get_int)(&rd, &(rec[0]));

  for( i = 0; i < l; i++ ) {
    c = (*rd.read_record)(&rd, rec);

    j = (*rd.get_int)(&rd, &(rec[0])) - 1;
    k = (*rd.get_int)(&rd, &(rec[1])) - 1;
    m = (*rd.get_int)(&rd, &(rec[2])) - 1;

    if (j < reax->num_atom_types && k < reax->num_atom_types &&
        m < reax->num_atom_types) {
      cnt = reax->thbp[j][k][m].cnt;
      reax->thbp[j][k][m].cnt++;
      reax->thbp[m][k][j].cnt++;

      val = (*rd.get_real)(&rd, &(rec[3]));
      reax->thbp[j][k][m].prm[cnt].theta_00 = val;
      reax->thbp[m][k][j].prm[cnt].theta_00 = val;

      val = (*rd.get_real)(&rd, &(rec[4]));
      reax->thbp[j][k][m].prm[cnt].p_val1 = val;
      reax->thbp[m][k][j].prm[cnt].p_val1 = val;

      val = (*rd.get_real)(&rd, &(rec[5]));
      reax->thbp[j][k][m].prm[cnt].p_val2 = val;
      reax->thbp[m][k][j].prm[cnt].p_val2 = val;

      val = (*rd.get_real)(&rd, &(rec[6]));
      reax->thbp[j][k][m].prm[cnt].p_coa1 = val;
      reax->thbp[m][k][j].prm[cnt].p_coa1 = val;

      val = (*rd.get_real)(&rd, &(rec[7]));
      reax->thbp[j][k][m].prm[cnt].p_val7 = val;
      reax->thbp[m][k][j].prm[cnt].p_val7 = val;

      val = (*rd.get_real)(&rd, &(rec[8]));
      reax->thbp[j][k][m].prm[cnt].p_pen1 = val;
      reax->thbp[m][k][j].prm[cnt].p_pen1 = val;

      val = (*rd.get_real)(&rd, &(rec[9]));
      reax->thbp[j][k][m].prm[cnt].p_val4 = val;
      reax->thbp[m][k][j].prm[cnt].p_val4 = val;
    }
  }

  /* clear all entries first */
  for( i = 0; i < reax->num_atom_types; ++i )
    for( j = 0; j < reax->num_atom_types; ++j )
      for( k = 0; k < reax->num_atom_types; ++k )
        for( m = 0; m < reax->num_atom_types; ++m ) {
          reax->fbp[i][j][k][m].cnt = 0;
          tor_flag[i][j][k][m] = 0;
        }

  /* next line is number of 4-body params and some comments */
  c = (*rd.read_record)(&rd, rec);
  l = (*rd.get_int)(&rd, &(rec[0]));

  for( i = 0; i < l; i++ ) {
    c = (*rd.read_record)(&rd, rec);

    j = (*rd.get_int)(&rd, &(rec[0])) - 1;
    k = (*rd.get_int)(&rd, &(rec[1])) - 1;
    m = (*rd.get_int)(&rd, &(rec[2])) - 1;
    n = (*rd.get_int)(&rd, &(rec[3])) - 1;

    if (j >= 0 && n >= 0) { // this means the entry is not in compact form
      if (j < reax->num_atom_types && k < reax->num_atom_types &&
          m < reax->num_atom_types && n < reax->num_atom_types) {
        tor_flag[j][k][m][n] = 1;
        tor_flag[n][m][k][j] = 1;

        reax->fbp[j][k][m][n].cnt = 1;
        reax->fbp[n][m][k][j].cnt = 1;

        val = (*rd.get_real)(&rd, &(rec[4]));
        reax->fbp[j][k][m][n].prm[0].V1 = val;
        reax->fbp[n][m][k][j].prm[0].V1 = val;

        val = (*rd.get_real)(&rd, &(rec[5]));
        reax->fbp[j][k][m][n].prm[0].V2 = val;
        reax->fbp[n][m][k][j].prm[0].V2 = val;

        val = (*rd.get_real)(&rd, &(rec[6]));
        reax->fbp[j][k][m][n].prm[0].V3 = val;
        reax->fbp[n][m][k][j].prm[0].V3 = val;

        val = (*rd.get_real)(&rd, &(rec[7]));
        reax->fbp[j][k][m][n].prm[0].p_tor1 = val;
        reax->fbp[n][m][k][j].prm[0].p_tor1 = val;

        val = (*rd.get_real)(&rd, &(rec[8]));
        reax->fbp[j][k][m][n].prm[0].p_cot1 = val;
        reax->fbp[n][m][k][j].prm[0].p_cot1 = val;
      }
    }
    else { /* This means the entry is of the form 0-X-Y-0 */
      if( k < reax->num_atom_types && m < reax->num_atom_types )
        for( p = 0; p < reax->num_atom_types; p++ )
          for( o = 0; o < reax->num_atom_types; o++ ) {
            reax->fbp[p][k][m][o].cnt = 1;
            reax->fbp[o][m][k][p].cnt = 1;

            if (tor_flag[p][k][m][o] == 0) {
              reax->fbp[p][k][m][o].prm[0].V1 = (*rd.get_real)(&rd, &(rec[4]));
              reax->fbp[p][k][m][o].prm[0].V2 = (*rd.get_real)(&rd, &(rec[5]));
              reax->fbp[p][k][m][o].prm[0].V3 = (*rd.get_real)(&rd, &(rec[6]));
              reax->fbp[p][k][m][o].prm[0].p_tor1 = (*rd.get_real)(&rd, &(rec[7]));
              reax->fbp[p][k][m][o].prm[0].p_cot1 = (*rd.get_real)(&rd, &(rec[8]));
            }

            if (tor_flag[o][m][k][p] == 0) {
              reax->fbp[o][m][k][p].prm[0].V1 = (*rd.get_real)(&rd, &(rec[4]));
              reax->fbp[o][m][k][p].prm[0].V2 = (*rd.get_real)(&rd, &(rec[5]));
              reax->fbp[o][m][k][p].prm[0].V3 = (*rd.get_real)(&rd, &(rec[6]));
              reax->fbp[o][m][k][p].prm[0].p_tor1 = (*rd.get_real)(&rd, &(rec[7]));
              reax->fbp[o][m][k][p].prm[0].p_cot1 = (*rd.get_real)(&rd, &(rec[8]));
            }
          }
    }
  }



  /* next line is number of hydrogen bond params and some comments */
  c = (*rd.read_record)(&rd, rec);
  l = (*rd.get_int)(&rd, &(rec[0]));

  for( i = 0; i < l; i++ ) {
    c = (*rd.read_record)(&rd, rec);

    j = (*rd.get_int)(&rd, &(rec[0])) - 1;
    k = (*rd.get_int)(&rd, &(rec[1])) - 1;
    m = (*rd.get_int)(&rd, &(rec[2])) - 1;


    if( j < reax->num_atom_types && m < reax->num_atom_types ) {
      val = (*rd.get_real)(&rd, &(rec[3]));
      reax->hbp[j][k][m].r0_hb = val;

      val = (*rd.get_real)(&rd, &(rec[4]));
      reax->hbp[j][k][m].p_hb1 = val;

      val = (*rd.get_real)(&rd, &(rec[5]));
      reax->hbp[j][k][m].p_hb2 = val;

      val = (*rd.get_real)(&rd, &(rec[6]));
      reax->hbp[j][k][m].p_hb3 = val;
    }
  }

cleanup_tor_flag:
  /* deallocate tor_flag */
  for( i = 0; i < reax->num_atom_types; i++ ) {
    for( j = 0; j < reax->num_atom_types; j++ ) {
      for( k = 0; k < reax->num_atom_types; k++ ) {
        free( tor_flag[i][j][k] );
      }
      free( tor_flag[i][j] );
    }
    free( tor_flag[i] );
  }
  free( tor_flag );

cleanup_reader:
  /* deallocate helper storage */
  free(rec);
  (*rd.destroy)(&rd);

  /* close file */
  fclose(fp);

end:
  if (error_code) {
      MPI_Abort( comm, error_code );
  }
  else {
    return SUCCESS;
  }
}
