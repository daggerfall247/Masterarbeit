
/*******************************************************************************
 *
 * File config_IO.c
 *
 * Copyright (C) 2013 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines for the reading and writing of the configurations
 *
 * Externally accessible functions:
 *
 * void writeConfig(char *outfile)
 *      Writes the current configuration of gauge fields into the file whose
 *      name is handed to the function via the string `outfile'.
 *
 * void readConfig(char *infile)
 *      Reads the gauge field from the file whose name is handed to the
 *      function via the string infile and stores it into the standard gauge
 *      field.
 *
 *******************************************************************************/

#define CONFIG_IO_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "misc.h"
#include "headers.h"
#include "modules.h"

#define DEBUG_IO 0

#if (MASTER_FIELD == 0)
void writeConfig(char *outfile)
{
    FILE *fout = NULL;
    uint64_t ii, in, iend;
    uint64_t icheck, icheck2;
    stdint_t lswrite[DIM], info[2];
    double plaq;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("write_config -- in");
    logging("\nWriting configuration %s:\n", outfile);
    t1 = getTime();

#if (DIM == 2)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
#elif (DIM == 3)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
    lswrite[2] = (stdint_t)(LENGS2);
#elif (DIM == 4)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
    lswrite[2] = (stdint_t)(LENGS2);
    lswrite[3] = (stdint_t)(LENGS3);
#endif

    iend = endianness();
    plaq = plaquette();
    icheck = 0;

    fout = fopen(outfile, "wb");
    error(fout == NULL, "write_config [config_IO.c]",
          "Unable to create output file %s!", outfile);

    info[0] = DIM;
    info[1] = SUN;

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lswrite);
        bswap_double(1, &plaq);
    }

    icheck = fwrite(info, sizeof(stdint_t), 2, fout);
    icheck += fwrite(lswrite, sizeof(stdint_t), DIM, fout);
    icheck += fwrite(&plaq, sizeof(double), 1, fout);
    error(icheck != DIM + 3, "write_config [config_IO.c]", "Write error!");

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "write_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;
    for (in = 0; in < VOL; in++)
    {
        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_sun_dble_array(zw, *pu[in][ii]);

        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }
        icheck += fwrite(buff, sizeof(double), SUNVOL * DIM, fout);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "write_config [config_IO.c]", "Write error!");

    fclose(fout);

    free(buff);

    checkpoint("write_config -- out");
    t2 = getTime();
    logging("\nConfiguration exported (took %.3e sec)!\n", t2 - t1);
}
#elif (MASTER_FIELD == 1)

void writeHeaderToConfig(char *outfile)
{
    FILE *fout = NULL;
    stdint_t lswrite[DIM], info[2];
    double t1, t2;
    int iend, icheck;

    checkpoint("write_header_config -- in");
    logging("\nWriting configuration header %s:\n", outfile);
    t1 = getTime();

#if (DIM == 2)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
#elif (DIM == 3)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
    lswrite[2] = (stdint_t)(LENGS2_MF);
#elif (DIM == 4)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
    lswrite[2] = (stdint_t)(LENGS2_MF);
    lswrite[3] = (stdint_t)(LENGS3_MF);
#endif
    iend = endianness();
    icheck = 0;

    fout = fopen(outfile, "wb");
    error(fout == NULL, "write_config [config_IO.c]", "Unable to create output file %s!", outfile);

    info[0] = DIM;
    info[1] = SUN;

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lswrite);
    }

    icheck = fwrite(info, sizeof(stdint_t), 2, fout);
    icheck += fwrite(lswrite, sizeof(stdint_t), DIM, fout);

    error(icheck != DIM + 2, "write_config [config_IO.c]", "Write error!");
    checkpoint("write_header_config -- out");
    t2 = getTime();
    logging("\nConfig Header exported (took %.3e sec)!\n", t2 - t1);
    fclose(fout);
}

void writePlaquetteToConfig(char *outfile, double plaq)
{
    FILE *fout = NULL;
    int iend;
    long int icheck;
    double t1, t2;
    long int offset_header;

    checkpoint("write_plaquette_to_config -- in");
    logging("\nWriting plaquette to configuration %s:\n", outfile);
    t1 = getTime();

    iend = endianness();
    icheck = 0;

    fout = fopen(outfile, "rb+");
    error(fout == NULL, "write_plaquette_to_config [config_IO.c]", "Unable to open output file %s!", outfile);

    if (iend == BIG_ENDIAN)
    {
        bswap_double(1, &plaq);
    }

    offset_header = 8UL + DIM * 4UL;
    fseek(fout, offset_header, SEEK_SET);
    icheck += fwrite(&plaq, sizeof(double), 1, fout);

    error(icheck != 1, "write_plaquette_to_config [config_IO.c]", "Write error!");
    checkpoint("write_plaquette_to_config -- out");
    t2 = getTime();
    logging("\nPlaquette exported (took %.3e sec)\n", t2 - t1);
    fclose(fout);
}

void writeConfig(char *outfile)
{
    FILE *fout = NULL;
    uint64_t ii, in, kk, iend, buffer_length;
    uint64_t icheck, icheck2;
    uint64_t pos;
    long int offset_header;
    double *zw, *buff = NULL;
    double t1, t2;

#if (DIM == 2)
    buffer_length = LENGS1;
#elif (DIM == 3)
    buffer_length = LENGS2;
#elif (DIM == 4)
    buffer_length = LENGS3;
#endif

    checkpoint("write_config -- in");
    logging("\nWriting configuration %s:\n", outfile);
    t1 = getTime();

    iend = endianness();
    icheck = 0;

    fout = fopen(outfile, "rb+");
    error(fout == NULL, "write_config [config_IO.c]", "Unable to create output file %s!", outfile);

    buff = malloc(buffer_length * SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "write_config [config_IO.c]", "Unable to allocate buffers!");

    offset_header = 16UL + 4UL * (unsigned long)(DIM);

    for (in = 0; in < VOL; in += buffer_length)
    {
        for (kk = 0, zw = buff; kk < buffer_length; kk++)
        {
            for (ii = 0; ii < DIM; ii++, zw += SUNVOL)
            {
                mk_sun_dble_array(zw, *pu[i[in + kk]][ii]);
            }
        }
        if (iend == BIG_ENDIAN)
        {
            bswap_double(buffer_length * SUNVOL * DIM, buff);
        }

        pos = index_mf[i[in]];
        pos *= (uint64_t)(sizeof(double) * SUNVOL * DIM);
        fseek(fout, offset_header + pos, SEEK_SET);
        icheck += fwrite(buff, sizeof(double), buffer_length * SUNVOL * DIM, fout);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "write_config [config_IO.c]", "Write error!");

    fclose(fout);

    free(buff);

    checkpoint("write_config -- out");
    t2 = getTime();
    logging("\nConfiguration exported (took %.3e sec)!\n", t2 - t1);
}

void prepareConfig(char *outfile)
{
    uint64_t n;
    writeHeaderToConfig(outfile);
    writePlaquetteToConfig(outfile, 0.0);
    for (n = 0; n < VOL_SL; n++)
    {
        updateArrayOfIndexMF(sl[n]);
        writeConfig(outfile);
    }
}
#endif

#if (MASTER_FIELD == 0)
/*
void readConfig_MarcStyle(const char *filename)
{
    char string1[1000];

    FILE *fd;

    if ((fd = fopen(filename, "r")) == NULL)
    {
        fprintf(stderr, "Error: void read_gauge_field(...!\n");
        fprintf(stderr, "  Filename does not exist.\n");
        exit(EXIT_FAILURE);
    }

    while (1)
    {
        if (fgets(string1, 995, fd) == NULL)
        {
            fprintf(stderr, "Error: void read_gauge_field(...!\n");
            fprintf(stderr, "  Wrong file format.\n");
            exit(EXIT_FAILURE);
        }

        if (strcmp(string1, "BEGIN\n") == 0)
            break;
    }

    uint64_t in, mu;

    for (in = 0; in < VOL; in++)
    {
        for (mu = 0; mu < 4; mu++)
        {
            double h[4];

            if (fread(h, sizeof(double), 4, fd) != 4)
            {
                fprintf(stderr, "Error: void read_gauge_field(...!\n");
                fprintf(stderr, "  File too short!\n");
                exit(EXIT_FAILURE);
            }

            if (fabs(h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3] - 1.0) > 0.0000000001)
            {
                fprintf(stderr, "Error: void read_gauge_field(...!\n");
                fprintf(stderr, "  in = %2d, mu = %2d   -->   h^2 = %.6lf.\n", in, mu, h[0] * h[0] + h[1] * h[1] + h[2] * h[2] + h[3] * h[3]);
                exit(EXIT_FAILURE);
            }

            mk_dble_array_sun(h, *pu[in][mu]);
        }
    }

    double h[4];

    if (fread(h, sizeof(double), 4, fd) != 0)
    {
        fprintf(stderr, "Error: void read_gauge_field(...!\n");
        fprintf(stderr, "  File too long!\n");
        exit(EXIT_FAILURE);
    }

    fclose(fd);
}
*/

void readConfigForThermalisation(char *infile)
{
    FILE *fin = NULL;
    uint64_t ileng[DIM], ileng_half[DIM], VOL_small;
    uint64_t volumeFactor[DIM], volumeFactor_small[DIM];
    int N = pow(2, DIM);
    int mirror[N][DIM];
    uint64_t ii, kk, jj, in, iend, index;
    uint64_t icheck, icheck2;
    uint64_t coord[DIM];
    stdint_t lscheck[DIM], info[2];
    sun_mat U[DIM], U_dag[DIM];
    double plaq, plaq0, eps;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("read_config_thermalisation -- in");
    logging("\nReading configuration for thermalisation %s:\n", infile);
    t1 = getTime();

    VOL_small = VOL / pow(2, DIM);
#if (DIM == 2)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng_half[0] = LENGT / 2;
    ileng_half[1] = LENGS1 / 2;
    volumeFactor[0] = LENGS1;
    volumeFactor[1] = 1;
    volumeFactor_small[0] = LENGS1 / 2;
    volumeFactor_small[1] = 1;
#elif (DIM == 3)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
    ileng_half[0] = LENGT / 2;
    ileng_half[1] = LENGS1 / 2;
    ileng_half[2] = LENGS2 / 2;
    volumeFactor[0] = LENGS1 * LENGS2;
    volumeFactor[1] = LENGS2;
    volumeFactor[2] = 1;
    volumeFactor_small[0] = LENGS1 * LENGS2 / 4;
    volumeFactor_small[1] = LENGS2 / 2;
    volumeFactor_small[2] = 1;
#elif (DIM == 4)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
    ileng[3] = LENGS3;
    ileng_half[0] = LENGT / 2;
    ileng_half[1] = LENGS1 / 2;
    ileng_half[2] = LENGS2 / 2;
    ileng_half[3] = LENGS3 / 2;
    volumeFactor[0] = LENGS1 * LENGS2 * LENGS3;
    volumeFactor[1] = LENGS2 * LENGS3;
    volumeFactor[2] = LENGS3;
    volumeFactor[3] = 1;
    volumeFactor_small[0] = LENGS1 * LENGS2 * LENGS3 / 8;
    volumeFactor_small[1] = LENGS2 * LENGS3 / 4;
    volumeFactor_small[2] = LENGS3 / 2;
    volumeFactor_small[3] = 1;
#endif

    error(pu[0][0] == NULL, "read_config_thermalisation [config_IO.c]", "Fields are not allocated!");

    iend = endianness();
    icheck = 0;

    fin = fopen(infile, "rb");
    error(fin == NULL, "read_config_thermalisation [config_IO.c]",
          "Unable to read input file %s!", infile);
    icheck = fread(info, sizeof(stdint_t), 2, fin);
    icheck += fread(lscheck, sizeof(stdint_t), DIM, fin);
    icheck += fread(&plaq0, sizeof(double), 1, fin);
    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lscheck);
        bswap_double(1, &plaq0);
    }
    error(icheck != DIM + 3, "read_config_thermalisation [config_IO.c]", "Read error!");
    error((info[0] != DIM) || (info[1] != SUN), "read_config_thermalisation [config_IO.c]", "Incompatible parameters!");
    for (ii = 0; ii < DIM; ii++)
    {
        error(lscheck[ii] != ileng_half[ii], "read_config_thermalisation [config_IO.c]", "Incompatible lattice size!");
    }
    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "read_config_thermalisation [config_IO.c]", "Unable to allocate buffers!");

    for (ii = 0; ii < N; ii++)
    {
        //logging("mirror %d :", ii);
        for (kk = 0; kk < DIM; kk++)
        {
            mirror[ii][kk] = (ii / (int)pow(2, kk)) % 2;
            //logging(" %d", mirror[ii][kk]);
        }
        //logging("\n\n");
    }

    icheck = 0;

    for (in = 0; in < VOL_small; in++) // for each original lattice site
    {
        icheck += fread(buff, sizeof(double), SUNVOL * DIM, fin);
        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL) // temporarily store all links
        {
            mk_dble_array_sun(zw, U[ii]);
            sun_dag(U_dag[ii], U[ii]);
        }

        index = in;
        for (ii = 0; ii < DIM; ii++) // calculate coordinates of current index
        {
            coord[ii] = index / volumeFactor_small[ii];
            index = index % volumeFactor_small[ii];
        }

        //logging("coordinates of index %d: ", in);
        for (ii = 0; ii < DIM; ii++)
        {
            //logging("%d ", coord[ii]);
        }
        //logging("\n");

        for (ii = 0; ii < DIM; ii++) // for each link
        {
            for (kk = 0; kk < N; kk++) // for each sublattice..
            {
                index = 0;
                for (jj = 0; jj < DIM; jj++) // .. calculate the index of the site, and ..
                {
                    if (mirror[kk][jj] && jj == ii)
                    {
                        index += volumeFactor[jj] * (ileng[jj] - coord[jj] - 1);
                    }
                    else if (mirror[kk][jj] && jj != ii)
                    {
                        if (coord[jj] == 0)
                        {
                            index += volumeFactor[jj] * (coord[jj] + ileng_half[jj]);
                        }
                        else if (coord[jj] != 0)
                        {
                            index += volumeFactor[jj] * (ileng[jj] - coord[jj]);
                        }
                    }
                    else
                    {
                        index += volumeFactor[jj] * coord[jj];
                    }
                }
                //logging("index for %d-dir link in mirrored sublattice %d: %d\n", ii, kk, index);

                // .. place the links into container.
                if (mirror[kk][ii])
                {
                    *pu[index][ii] = U_dag[ii];
                }
                else
                {
                    *pu[index][ii] = U[ii];
                }
            }
        }
    }

    icheck2 = DIM * SUNVOL * VOL / N;

    error(icheck != icheck2, "read_config_thermalisation [config_IO.c]", "Read error!");

    fclose(fin);

    plaq = plaquette();
    eps = sqrt((double)(NPLAQ * VOL)) * DBL_EPSILON;
    printf("read plaq %.6e, calc plaq %.6e\n", plaq0, plaq);
    error(fabs(plaq - plaq0) > eps, "read_config_thermalisation [config_IO.c]", "Plaquette test failed!");

    free(buff);

    checkpoint("read_config_thermalisation -- out");
    t2 = getTime();
    logging("\nConfiguration imported for thermalisation(took %.3e sec)!\n", t2 - t1);
}

void readConfig(char *infile)
{
    FILE *fin = NULL;
    uint64_t ii, in, iend;
    uint64_t icheck, icheck2;
    uint64_t ileng[DIM];
    stdint_t lscheck[DIM], info[2];
    double plaq, plaq0, eps;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("read_config -- in");
    logging("\nReading configuration %s:\n", infile);
    t1 = getTime();

#if (DIM == 2)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
#elif (DIM == 3)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
#elif (DIM == 4)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
    ileng[3] = LENGS3;
#endif

    error(pu[0][0] == NULL, "read_config [config_IO.c]", "Fields are not allocated!");

    iend = endianness();
    icheck = 0;

    fin = fopen(infile, "rb");
    error(fin == NULL, "read_config [config_IO.c]",
          "Unable to read input file %s!", infile);
    icheck = fread(info, sizeof(stdint_t), 2, fin);
    icheck += fread(lscheck, sizeof(stdint_t), DIM, fin);
    icheck += fread(&plaq0, sizeof(double), 1, fin);
    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lscheck);
        bswap_double(1, &plaq0);
    }
    error(icheck != DIM + 3, "read_config [config_IO.c]", "Read error!");
    error((info[0] != DIM) || (info[1] != SUN), "read_config [config_IO.c]", "Incompatible parameters!");
    for (ii = 0; ii < DIM; ii++)
    {
        error(lscheck[ii] != ileng[ii], "read_config [config_IO.c]", "Incompatible lattice size!");
    }

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "read_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;
    for (in = 0; in < VOL; in++)
    {
        icheck += fread(buff, sizeof(double), SUNVOL * DIM, fin);
        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_dble_array_sun(zw, *pu[in][ii]);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "read_config [config_IO.c]", "Read error!");

    fclose(fin);

    plaq = plaquette();
    eps = sqrt((double)(NPLAQ * VOL)) * DBL_EPSILON;
    printf("read plaq %.6e, calc plaq %.6e\n", plaq0, plaq);
    error(fabs(plaq - plaq0) > eps, "read_config [config_IO.c]", "Plaquette test failed!");

    free(buff);

    checkpoint("read_config -- out");
    t2 = getTime();
    logging("\nConfiguration imported (took %.3e sec)!\n", t2 - t1);
}
#elif (MASTER_FIELD == 1)
void readConfig(char *infile)
{
    FILE *fin = NULL;
    uint64_t ii, in, iend;
    uint64_t icheck, icheck2;
    uint64_t pos;
    long int offset_header;
    uint64_t ileng_mf[DIM];
    stdint_t info[2], lscheck_mf[DIM];
    double plaq0;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("read_config -- in");
    logging("\nReading configuration %s:\n", infile);
    t1 = getTime();

/*added master-field extents*/
#if (DIM == 2)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
#elif (DIM == 3)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
    ileng_mf[2] = LENGS2_MF;
#elif (DIM == 4)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
    ileng_mf[2] = LENGS2_MF;
    ileng_mf[3] = LENGS3_MF;
#endif

    error(pu[0][0] == NULL, "read_config [config_IO.c]", "Fields are not allocated!");

    iend = endianness();
    icheck = 0;

    fin = fopen(infile, "rb");
    error(fin == NULL, "read_config [config_IO.c]", "Unable to read input file %s!", infile);
    icheck = fread(info, sizeof(stdint_t), 2, fin);
    icheck += fread(lscheck_mf, sizeof(stdint_t), DIM, fin); /*added master_field extent values*/
    icheck += fread(&plaq0, sizeof(double), 1, fin);

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lscheck_mf); /*mf*/
        bswap_double(1, &plaq0);
    }

    error(icheck != DIM + 3, "read_config [config_IO.c]", "Read error! (Master-field)");
    error((info[0] != DIM) || (info[1] != SUN), "read_config [config_IO.c]", "Incompatible parameters!");
    for (ii = 0; ii < DIM; ii++)
    {
        error(lscheck_mf[ii] != ileng_mf[ii], "read_config [config_IO.c]", "Incompatible master_field size!"); /*added master-field size check*/
    }

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "read_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;

    offset_header = 16UL + 4UL * (unsigned long)(DIM);

    for (in = 0; in < VOL2; in++)
    {
        pos = index_mf[in];
        pos *= (uint64_t)(SUNVOL * DIM * sizeof(double));
        fseek(fin, offset_header + pos, SEEK_SET); // TODO
        icheck += fread(buff, sizeof(double), SUNVOL * DIM, fin);

        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
        {
            mk_dble_array_sun(zw, *pu[in][ii]);
        }
    }

    icheck2 = DIM * SUNVOL * VOL2;
    error(icheck != icheck2, "read_config [config_IO.c]", "Read error!");

    fclose(fin);

    // plaq = plaquette();
    // eps = sqrt((double)(NPLAQ * VOL)) * DBL_EPSILON;
    // error(fabs(plaq - plaq0) > eps, "read_config [config_IO.c]", "Plaquette test failed!");

    free(buff);

    checkpoint("read_config -- out");
    t2 = getTime();
    logging("\nConfiguration imported (took %.3e sec)!\n", t2 - t1);
}

void readConfigForThermalisation(char *infile)
{
}

#endif