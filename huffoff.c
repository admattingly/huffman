/* huffoff.c */

/*
 * Use the entropy encoding feature of the CMPSC instruction to expand a file
 * compressed using entropy encoding only - no Lempel-Ziv dictionary of multi-character
 * sequences.
 *
 * To compile under z/OS UNIX:
 *
 * $ xlc -qasm -qasmlib=sys1.maclib -qin=all:nostp -ohuffoff huffoff.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FALSE               0
#define TRUE                1

/* compression modes */
#define MODE_COPY           0
#define MODE_ENTROPY        1

/* tuning parameters */
#define BUFFER_SIZE         65536   // size of buffer for reading/writing files
#define MINIMUM_SIZE        1024    // input file size below which compression is unnecessary
#define MAXBITS             16      // maximum codeword length permitted by CMPSC

/* CMPSC option bits (GPR0) */
#define SC_EXPAND       0x100
#define SC_FORMAT1      0x200
#define SC_SYMTRAN      0x10000
#define SC_ZEROPAD      0x20000
#define SC_ORDER        0x40000
#define SC_ENTROPY      0x80000

/* compressed file header */
#define HUFF_EYECATCHER     "HUFF"
#define HUFF_ECLENGTH       4
#define HUFF_VERSION        1
#pragma pack(packed)
typedef struct _HUFF_HEADER {
    char            eyecatcher[HUFF_ECLENGTH];
    unsigned short  version;
    unsigned short  mode;
} HUFF_HEADER, *PHUFF_HEADER;

typedef struct _HUFF_EE_EXTEND {
    short           ed[MAXBITS];
    unsigned char   alphabet[256];
} HUFF_EE_EXTEND, *PHUFF_EE_EXTEND;
#pragma pack(reset)

/* package structure for Package Merge algorithm */
#pragma pack(packed)
typedef struct _PACKAGE {
    int     freq;
    int     symbol;
} PACKAGE, *PPACKAGE;
#pragma pack(reset)

/* private function prototypes */
static int buildDictionary(unsigned char *buffer, unsigned char **pedict, unsigned int *pedictlen);
static int copy(FILE *fin, FILE *fout);
static int expand(FILE *fin, FILE *fout);
static void free4k(void *addr, unsigned int size);
static void *malloc4k(unsigned int size);
static int openFiles(char *infile, char *outfile, FILE **pfin, FILE **pfout, int *pmode);
static int parseCommandLine(int argc, char **argv, char **pinfile, char **poutfile);
static void printParms(int              callno,
                       char             *when,
                       unsigned long    gpr0,
                       unsigned long    gpr1,
                       int              cc,
                       unsigned char    *fopaddr,
                       unsigned long    foplen,
                       unsigned char    *sopaddr,
                       unsigned long    soplen);

int main(int argc, char **argv)
{
    int             compressMode, rc;
    char            *infile, *outfile;
    FILE            *fin, *fout;

    rc = 8;
    if (parseCommandLine(argc, argv, &infile, &outfile)) {
        if (openFiles(infile, outfile, &fin, &fout, &compressMode)) {
            if (compressMode == MODE_ENTROPY) {
                if (expand(fin, fout)) {
                    rc = 0;
                }
            }
            else {
                if (copy(fin, fout)) {
                    rc = 0;
                }
            }
            fclose(fin);
            fclose(fout);
        }
    }
    return rc;
}

int buildDictionary(unsigned char *buffer, unsigned char **pedict, unsigned int *pedictlen)
{
    int             i, result;
    const int       n = 256;            // dictionary has 256 entries
    unsigned char   *edict;             // compression dictionary
    unsigned int    edictlen;
    unsigned char   entry[8];           // dictionary entry
    PHUFF_EE_EXTEND phee;               // ED and descending-frequency-sequence alphabet

    result = TRUE;
    edict  = NULL;

    /* allocate space for dictionary on 4KB page boundary, as required by CMPSC */
    edictlen = n * 8 + MAXBITS * sizeof(short);  // room for classic dictionary and ED
    edict = malloc4k(edictlen);
    if (edict == NULL) {
        result = FALSE;
        fprintf(stderr, "ERROR: malloc4k() of expansion dictionary failed.");
    }

    /* classic dictionary comprises 256-character alphabet in descending frequency sequence */
    phee = (PHUFF_EE_EXTEND)(void *)buffer;
    memset(entry, 0x00, 8);
    entry[0] = 1;   // PSL = 0 and CSL = 1 for every entry (i.e. unpreceded entry - see POP Appendix D)
    for (i = 0; i < 256; i++) {
        entry[1] = phee->alphabet[i];
        memcpy(edict + 8 * i, entry, 8);
    }

    /* copy entity descriptor after classic dictionary */
    memcpy(edict + 8 * n, phee->ed, MAXBITS * sizeof(short));

    *pedict    = edict;
    *pedictlen = edictlen;

    return result;
}

int copy(FILE *fin, FILE *fout)
{
    int             inbytes, outbytes, result;
    unsigned char   *buffer;
    PHUFF_HEADER    phh;

    result = TRUE;

    buffer = malloc(BUFFER_SIZE);
    if (buffer == NULL) {
        result = FALSE;
        perror("ERROR: malloc() of I/O buffer");
    }
    else {
        /* write the compressed file header */
        phh = (PHUFF_HEADER)(void *)buffer;
        memcpy(phh->eyecatcher, HUFF_EYECATCHER, HUFF_ECLENGTH);
        phh->version = HUFF_VERSION;
        phh->mode    = MODE_COPY;
        outbytes = fwrite(buffer, 1, sizeof(HUFF_HEADER), fout);
        if (outbytes != sizeof(HUFF_HEADER)) {
            result = FALSE;
            perror("ERROR: fwrite() of compressed file header");
        }
        else {
            /* copy the input file to the output file */
            while (result && !feof(fin)) {
                inbytes = fread(buffer, 1, BUFFER_SIZE, fin);
                if (inbytes > 0) {
                    outbytes = fwrite(buffer, 1, inbytes, fout);
                    if (outbytes != inbytes) {
                        result = FALSE;
                        perror("ERROR: fwrite() writing what was read");
                    }
                }
            }
        }
    }
    return result;
}

int expand(FILE *fin, FILE *fout)
{
    int             inbytes, outbytes, intotal, outtotal, ncalls, result;
    unsigned char   *inbuffer, *outbuffer, *edict;
    unsigned int    edictlen;
    int             cbn, cc, offset, foplen, soplen;
    unsigned char   *fopaddr, *sopaddr;
    long            gpr0, gpr1;

    result    = TRUE;
    inbuffer  = NULL;
    outbuffer = NULL;
    intotal   = sizeof(HUFF_HEADER);    // header must have already been read
    outtotal  = 0;

    inbuffer = malloc(BUFFER_SIZE);
    if (inbuffer == NULL) {
        result = FALSE;
        perror("ERROR: malloc() of input buffer");
    }
    else {
        outbuffer = malloc(BUFFER_SIZE);
        if (outbuffer == NULL) {
            result = FALSE;
            perror("ERROR: malloc() of output buffer");
        }
    }

    if (result) {
        /* read the entropy header */
        inbytes = fread(inbuffer, 1, sizeof(HUFF_EE_EXTEND), fin);
        intotal += inbytes;
        if (inbytes != sizeof(HUFF_EE_EXTEND)) {
            result = FALSE;
            perror("ERROR: fread() of entropy header extension\n");
        }
        else {
            result = buildDictionary(inbuffer, &edict, &edictlen);
        }
    }

    /* perform expansion operation */
    if (result) {
        cbn     = 0;    // compressed input starts on a byte boundary
        ncalls  = 0;
        fopaddr = outbuffer;
        foplen  = BUFFER_SIZE;
        while (result && !feof(fin)) {
            inbytes = fread(inbuffer, 1, BUFFER_SIZE, fin);
            if (inbytes > 0) {
                intotal += inbytes;
                sopaddr = inbuffer;
                soplen  = inbytes;
                cc      = 3;
                offset  = 2048;     // ED starts at 8 * 256 = 2048
                while (cc == 3) {
                    gpr0 = SC_ENTROPY | SC_EXPAND;
                    gpr1 = ((unsigned long)edict) + (((offset % 65536) / 128) << 3) + cbn;
                    ncalls++;
//                        printParms(ncalls, "Before", gpr0, gpr1, cc, fopaddr, foplen, sopaddr, soplen);
                    __asm(" CMPSC %[out],%[in]\n"
                          " IPM   %[cc]       \n"
                          " SRL   %[cc],28      "
                          :       "+NR:r0"(gpr0),
                                  "+NR:r1"(gpr1),
                            [cc]      "=r"(cc),
                            [in]   "+RP:e"(sopaddr),
                                      "+r"(soplen),
                            [out]  "+RP:e"(fopaddr),
                                      "+r"(foplen)
                          :
                          : );
//                        printParms(ncalls, "After ", gpr0, gpr1, cc, fopaddr, foplen, sopaddr, soplen);
                    /* recover CBN from GRP1 */
                    cbn = gpr1 & 0x07;
                    switch (cc) {
                    case 0: /* input exhausted - read more input (if any) */
                        if (soplen > 0) {   /* incomplete symbol at end of input buffer */
                            printf("WARNING: cc=0 and second op len > 0\n");
                            memmove(inbuffer, sopaddr, soplen);
                        }
                        /* fall out of inner CMPSC loop when cc = 0, so outer loop can read more input */
                        break;
                    case 1: /* output buffer exhausted */
                        outbytes = fwrite(outbuffer, 1, fopaddr - outbuffer, fout);
                        outtotal += outbytes;
                        if (outbytes != fopaddr - outbuffer) {
                            result = FALSE;
                            perror("ERROR: fwrite() to output file");
                        }
                        if (foplen > 0) {
                            printf("WARNING: cc=1 and first op len > 0\n");
                            memmove(outbuffer, fopaddr, foplen);   /* in case of incomplete symbol */
                        }
                        fopaddr = outbuffer;
                        foplen = BUFFER_SIZE;
                        cc = 3;
                        break;
                    default:
                        break;
                    }
                }
            }
        }
        if (result) {
            if (cc == 0 && fopaddr > outbuffer) {  /* some output left to write */
                outbytes = fwrite(outbuffer, 1, fopaddr - outbuffer, fout);
                outtotal += outbytes;
                if (outbytes != fopaddr - outbuffer) {
                    result = FALSE;
                    perror("ERROR: fwrite() to output file");
                }
            }
        }
    }

    /* produce a report on space savings */
    if (result) {
        printf("Input file size:  %10d bytes\n", intotal);
        printf("Output file size: %10d bytes\n", outtotal);
    }

    /* tidy up */
    if (inbuffer != NULL) {
        free(inbuffer);
    }
    if (outbuffer != NULL) {
        free(outbuffer);
    }
    if (edict != NULL) {
        free4k(edict, edictlen);
    }

    return result;
}

void free4k(void *addr, unsigned int size)
{
    int     rc;

    __asm (" SYSSTATE ARCHLVL=2                      \n"
           " STORAGE RELEASE,ADDR=(1),LENGTH=(%[len])  "
           :        "+NR:r1"(addr),
                   "=NR:r15"(rc)
           : [len]       "r"(size)
           : "r0", "r1", "r14", "r15");
    if (rc != 0) {
        fprintf(stderr, "ERROR: STORAGE RELEASE failed, rc=%d\n", rc);
    }
}

void *malloc4k(unsigned int size)
{
    int     rc;
    void    *addr;

    addr = NULL;
    __asm(" SYSSTATE ARCHLVL=2                       \n"
          " STORAGE OBTAIN,LENGTH=(%[len]),BNDRY=PAGE  "
         :        "+NR:r1"(addr),
                 "=NR:r15"(rc)
         : [len]       "r"(size)
         : "r0", "r1", "r14", "r15");
    if (rc != 0) {
        fprintf(stderr, "ERROR: STORAGE OBTAIN failed, rc=%d\n", rc);
    }
    return addr;
}

int openFiles(char *infile, char *outfile, FILE **pfin, FILE **pfout, int *pmode)
{
    int             result, inbytes;
    HUFF_HEADER     hdr;

    result = TRUE;
    *pfin  = NULL;
    *pfout = NULL;

    *pfin = fopen(infile, "rb");
    if (*pfin == NULL) {
        result = FALSE;
        perror("ERROR: fopen() of input file");
    }
    else {
        /* read and validate the input file header */
        inbytes = fread(&hdr, 1, sizeof(HUFF_HEADER), *pfin);
        if (inbytes != sizeof(HUFF_HEADER)) {
            result = FALSE;
            perror("ERROR: fread() of input file header");
        }
        else {
            if (memcmp(hdr.eyecatcher, HUFF_EYECATCHER, HUFF_ECLENGTH) != 0 ||
                hdr.version != HUFF_VERSION) {
                result = FALSE;
                fprintf(stderr, "ERROR: Input file contains invalid header\n");
            }
        }
    }

    /* find the compression mode used by the input file */
    if (result) {
        switch (hdr.mode) {
        case MODE_COPY:
        case MODE_ENTROPY:
            *pmode = hdr.mode;
            break;
        default:
            result = FALSE;
            fprintf(stderr, "ERROR: Input file specifies unsupported compression mode, %d\n", hdr.mode);
            break;
        }
    }

    /* open the output file */
    if (result) {
        *pfout = fopen(outfile, "wb");
        if (*pfout == NULL) {
            perror("ERROR: fopen() of output file");
            fclose(*pfin);
            *pfin = NULL;
        }
        else {
            result = TRUE;
        }
    }
    return result;
}

int parseCommandLine(int argc, char **argv, char **pinfile, char **poutfile)
{
    int     result;

    result = FALSE;

    if (argc == 3) {
        result = TRUE;
        *pinfile  = argv[1];
        *poutfile = argv[2];
    }
    else {
        fprintf(stderr, "ERROR: Invalid command-line parameters\n");
        fprintf(stderr, "\nUsage: huffon infile outfile\n\n");
    }
    return result;
}

void printParms(int              callno,
                char             *when,
                unsigned long    gpr0,
                unsigned long    gpr1,
                int              cc,
                unsigned char    *fopaddr,
                unsigned long    foplen,
                unsigned char    *sopaddr,
                unsigned long    soplen)
{
    printf("%u: %s call to CMPSC: cc=%u\n", callno, when, cc);
    printf("  GPR 0:    %0*X\n", 2 * sizeof(long), gpr0);
    printf("    Expand:   %s\n", (gpr0 & SC_EXPAND) ? "TRUE" : "False");
    printf("    Format1:  %s\n", (gpr0 & SC_FORMAT1) ? "TRUE" : "False");
    printf("    CDSS:     %u\n", (gpr0 & 0xf000) >> 12);
    printf("    SymTran:  %s\n", (gpr0 & SC_SYMTRAN) ? "TRUE" : "False");
    printf("    ZeroPad:  %s\n", (gpr0 & SC_ZEROPAD) ? "TRUE" : "False");
    printf("    OrderPr:  %s\n", (gpr0 & SC_ORDER) ? "TRUE" : "False");
    printf("    Entropy:  %s\n", (gpr0 & SC_ENTROPY) ? "TRUE" : "False");
    printf("  GPR 1:    %0*X\n", 2 * sizeof(long), gpr1);
    printf("    CmpBit#:  %u\n", gpr1 & 0x07);
    printf("    Offset:   %u bytes\n", ((gpr1 & 0xff8) >> 3) << 7);
    printf("    DictOrg@: %0*X\n", 2 * sizeof(long), gpr1 - (gpr1 & 0xfff));
    printf("  Op1 Addr: %0*X\n", 2 * sizeof(long), (unsigned long)fopaddr);
    printf("  Op1 Len:  %u\n", foplen);
    printf("  Op2 Addr: %0*X\n", 2 * sizeof(long), (unsigned long)sopaddr);
    printf("  Op2 Len:  %u\n", soplen);
}
