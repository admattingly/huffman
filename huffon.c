/* huffon.c */

/*
 * Use the entropy encoding feature of the CMPSC instruction to compress a file
 * using entropy encoding only - no Lempel-Ziv dictionary of multi-character
 * sequences.
 *
 * To compile under z/OS UNIX:
 *
 * $ xlc -qasm -qasmlib=sys1.maclib -qin=all:nostp -ohuffon huffon.c
 */

#define _LARGE_TIME_API
#define _POSIX_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

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
#pragma pack(reset)

/* package structure for Package Merge algorithm */
#pragma pack(packed)
typedef struct _PACKAGE {
    int     freq;
    int     symbol;
} PACKAGE, *PPACKAGE;
#pragma pack(reset)

/* private function prototypes */
static int buildDictionary(unsigned char *buffer, int inbytes, unsigned char **pcdict, unsigned int *pcdictlen, unsigned char *alphabet);
static int compareFreq(const void *element1, const void *element2);
static int compress(FILE *fin, FILE *fout);
static int copy(FILE *fin, FILE *fout);
static void free4k(void *addr, unsigned int size);
static void *malloc4k(unsigned int size);
static int openFiles(char *infile, char *outfile, FILE **pfin, FILE **pfout, int *pmode);
static int packageMerge(int *freq, int n, short *ed);
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
                if (compress(fin, fout)) {
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

int buildDictionary(unsigned char *buffer, int inbytes, unsigned char **pcdict, unsigned int *pcdictlen, unsigned char *alphabet)
{
    int             i, j, k, result, sum;
    const int       n = 256;            // dictionary has 256 entries
    int             freq[n + 1];        // need a slot for the EOR marker
    short           ed[MAXBITS];        // entity descriptor
    unsigned char   *cdict;             // compression dictionary
    unsigned int    cdictlen;
    PPACKAGE        symbol;
    unsigned short  *stt;               // symbol translation table entries
    unsigned short  codeword;

    result = TRUE;
    symbol = NULL;
    cdict  = NULL;

    /* compute character frequencies, based on first buffer load read from input file */
    memset(freq, 0x00, sizeof(freq));   // zero out the frequency array
    for (i = 0; i < inbytes; i++) {
        freq[(int)buffer[i]]++;
    }
    freq[n] = -1;   // assign lowest frequency to EOR marker

    /* print out character frequencies */
//    for (i = 0; i < n; i++) {
//        printf("Symbol: %3d, Frequency: %5d\n", i, freq[i]);
//    }

    /* run a package merge to produce an entropy descriptor */
    result = packageMerge(freq, n + 1, ed);

    /* print out entity descriptor */
    if (result) {
        sum = 0;
        printf("Entity descriptor:\n");
        printf("Bits Count\n");
        printf("---- -----\n");
        for (i = 0; i < MAXBITS; i++) {
            printf("%4d %5d\n", i + 1, ed[i]);
            sum += ed[i];
        }
        printf("     =====\n");
        printf("     %5d\n", sum);
    }

    /* allocate space for dictionary on 4KB page boundary, as required by CMPSC */
    cdictlen = n * 8 + sizeof(ed) + n * 2;  // room for classic dictionary, ED and STT
    cdict = malloc4k(cdictlen);
    if (cdict == NULL) {
        result = FALSE;
        fprintf(stderr, "ERROR: malloc4k() of compression dictionary failed.");
    }

    /* classic dictionary comprises 256-character alphabet where none has children */
    memset(cdict, 0x00, 8 * n);

    /* copy entity descriptor after classic dictionary */
    memcpy(cdict + 8 * n, ed, sizeof(ed));

    /* sort symbols by frequency for assignment of codewords */
    symbol = (PPACKAGE)malloc(n * sizeof(PACKAGE));
    if (symbol == NULL) {
        result = FALSE;
        perror("ERROR: malloc() of symbols");
    }
    else {
        for (i = 0; i < n; i++) {
            symbol[i].symbol = i;
            symbol[i].freq   = freq[i];
        }
        /* sort into ascending frequency sequence */
        qsort(symbol, n, sizeof(PACKAGE), compareFreq);
    }

    /* assign codewords */
    stt = (unsigned short *)(void *)(cdict + 8 * n + sizeof(ed));
    codeword = 2;
    i = n - 1;  // index of most frequent symbol in 'symbol' array
    for (j = 0; j < MAXBITS; j++) {
        if (ed[j] > 0) {
            for (k = 0; k < ed[j]; k++) {
                codeword--;
                if (i >= 0) {   // EOR marker is accounted for in ED, but doesn't appear in STT
                    stt[symbol[i].symbol] = (unsigned short)(codeword << (MAXBITS - j - 1));
//                    printf("Symbol: %d, Codeword: %04X\n", symbol[i].symbol, stt[symbol[i].symbol]);
                    /* we also need a descending-frequency-sequence alphabet for expansion */
                    alphabet[n-1-i] = (unsigned char)symbol[i].symbol;
                    i--;
                }
            }
        }
        codeword *= 2;
    }

    /* tidy up */
    if (symbol != NULL) {
        free(symbol);
    }
    if (!result) {
        if (cdict != NULL) {
            free4k(cdict, cdictlen);
            cdict = NULL;   // in case caller tries to free it again
        }
    }

    *pcdict    = cdict;
    *pcdictlen = cdictlen;

    return result;
}

int compareFreq(const void *element1, const void *element2)
{
    const PACKAGE   *p1 = element1;
    const PACKAGE   *p2 = element2;

    /* ascending frequency order and within that, leaves come before parents */
    return (p1->freq < p2->freq ? -1 : (p1->freq > p2->freq ? 1 :
            (p1->symbol > p2->symbol ? -1 : (p1->symbol < p2->symbol ? 1 : 0))));
}

int compress(FILE *fin, FILE *fout)
{
    int             inbytes, outbytes, intotal, outtotal, ncalls, result;
    unsigned char   *inbuffer, *outbuffer, *cdict, alphabet[256];
    unsigned int    cdictlen;
    PHUFF_HEADER    phh;
    int             cbn, cc, offset, foplen, soplen;
    unsigned char   *fopaddr, *sopaddr;
    long            gpr0, gpr1;
    const short     eor = 0;

    result    = TRUE;
    inbuffer  = NULL;
    outbuffer = NULL;
    intotal   = 0;
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
        /* write the compressed file header */
        phh = (PHUFF_HEADER)(void *)inbuffer;
        memcpy(phh->eyecatcher, HUFF_EYECATCHER, HUFF_ECLENGTH);
        phh->version = HUFF_VERSION;
        phh->mode    = MODE_ENTROPY;
        outbytes = fwrite(inbuffer, 1, sizeof(HUFF_HEADER), fout);
        outtotal += outbytes;
        if (outbytes != sizeof(HUFF_HEADER)) {
            result = FALSE;
            perror("ERROR: fwrite() of compressed file header");
        }
    }

    /* perform compression operation */
    if (result) {
        cbn     = 0;    // compressed output starts on a byte boundary
        ncalls  = 0;
        fopaddr = outbuffer;
        foplen  = BUFFER_SIZE;
        while (result && !feof(fin)) {
            inbytes = fread(inbuffer, 1, BUFFER_SIZE, fin);
            if (inbytes > 0) {
                if (intotal == 0) {
                    result = buildDictionary(inbuffer, inbytes, &cdict, &cdictlen, alphabet);
                    if (result) {
                        /* write out the ED to the output file */
                        outbytes = fwrite(cdict + 2048, 1, MAXBITS * sizeof(short), fout);
                        outtotal += outbytes;
                        if (outbytes != MAXBITS * sizeof(short)) {
                            result = FALSE;
                            perror("ERROR: fwrite() of ED");
                        }
                        else {
                            /* write out the descending-frequency-sequence alphabet */
                            outbytes = fwrite(alphabet, 1, 256, fout);
                            outtotal += outbytes;
                            if (outbytes != 256) {
                                result = FALSE;
                                perror("ERROR: fwrite() of alphabet");
                            }
                        }
                    }
                }
                if (result) {   /* we have valid compression dictionary */
                    intotal += inbytes;
                    sopaddr = inbuffer;
                    soplen  = inbytes;
                    cc      = 3;
                    offset  = 2048;     // ED starts at 8 * 256 = 2048
                    while (cc == 3) {
                        gpr0 = SC_ENTROPY;
                        gpr1 = ((unsigned long)cdict) + (((offset % 65536) / 128) << 3) + cbn;
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
        }
        if (result) {
            if (cc == 0 && fopaddr > outbuffer) {  /* some output left to write */
                outbytes = fwrite(outbuffer, 1, fopaddr - outbuffer, fout);
                outtotal += outbytes;
                if (outbytes != fopaddr - outbuffer) {
                    result = FALSE;
                    perror("ERROR: fwrite() to output file");
                }
                else {
                    if (cbn > 0) { /* incomplete byte at end of compressed output */
                        printf("INFO: Incomplete byte at end of output buffer - writing one extra byte\n");
                        outbytes = fwrite(fopaddr, 1, 1, fout);
                        outtotal += outbytes;
                        if (outbytes != 1) {
                            result = FALSE;
                            perror("ERROR: fwrite() of incomplete byte to output file");
                        }
                    }
                }
            }
            if (cbn > 0) { /* write an end-of-file marker for compression with entropy-encoding*/
                printf("INFO: Compressing with entropy encoding and last byte incomplete - writing X'0000' to mark end of file\n");
                outbytes = fwrite(&eor, 1, 2, fout);
                outtotal += outbytes;
                if (outbytes != 2) {
                    result = FALSE;
                    perror("ERROR: fwrite() of End-Of-Record marker to output file");
                }
            }
            else {
                printf("INFO: Compressing with entropy encoding and last byte complete - do nothing\n");
            }
        }
    }

    /* produce a report on space savings */
    if (result) {
        printf("Input file size:  %10d bytes\n", intotal);
        printf("Output file size: %10d bytes\n", outtotal);
        printf("Space reduction:  %10.1f %%\n", 100.0 * (1.0 - ((double)outtotal) / ((double)intotal)));
    }

    /* tidy up */
    if (inbuffer != NULL) {
        free(inbuffer);
    }
    if (outbuffer != NULL) {
        free(outbuffer);
    }
    if (cdict != NULL) {
        free4k(cdict, cdictlen);
    }

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
    int             rc, result;
    struct stat64   st;

    result = FALSE;
    *pmode = MODE_COPY;

    /* get size of input file */
    rc = stat64(infile, &st);
    if (rc != 0) {
        perror("ERROR: stat() of input file");
    }
    else {
        if (st.st_size >= MINIMUM_SIZE) {
            *pmode = MODE_ENTROPY;  /* file is large enough for entropy encoding to by worthwhile */
        }
        *pfin = fopen(infile, "rb");
        if (*pfin == NULL) {
            perror("ERROR: fopen() of input file");
        }
        else {
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
    }
    return result;
}

int packageMerge(int *freq, int n, short *ed)
{
    int         count, i, j, result;
    PPACKAGE    pkg[MAXBITS];
    int         pkglen[MAXBITS], solution[MAXBITS];
    int         *bitlength;

    result = TRUE;

    /* load the first package list */
    pkglen[0] = n;
    pkg[0] = (PPACKAGE)malloc(n * sizeof(PACKAGE));
    if (pkg[0] == NULL) {
        result = FALSE;
        perror("ERROR: malloc() of pkg[0]");
    }
    else {
        for (i = 0; i < n; i++) {
            pkg[0][i].freq   = freq[i];
            pkg[0][i].symbol = i;
        }
        /* put list in ascending frequency sequence */
        qsort(pkg[0], n, sizeof(PACKAGE), compareFreq);
    }

    /* construct remaining package lists */
    for (j = 1; j < MAXBITS && result; j++) {
        pkglen[j] = pkglen[j-1] / 2 + n;
        pkg[j] = (PPACKAGE)malloc(pkglen[j] * sizeof(PACKAGE));
        if (pkg[j] == NULL) {
            result = FALSE;
            perror("ERROR: malloc() of pkg[j]");
        }
        else {
            /* pair previous list's packages into parents */
            for (i = 0; i < pkglen[j-1] / 2; i++) {
                pkg[j][i].freq   = pkg[j-1][2*i].freq + pkg[j-1][2*i+1].freq;
                pkg[j][i].symbol = -1;  // a parent node
            }
            /* add original leaf nodes */
            for (i = 0; i < n; i++) {
                pkg[j][pkglen[j]-n+i] = pkg[0][i];
            }
            /* put list in ascending frequency sequence */
            qsort(pkg[j], pkglen[j], sizeof(PACKAGE), compareFreq);
        }
    }

    /* compute solution lengths */
    if (result) {
        solution[MAXBITS-1] = 2 * n - 2;
        for (j = MAXBITS - 2; j >= 0; j--) {
            count = 0;
            for (i = 0; i < solution[j+1]; i++) {
                if (pkg[j+1][i].symbol == -1) {
                    count++;
                }
            }
            solution[j] = 2 * count;
        }
    }

    /* calculate bit length for each symbol */
    if (result) {
        bitlength = (int *)malloc(n * sizeof(int));
        if (bitlength == NULL) {
            result = FALSE;
            perror("ERROR: malloc() of bitlength[]");
        }
        else {
            memset(bitlength, 0x00, n * sizeof(int));
            for (j = 0; j < MAXBITS; j++) {
                for (i = 0; i < solution[j]; i++) {
                    if (pkg[j][i].symbol >= 0) {
                        bitlength[pkg[j][i].symbol]++;
                    }
                }
            }
        }
    }

    /* compute entropy descriptor from symbol bit lengths */
    if (result) {
        memset(ed, 0x00, MAXBITS * sizeof(short));
        for (i = 0; i < n; i++) {
            ed[bitlength[i]-1]++;
        }
    }

    /* tidy up */
    for (j = 0; j < MAXBITS; j++) {
        if (pkg[j] != NULL) {
            free(pkg[j]);
        }
    }
    if (bitlength != NULL) {
        free(bitlength);
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
