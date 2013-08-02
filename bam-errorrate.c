#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <inttypes.h>
#include <time.h>
#include <unistd.h>

#define TOTAL 0
#define MATCH 1
#define MISMATCH 2
#define AMBIGUOUS 3
#define INSERTION 4
#define DELETION 5
#define READ_COUNT 3
#define TYPE_COUNT 6
#define NDIV 1000000000.0
#define DEF_SEQ_LEN 16384
#define STATUS_MOD 2500000

int main(int argc, char const *argv[]) {
    unsigned int i, j;
    struct timespec ts;
    double start_time;
    int positions = DEF_SEQ_LEN;
    uint64_t count;
    uint64_t* read_counts[READ_COUNT][TYPE_COUNT];

    clock_gettime(CLOCK_REALTIME, &ts);
    start_time = ((double) ts.tv_sec) + (ts.tv_nsec/NDIV);

    count = 0;

// auto resizing of read_counts arrays if read lengths are inconsistent
    for (i = 0; i < READ_COUNT; i++) {
        for (j = 0; j < TYPE_COUNT; j++) {
            read_counts[i][j] = (uint64_t*) calloc(positions, sizeof(uint64_t));
        }
    }

    while (1) {
        char *line, *loc;
        char *cigar, *seq, *qual, *md, *ecigar, *emd;
        long flag;
        int spos, mpos, qpos, mod;
        int read_end, rev_complement;
        unsigned int cpos;
        unsigned int ecigar_len, ecigar_maxlen, emd_len, emd_maxlen;
        int hardclipwarned;
        int badcigardelwarned;
        int line_size;
        size_t nbytes = DEF_SEQ_LEN;
        // unsigned int seq_length;

        count++;
        if (!(count % STATUS_MOD)) {
            clock_gettime(CLOCK_REALTIME, &ts);
            double cur_time = ((double) ts.tv_sec) + (ts.tv_nsec/NDIV);
            // fprintf(stderr, "%" PRIu64 ": %-8.2f reads/sec, %8f seconds elapsed\n", count, (count/(cur_time - start_time)), cur_time - start_time);
            fprintf(stderr, "%12"PRIu64": %12.2lf reads/sec, %8.2lf seconds elapsed\n", count, count/(cur_time - start_time), cur_time - start_time);
        }

        // grab line
        line = (char*) malloc(nbytes+1);
        line_size = getline(&line, &nbytes, stdin);
        if (line_size == -1) {
            fprintf(stderr, "EOF\n");
            free(line);
            break;
        }

        // chomp
        if (line[line_size-1] == '\n') {
            line[line_size-1] = '\0';
        }

        // skip headers
        if (line[0] == '@') {
            free(line);
            continue;
        }

        // split by fields:
        // qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags
        flag = 0; cigar = NULL; seq = NULL; qual = NULL; md = NULL;
        loc = strtok(line,"\t");
        i = 0;
        while (1) {
            i++;
            loc = strtok(NULL, "\t");
            if (loc == NULL) {
                break;
            } else if (i == 1) {
                char * next;
                flag = strtol(loc, &next, 10);
                if (loc == next) {
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "Could not parse flag out of SAM line.\n"); return 1;
                }
            } else if (i == 5) {
                cigar = loc;
            } else if (i == 9) {
                seq = loc;
            } else if (i == 10) {
                qual = loc;
            } else if (i >= 11) {
                // search for MD tag
                loc = strstr(loc,"MD:Z:");
                if (loc != NULL) {
                    md = loc+5;
                }
            }
        }

        //if (!flag) {
        //    fprintf(stderr, "ERROR: %s\n", line);
        //    fprintf(stderr, "Could not parse flag out of SAM line.\n"); return 1;
        //} 

        // We parse the flag first because unmapped reads may not have an MD tag...
        // We skip unmapped reads.
        if (flag & 4) {
            free(line);
            continue; // Skip unmapped reads
        } else if (flag & 256) {
            free(line);
            continue; // Indicates "Secondary alignment"
            // We skip for consistency with the original tool
        } else if (flag & 1024) {
            free(line);
            continue; // Indicates "PCR or optical duplicate"
            // We skip for consistency with the original tool
        }

        // bsmap specific hack
        if (!strcmp(cigar, "0M") && !strcmp(seq, "*") && !strcmp(qual, "*")) {
            free(line);
            continue;
        }

        // Determine the read end
        read_end = 0;
        if (flag & 1) {
            if (flag & 64) {
                read_end = 1;
            } else if (flag & 128) {
                read_end = 2;
            }
        }

        if (!cigar) {
            fprintf(stderr, "ERROR: %s\n", line);
            fprintf(stderr, "Could not parse CIGAR out of SAM line.\n"); return 1;
        } else if (!seq) {
            fprintf(stderr, "ERROR: %s\n", line);
            fprintf(stderr, "Could not parse seq out of SAM line.\n"); return 1;
        } else if (!md) {
            fprintf(stderr, "ERROR: %s\n", line);
            fprintf(stderr, "Could not parse MD out of SAM line.\n"); return 1;
        }

        // TODO: Be VERY careful if you remove the following two tests. Code marked with XXX WILL break if you do so. 
        if (strchr(cigar, 'P')) {
            fprintf(stderr, "ERROR: %s\n", line);
            fprintf(stderr, "P CIGAR operation not implemented\n"); return 1;
        /*
         *} else if (strchr(cigar,'N')) {
         *    fprintf(stderr, "ERROR: %s\n", line);
         *    fprintf(stderr, "N CIGAR operation not implemented\n"); return 1;
         */
        }
        
        
        // parse the cigar string into an expanded cigar string
        // for instance 3M2D5M1S would become MMMDDMMMMMS
        // this makes it faster (or at least easier) to iterate over

        // seq_length = 0;
        ecigar_len = 0;
        ecigar_maxlen = DEF_SEQ_LEN;
        ecigar = (char*) malloc(sizeof(char)*ecigar_maxlen);
        loc = cigar;

        // i'm worried there might be a bug in the following loop...

        while(1) {
            unsigned long c;
            char op;
            char* next = loc;
            c = strtoul(loc, &next, 10);
            if (loc == next) {
                fprintf(stderr, "ERROR: %s\n", line);
                fprintf(stderr, "Error parsing CIGAR string.\n"); return 1;
            }

            op = *next;
            loc = next + 1;

            if (c == 0) {
                fprintf(stderr, "ERROR: %s\n", line);
                fprintf(stderr, "Error parsing CIGAR string.\n"); return 1;
            }

            if (op != 'M' &&
                op != 'I' &&
                op != 'D' &&
                op != 'N' &&
                op != 'S' &&
                op != 'H' &&
                op != 'P' &&
                op != '=' &&
                op != 'X'
            ) {
                fprintf(stderr, "ERROR: %s\n", line);
                fprintf(stderr, "Invalid OP found in CIGAR.\n"); return 1;
            }

            // XXX I think we can just ignore N, it seems to represent an intron and it isn't reflected in the MD tag
            if (op == 'N') {
                continue;
            }

            // fill in the ecigar
            for (i = 0; i < c; i++) {
                if (ecigar_len == (ecigar_maxlen - 1)) {
                    ecigar_maxlen += 10;
                    ecigar = (char*)realloc(ecigar, ecigar_maxlen);
                }
                ecigar[ecigar_len] = op;
                ecigar_len++;
            }

            // keep a running sequence length
            //if (op != 'D') {
            //    seq_length += c;
            //}

            if (*loc == '\0') {
                ecigar[ecigar_len] = '\0';
                break;
            }
        }

        if (ecigar_len != strlen(ecigar)) {
            fprintf(stderr, "ERROR: %s\n", line);
            fprintf(stderr, "ecigar_len != strlen(ecigar)\n"); return 1;
        }

        // we also turn the MD tag into an expanded MD tag

        emd_len = 0;
        emd_maxlen = DEF_SEQ_LEN;
        emd = (char*) malloc(sizeof(char)*emd_maxlen);
        loc = md;
        int should_be_num = 1;

        while(1) {
            unsigned long c;
            char op;
            char* next;

            if (should_be_num) {
                c = strtoul(loc, &next, 10);
                if (loc == next) {
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "Error parsing MD tag.\n"); return 1;
                }
                op = '=';

                // make sure a 0 from strtoul is really a 0
                if (c == 0 && *loc != '0') {
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "Expected nonzero match in MD tag.\n"); return 1;
                }
                loc = next;

                should_be_num = 0;
            } else {
                if (*loc == '^') {
                    loc++;
                    op = 'D';
                } else {
                    op = 'X';
                }

                c = 0;
                while(1) {
                    if (*loc == '\0') {
                        fprintf(stderr, "ERROR: %s\n", line);
                        fprintf(stderr, "MD tag ended on a base.\n"); return 1;
                    } else if (
                        *loc == '0' ||
                        *loc == '1' ||
                        *loc == '2' ||
                        *loc == '3' ||
                        *loc == '4' ||
                        *loc == '5' ||
                        *loc == '6' ||
                        *loc == '7' ||
                        *loc == '8' ||
                        *loc == '9'
                    ) {
                        break;
                    }
                    loc++;
                    c++;
                }

                // this is unusual but would probably happen if the MD tag was like 5^5
                if (c == 0) {
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "Expected a base but didn't find one.\n"); return 1;
                }
                should_be_num = 1;
            }

            for (i = 0; i < c; i++) {
                if (emd_len == (emd_maxlen - 1)) {
                    emd_maxlen += 10;
                    emd = (char*)realloc(emd, emd_maxlen);
                }
                emd[emd_len] = op;
                emd_len++;
            }

            if (*loc == '\0') {
                emd[emd_len] = '\0';
                break;
            }
        }

        // now we're ready to do the counts

        rev_complement = flag & 16;

        // NOTE: we made changes to the pileup version; go compare it
        // TODO: not sure if hard clipping is implemented correctly
        
        
        // The positions within the cigar string, md tag, read, and sequence do not match up with each other.
        //      For instance, if we encounter a deletion in the cigar string, we *don't* increment our position
        //      within the sequence or the read. If the sequence has a 10bp hard clipp, then the position within
        //      the sequence will be 10 off from the position in the actual read.
        // spos: the position within the sequence from the sam file
        // qpos: the position in terms of the actual read; this is the position used for recording counts
        // cpos: position within the expander cigar string
        // mpos: position within the expanded MD tag
        spos = 0;
        mpos = 0;
        if (rev_complement) {
            qpos = strlen(seq) - 1;
            mod = -1;
        } else {
            qpos = 0;
            mod = 1;
        }
        hardclipwarned    = 0;
        badcigardelwarned = 0;

        for (cpos = 0; cpos < ecigar_len; cpos++) {
            if (ecigar[cpos] == 'S') {
                qpos += mod;
                spos++;
            } else if (ecigar[cpos] == 'H') {
                qpos += mod;
                // TODO not sure if hard clipping is implemented correctly
                if (hardclipwarned == 0) {
                    hardclipwarned = 1;
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "Support for hard clipping untested\n");
                }
            } else if (ecigar[cpos] == 'I') {
                read_counts[read_end][INSERTION][qpos]++;
                read_counts[read_end][TOTAL][qpos]++;
                qpos += mod;
                spos++;
            } else if (ecigar[cpos] == 'D') {
                if (qpos - mod < 0 || qpos - mod >= 16384) {
                    if (badcigardelwarned == 0) {
                        badcigardelwarned = 1;
                        fprintf(stderr, "ERROR: %s\n", line);
                        fprintf(stderr, "Detected a deletion at the beginning of a CIGAR string.\n");
                        fprintf(stderr, "This might indicate a poorly constructed bam.\n");
                        fprintf(stderr, "Only warning about this once. Diagnostic info follows:\n");
                        fprintf(stderr, "(spos=%d, qpos=%d, cpos=%d, mpos=%d)\n", spos, qpos, cpos, mpos);
                        fprintf(stderr, "(qpos - mod out of bounds; %d - %d = %d)\n", qpos, mod, qpos-mod);
                    }
                } else {
                    read_counts[read_end][DELETION][qpos-mod]++; // TODO use previous qpos for consistency?
                }
                if (!(emd[mpos] == 'D')) {
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "CIGAR indicated a deletion while MD did not.\n"); return 1;
                }
                mpos++;
            } else if (toupper(seq[spos]) == 'N') {
                read_counts[read_end][AMBIGUOUS][qpos]++;
                read_counts[read_end][TOTAL][qpos]++;
                qpos += mod;
                spos++;
                mpos++;
            } else if (emd[mpos] == 'X') {
                read_counts[read_end][MISMATCH][qpos]++;
                read_counts[read_end][TOTAL][qpos]++;
                // We know CIGAR op isn't I, D, S, or H (N and P are not implemented)
                // So CIGAR op must be M, =, or X. Since MD indicates a mismatch, die if CIGAR indicates a match
                // TODO replace $sam{qname} from perl version with actual variable
                if (ecigar[cpos] == '=') { // XXX
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "MD indicated a mismatch while CIGAR indicated a match in $sam{qname}\n"); return 1;
                }
                qpos += mod;
                spos++;
                mpos++;
            } else if (emd[mpos] == '=') {
                read_counts[read_end][MATCH][qpos]++;
                read_counts[read_end][TOTAL][qpos]++;
                // We know CIGAR op isn't I, D, S, or H (N and P are not implemented)
                // So CIGAR op must be M, =, or X. Since MD indicates a mismatch, die if CIGAR indicates a match
                // TODO replace $sam{qname} from perl version with actual variable
                if (ecigar[cpos] == 'X') { // XXX
                    fprintf(stderr, "ERROR: %s\n", line);
                    fprintf(stderr, "MD indicated a match while CIGAR indicated a mismatch in $sam{qname}\n"); return 1;
                }
                qpos += mod;
                spos++;
                mpos++;
            } else {
                // fprintf(stderr, "c %c, m %c, s %c, %s, %s, %s, %s\n", ecigar[cpos], emd[mpos], seq[spos], cigar, ecigar, md, emd);
                fprintf(stderr, "ERROR: %s\n", line);
                fprintf(stderr, "MD tag was expected to indicate a match or mismatch but it indicated neither\n"); return 1;
            }

        }

        free(line);
        free(ecigar);
        free(emd);
    }

    int re, pos;

    if (argc > 0) {
        printf("# bam-errorate; command executed: %s\n", argv[0]);
    }
    printf("read_end\tposition\ttotal\tmatch\terror\terror_rate\tmismatch\tmismatch_rate\tambiguous\tambiguous_rate\tinsertion\tinsertion_rate\tdeletion\tdeletion_rate\n");

    uint64_t total_read_count;
    for (pos = positions-1; pos >= 0; pos--) {
        total_read_count = 0;
        for (re = 0; re < READ_COUNT; re++) {
            for (i = 0; i < TYPE_COUNT; i++) {
                total_read_count += read_counts[re][i][pos];
            }
        }
        positions = pos;
        if (!(total_read_count == 0)) {
            break;
        }
    }
    for (re = 0; re < READ_COUNT; re++) {
        uint64_t sum_total, sum_total_plusdel, sum_position_count, sum_match, sum_mismatch, sum_ambiguous, sum_insertion, sum_deletion, sum_error;

        total_read_count = 0;

        for (i = 0; i < TYPE_COUNT; i++) {
            for (pos = 0; pos < positions; pos++) {
                total_read_count += read_counts[re][i][pos];
            }
        }

        if (total_read_count == 0) {
            continue; // we found nothing for this read end
        }

        sum_total          = 0;
        sum_total_plusdel  = 0;
        sum_position_count = 0;
        sum_match          = 0;
        sum_mismatch       = 0;
        sum_ambiguous      = 0;
        sum_insertion      = 0;
        sum_deletion       = 0;
        sum_error          = 0;

        for (pos = 0; pos <= positions; pos++) {
            uint64_t total, match, mismatch, ambiguous, insertion, deletion, error;
            double error_rate, mismatch_rate, ambiguous_rate, insertion_rate, deletion_rate;
            uint64_t position_count;
            position_count = read_counts[re][TOTAL][pos];
            if (position_count == 0) {
                printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", re, pos, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
                continue;
            }

            match     = read_counts[re][MATCH][pos];
            mismatch  = read_counts[re][MISMATCH][pos];
            ambiguous = read_counts[re][AMBIGUOUS][pos];
            insertion = read_counts[re][INSERTION][pos];
            deletion  = read_counts[re][DELETION][pos];

            total         = match + mismatch + ambiguous + insertion;
            if (total != position_count) {
                fprintf(stderr, "Position count did not match total = match + mismatch + amiguous + insertion.\n"); return 1;
            }
            // total_plusdel = match + mismatch + ambiguous + insertion + deletion;
            error         = mismatch + ambiguous + insertion + deletion;

            error_rate     = (double) error / (double) total;
            mismatch_rate  = (double) mismatch / (double) total;
            ambiguous_rate = (double) ambiguous / (double) total;
            insertion_rate = (double) insertion / (double) total;
            deletion_rate  = (double) deletion / (double) total;

            sum_total          += total;
            sum_match          += match;
            sum_mismatch       += mismatch;
            sum_ambiguous      += ambiguous;
            sum_insertion      += insertion;
            sum_deletion       += deletion;
            sum_error          += error;

            printf("%d\t%d\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\n",
                re, pos, total, match, error, error_rate, mismatch, mismatch_rate, ambiguous, ambiguous_rate, insertion, insertion_rate, deletion, deletion_rate);
        }

        double error_rate, mismatch_rate, ambiguous_rate, insertion_rate, deletion_rate;
        error_rate     = (double) sum_error / (double) sum_total;
        mismatch_rate  = (double) sum_mismatch / (double) sum_total;
        ambiguous_rate = (double) sum_ambiguous / (double) sum_total;
        insertion_rate = (double) sum_insertion / (double) sum_total;
        deletion_rate  = (double) sum_deletion / (double) sum_total;

        printf("%d\t%s\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\t%"PRIu64"\t%.17lf\n",
                re, "SUM", sum_total, sum_match, sum_error, error_rate, sum_mismatch, mismatch_rate, sum_ambiguous, ambiguous_rate, sum_insertion, insertion_rate, sum_deletion, deletion_rate);
    }


    for (i = 0; i < READ_COUNT; i++) {
        for (j = 0; j < TYPE_COUNT; j++) {
            free(read_counts[i][j]);
        }
    }

    return 0;
}


// close $bam_fh;

// my @headers = qw/read_end position total match error error_rate mismatch mismatch_rate ambiguous ambiguous_rate insertion insertion_rate deletion deletion_rate/;
// print join("\t", @headers) . "\n";

