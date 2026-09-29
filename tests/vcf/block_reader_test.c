/* Replay consumption and EOF refill at the parser's publication boundary.
 * Interpose unlock only in this translation unit; production needs no hook. */
#include <assert.h>
#include <pthread.h>
#include <stdio.h>

static int consume_on_unlock(pthread_mutex_t *mutex);
#define pthread_mutex_unlock consume_on_unlock
#include "vcf/block_reader.c"
#undef pthread_mutex_unlock

static block_reader *reader;
static bool recycled;

static int consume_on_unlock(pthread_mutex_t *mutex) {
    int result = pthread_mutex_unlock(mutex);
    if (reader == NULL || recycled || reader->n_full == 0) return result;
    /* No other threads run in this test. Model the consumer finishing the
     * published batch, then the reader reusing that slot for EOF, before
     * the parser resumes after unlocking. */
    batch *b = reader->full[reader->full_head];
    assert(b->n == 1 && b->recs[0]->marker.pos == 1000);
    ref_gt_rec_release(b->recs[0]);
    reader->full_head = (reader->full_head + 1) % BLOCK_READER_SLOTS;
    --reader->n_full;
    read_lines(reader, b);  /* have_line is false: refill with EOF */
    reader->read[(reader->read_head + reader->n_read++) % BLOCK_READER_SLOTS] = b;
    recycled = true;
    return result;
}

int main(void) {
    char header[] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tR1";
    char *headers[] = {header};
    const char *line = "20\t1000\t.\tA\tC\t.\tPASS\t.\tGT\t0|1";
    vcf_header h;
    vcf_header_init(&h, "memory", headers, 1, line, strlen(line), NULL);
    reader = util_malloc(sizeof *reader);
    *reader = (block_reader){0};
    reader->header = &h;
    reader->n_threads = 1;
    pthread_mutex_init(&reader->mutex, NULL);
    pthread_cond_init(&reader->changed, NULL);
    batch *b = &reader->slots[0];
    kputs(line, &b->lines[0]);
    index_chrom(line, strlen(line));
    b->n = 1;
    reader->read[0] = b;
    reader->n_read = 1;

    parse_batches(reader);
    assert(recycled);
    /* A parser that rereads b->n after publishing exits before forwarding EOF. */
    assert(reader->n_full == 1 && reader->n_read == 0);
    assert(block_reader_next(reader) == NULL);
    assert(block_reader_next(reader) == NULL);

    pthread_mutex_destroy(&reader->mutex);
    pthread_cond_destroy(&reader->changed);
    free(b->lines[0].s);
    free(reader);
    vcf_header_free(&h);
    puts("PASS parser publishes EOF after its previous batch is recycled");
    return 0;
}
