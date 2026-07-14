/* Dynamic-batching queue for A3b: worker threads keep their exact sequential tube-walk, but
 * each per-tube GPU call becomes bq_submit() -- the thread blocks until a batch of pending
 * requests from many threads is flushed together (one GPU discover+trace batch).
 *
 * The alignment SELECTION logic (alow=eant, alast, triples) is unchanged per thread; only the
 * GPU call is batched, so correctness reduces to "the flush computes each request's result".
 *
 * Flush fires when pending >= flush_n OR every registered worker is waiting (deadlock guard).
 */
#ifndef FGA_BATCH_QUEUE_H
#define FGA_BATCH_QUEUE_H

typedef struct batch_queue batch_queue;

/* Process a batch: for i in [0,n), read reqs[i], write results[i]. Runs under the lock-free
 * flush window (only one thread calls it at a time). ctx is user data (e.g., the gpu_ctx). */
typedef void (*bq_flush_fn)(void **reqs, void **results, int n, void *ctx);

/* capacity = max concurrent submitters (>= number of workers). flush_n = target batch size.
 * flush_us = max microseconds a partial batch waits before flushing (deadlock-free; the first
 * waiter to time out flushes whatever has accumulated). */
batch_queue *bq_create(int capacity, int flush_n, int flush_us, bq_flush_fn fn, void *ctx);
void         bq_destroy(batch_queue *q);

void bq_register(batch_queue *q);     /* a worker announces it will submit */
void bq_deregister(batch_queue *q);   /* a worker is done (may force a final flush) */

/* Submit one request; blocks until the batch containing it is flushed, then returns (result
 * has been filled by the flush fn). */
void bq_submit(batch_queue *q, void *req, void *result);

#endif
