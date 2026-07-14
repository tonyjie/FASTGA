/* Unit test for the A3b dynamic-batching queue: N threads each submit M requests sequentially
 * (blocking per request); the flush doubles each. Verifies every result, that batching really
 * happens (flush sizes > 1), and that the tail drains with no deadlock. Build: make batch_queue_test
 */
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include "batch_queue.h"

#define NTHREADS 64
#define NREQ     500
#define FLUSH_N  16

static long flush_calls = 0, flush_items = 0, flush_max = 0, flush_multi = 0;
static pthread_mutex_t stat_mtx = PTHREAD_MUTEX_INITIALIZER;

// req and result are int* ; result = 2*req
static void flush_fn(void **reqs, void **results, int n, void *ctx)
{ (void) ctx;
  for (int i = 0; i < n; i++)
    *((int *) results[i]) = 2 * (*((int *) reqs[i]));
  pthread_mutex_lock(&stat_mtx);
  flush_calls++; flush_items += n; if (n > flush_max) flush_max = n; if (n > 1) flush_multi++;
  pthread_mutex_unlock(&stat_mtx);
}

static batch_queue *Q;
static long bad = 0;
static pthread_mutex_t bad_mtx = PTHREAD_MUTEX_INITIALIZER;

static void *worker(void *arg)
{ int base = *(int *) arg;
  bq_register(Q);
  for (int i = 0; i < NREQ; i++)
    { int req = base * NREQ + i, res = -1;
      bq_submit(Q, &req, &res);
      if (res != 2 * req)
        { pthread_mutex_lock(&bad_mtx); bad++; pthread_mutex_unlock(&bad_mtx); }
    }
  bq_deregister(Q);
  return NULL;
}

int main(void)
{ pthread_t th[NTHREADS]; int id[NTHREADS];
  Q = bq_create(NTHREADS, FLUSH_N, 2000 /*us*/, flush_fn, NULL);
  for (int t = 0; t < NTHREADS; t++) { id[t] = t; pthread_create(&th[t], NULL, worker, &id[t]); }
  for (int t = 0; t < NTHREADS; t++) pthread_join(th[t], NULL);
  bq_destroy(Q);

  long expect = (long) NTHREADS * NREQ;
  printf("requests: %ld   wrong results: %ld\n", expect, bad);
  printf("flush calls: %ld   items: %ld (ave batch %.1f, max %ld, batches>1: %ld)\n",
         flush_calls, flush_items, (double) flush_items / flush_calls, flush_max, flush_multi);
  int ok = (bad == 0) && (flush_items == expect) && (flush_max > 1);
  printf("%s\n", ok ? "PASS: all correct, batching happened, no deadlock" : "FAIL");
  return ok ? 0 : 1;
}
