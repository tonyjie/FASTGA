/* See batch_queue.h. Generation-based dynamic batching. */
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include "batch_queue.h"

struct batch_queue {
  pthread_mutex_t mtx;
  pthread_cond_t  cv;
  void          **reqs;
  void          **results;
  int             cap;         // max slots
  int             flush_n;     // target batch size
  long            flush_ns;    // partial-batch timeout (ns)
  int             pending;     // requests in the current (unflushed) batch
  int             registered;  // live workers (informational)
  unsigned long   generation;  // bumped on each flush
  bq_flush_fn     fn;
  void           *ctx;
};

batch_queue *bq_create(int capacity, int flush_n, int flush_us, bq_flush_fn fn, void *ctx)
{ batch_queue *q = (batch_queue *) calloc(1, sizeof(batch_queue));
  if (q == NULL) return NULL;
  pthread_mutex_init(&q->mtx, NULL);
  pthread_cond_init(&q->cv, NULL);
  q->reqs    = (void **) malloc(sizeof(void *) * capacity);
  q->results = (void **) malloc(sizeof(void *) * capacity);
  q->cap     = capacity;
  q->flush_n = flush_n < 1 ? 1 : (flush_n > capacity ? capacity : flush_n);
  q->flush_ns = (long) (flush_us < 1 ? 1 : flush_us) * 1000L;
  q->pending = 0;
  q->registered = 0;
  q->generation = 0;
  q->fn = fn; q->ctx = ctx;
  return q;
}

void bq_destroy(batch_queue *q)
{ if (q == NULL) return;
  pthread_mutex_destroy(&q->mtx);
  pthread_cond_destroy(&q->cv);
  free(q->reqs); free(q->results); free(q);
}

void bq_register(batch_queue *q)
{ pthread_mutex_lock(&q->mtx); q->registered++; pthread_mutex_unlock(&q->mtx); }

//  Run the flush on the current batch (caller holds the lock). Fills all results, resets the
//  batch, bumps the generation, wakes everyone.
static void do_flush(batch_queue *q)
{ int n = q->pending;
  if (n > 0)
    q->fn(q->reqs, q->results, n, q->ctx);
  q->pending = 0;
  q->generation += 1;
  pthread_cond_broadcast(&q->cv);
}

void bq_deregister(batch_queue *q)
{ pthread_mutex_lock(&q->mtx);
  q->registered--;
  pthread_mutex_unlock(&q->mtx);
}

void bq_submit(batch_queue *q, void *req, void *result)
{ pthread_mutex_lock(&q->mtx);
  int slot = q->pending++;
  q->reqs[slot]    = req;
  q->results[slot] = result;
  if (q->pending >= q->flush_n)
    { do_flush(q);                        // batch full: flush now
      pthread_mutex_unlock(&q->mtx);
      return;
    }
  //  Partial batch: wait to be flushed, but with a timeout so the first straggler flushes
  //  whatever accumulated (deadlock-free, and drains the tail).
  unsigned long gen = q->generation;
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  ts.tv_nsec += q->flush_ns;
  if (ts.tv_nsec >= 1000000000L) { ts.tv_sec += ts.tv_nsec / 1000000000L; ts.tv_nsec %= 1000000000L; }
  while (q->generation == gen)
    { int rc = pthread_cond_timedwait(&q->cv, &q->mtx, &ts);
      if (q->generation != gen)          // someone flushed our batch
        break;
      if (rc == ETIMEDOUT)
        { do_flush(q);                    // we time out first -> flush the partial batch
          break;
        }
    }
  pthread_mutex_unlock(&q->mtx);
}
