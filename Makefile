DEST_DIR = ~/bin

CFLAGS = -O3 -Wall -Wextra -Wno-unused-result -fno-strict-aliasing

CC = gcc

ALL = FAtoGDB GDBtoFA GDBstat GDBshow GIXmake GIXshow GIXrm GIXmv GIXcp FastGA ALNshow ALNtoPAF ALNtoPSL ALNreset ALNplot ALNchain PAFtoALN PAFtoPSL ONEview FastKS ONEalnTEST ANOstat ANOshow ANOtoBED BEDtoANO

all: $(ALL)

libfastk.c: gene_core.c gene_core.h
libfastk.h: gene_core.h

GDB.c: gene_core.c gene_core.h
GDB.h: gene_core.h

FAtoGDB: FAtoGDB.c GDB.c GDB.h ANO.c ANO.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FAtoGDB FAtoGDB.c GDB.c ANO.c gene_core.c ONElib.c -lm -lz

GDBtoFA: GDBtoFA.c GDB.c GDB.h ANO.c ANO.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBtoFA GDBtoFA.c GDB.c ANO.c gene_core.c ONElib.c -lm -lz

GDBstat: GDBstat.c GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBstat GDBstat.c GDB.c gene_core.c ONElib.c -lpthread -lm -lz

GDBshow: GDBshow.c GDB.h GDB.c ANO.h ANO.c select.c select.h hash.c hash.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o GDBshow GDBshow.c ONElib.c ANO.c GDB.c select.c hash.c gene_core.c -lpthread -lm -lz

GIXmake: GIXmake.c MSDsort.c libfastk.c libfastk.h ONElib.c ONElib.h ANO.c ANO.h GDB.c GDB.h
	$(CC) $(CFLAGS) -DLCPs -o GIXmake GIXmake.c MSDsort.c libfastk.c ONElib.c ANO.c GDB.c gene_core.c -lpthread -lm -lz

GIXshow: GIXshow.c libfastk.c libfastk.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXshow GIXshow.c libfastk.c gene_core.c -lpthread -lm

GIXrm: GIXrm.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXrm GIXrm.c gene_core.c -lm

GIXmv: GIXxfer.c GDB.c GDB.h ANO.c ANO.h gene_core.c ONElib.c ONElib.h gene_core.h
	$(CC) $(CFLAGS) -DMOVE -o GIXmv GIXxfer.c GDB.c ANO.c ONElib.c gene_core.c -lm -lz

GIXcp: GIXxfer.c GDB.c GDB.h ANO.c ANO.h ONElib.c ONElib.h gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o GIXcp GIXxfer.c GDB.c ANO.c ONElib.c gene_core.c -lm -lz

FastGA: FastGA.c libfastk.c libfastk.h GDB.c GDB.h RSDsort.c align.c align.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FastGA FastGA.c RSDsort.c libfastk.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

# GPU-accelerated FastGA (-G offloads forward-strand Local_Alignment to the A100)
FastGA.gpu: FastGA.c RSDsort.c libfastk.c align.c GDB.c alncode.c gene_core.c ONElib.c gpu/fastga_gpu.cu gpu/fastga_gpu.h gpu/batch_queue.c gpu/batch_queue.h
	nvcc -O3 -arch=sm_80 --default-stream per-thread -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	$(CC) $(CFLAGS) -DGPU -I. -o FastGA.gpu FastGA.c RSDsort.c libfastk.c align.c GDB.c alncode.c gene_core.c ONElib.c gpu/batch_queue.c gpu/fastga_gpu.o -lpthread -lm -lz -L/usr/local/cuda/lib64 -lcudart -lstdc++

FastKS: FastKS.c libfastk.c libfastk.h GDB.c GDB.h RSDsort.c align.c align.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o FastKS FastKS.c RSDsort.c libfastk.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNshow: ALNshow.c align.h align.c GDB.c GDB.h select.c select.h hash.c hash.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNshow ALNshow.c align.c GDB.c alncode.c select.c hash.c gene_core.c ONElib.c -lpthread -lm -lz

ALNtoPAF: ALNtoPAF.c align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPAF ALNtoPAF.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNtoPSL: ALNtoPSL.c align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ALNtoPSL ALNtoPSL.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

extract_tasks: gpu/extract_tasks.c gpu/task_format.h align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -I. -o extract_tasks gpu/extract_tasks.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

extract_disc: gpu/extract_disc.c gpu/disc_format.h align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -I. -o extract_disc gpu/extract_disc.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

extract_trace: gpu/extract_trace.c gpu/trace_format.h align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -I. -o extract_trace gpu/extract_trace.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

extract_seeds: gpu/extract_seeds.c gpu/wave_harness.h align.h align.c GDB.c GDB.h alncode.c alncode.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -I. -o extract_seeds gpu/extract_seeds.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

# GPU discovery fidelity vs local divergence (validation B): reads a .disc, buckets by ref_diffs/len
disc_fidelity: gpu/disc_fidelity.c gpu/disc_format.h gpu/fastga_gpu.h gpu/fastga_gpu.cu
	nvcc -O3 -arch=sm_80 -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	$(CC) $(CFLAGS) -I. -Igpu -o gpu/disc_fidelity gpu/disc_fidelity.c gpu/fastga_gpu.o -L/usr/local/cuda/lib64 -lcudart -lstdc++

# A1: GPU trace-point emission validator (compares GPU trace to Compute_Alignment reference)
trace_validate: gpu/trace_validate.cu gpu/trace_format.h
	nvcc -O3 -arch=sm_80 -Igpu -o gpu/trace_validate gpu/trace_validate.cu

ALNreset: ALNreset.c GDB.c GDB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNreset ALNreset.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

ALNplot: ALNplot.c hash.c hash.h select.c select.h GDB.c GDB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNplot ALNplot.c GDB.c alncode.c select.c hash.c gene_core.c ONElib.c -lpthread -lm -lz

ALNchain: ALNchain.c GDB.c GDB.h ONElib.c ONElib.h alncode.c alncode.h
	$(CC) $(CFLAGS) -o ALNchain ALNchain.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

PAFtoALN: PAFtoALN.c alncode.c hash.c hash.h GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o PAFtoALN PAFtoALN.c GDB.c alncode.c hash.c gene_core.c ONElib.c -lpthread -lm -lz

PAFtoPSL: PAFtoPSL.c gene_core.c gene_core.h
	$(CC) $(CFLAGS) -o PAFtoPSL PAFtoPSL.c gene_core.c -lpthread -lm -lz

ONEview: ONEview.c ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ONEview ONEview.c ONElib.c -lm -lz

ONEalnTEST: ONEaln.c ONEaln.h GDB.c GDB.h ONElib.c ONElib.h align.c align.h alncode.c alncode.h
	$(CC) $(CFLAGS) -DTEST -o ONEalnTEST ONEaln.c GDB.c alncode.c align.c gene_core.c ONElib.c -lm -lz

ANOstat: ANOstat.c ANO.c ANO.h GDB.c GDB.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ANOstat ANOstat.c ANO.c GDB.c gene_core.c ONElib.c -lpthread -lm -lz

ANOshow: ANOshow.c GDB.h GDB.c ANO.h ANO.c select.c select.h hash.c hash.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ANOshow ANOshow.c ONElib.c ANO.c GDB.c select.c hash.c gene_core.c -lpthread -lm -lz

ANOtoBED: ANOtoBED.c GDB.c GDB.h ANO.c ANO.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o ANOtoBED ANOtoBED.c GDB.c ANO.c gene_core.c ONElib.c -lm -lz

BEDtoANO: BEDtoANO.c hash.c hash.h GDB.c GDB.h ANO.c ANO.h ONElib.c ONElib.h
	$(CC) $(CFLAGS) -o BEDtoANO BEDtoANO.c hash.c GDB.c ANO.c gene_core.c ONElib.c -lm -lz

clean:
	rm -f $(ALL)
	rm -fr *.dSYM
	rm -f FastGA.tar.gz

install:
	cp $(ALL) $(DEST_DIR)

package:
	make clean
	tar -zcf FastGA.tar.gz LICENSE README.md Makefile *.h *.c

deploy:
	macdeployqt ALNview.app -dmg

# A2: batched gpu_trace_batch library self-test
trace_lib_test: gpu/trace_lib_test.cu gpu/trace_format.h gpu/fastga_gpu.h gpu/fastga_gpu.cu
	nvcc -O3 -arch=sm_80 -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	nvcc -O3 -arch=sm_80 -Igpu -o gpu/trace_lib_test gpu/trace_lib_test.cu gpu/fastga_gpu.o

# A3b feasibility: time batched gpu_trace_batch over all alignments of a part
trace_bench: gpu/trace_bench.cu gpu/trace_format.h gpu/fastga_gpu.h gpu/fastga_gpu.cu
	nvcc -O3 -arch=sm_80 -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	nvcc -O3 -arch=sm_80 -Igpu -o gpu/trace_bench gpu/trace_bench.cu gpu/fastga_gpu.o

# A3b feasibility: time batched gpu_discover_batch over all discovery tasks of a part
disc_bench: gpu/disc_bench.cu gpu/disc_format.h gpu/fastga_gpu.h gpu/fastga_gpu.cu
	nvcc -O3 -arch=sm_80 -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	nvcc -O3 -arch=sm_80 -Igpu -o gpu/disc_bench gpu/disc_bench.cu gpu/fastga_gpu.o

# Clean wave-vs-wave CPU baseline: Compute_Alignment at N threads on the same .trace tasks
cpu_trace_bench: gpu/cpu_trace_bench.c gpu/trace_format.h align.c align.h GDB.c gene_core.c ONElib.c
	$(CC) $(CFLAGS) -fopenmp -I. -Igpu -o gpu/cpu_trace_bench gpu/cpu_trace_bench.c align.c GDB.c alncode.c gene_core.c ONElib.c -lpthread -lm -lz

# On-device 2-bit decompression bench (the #2 hotspot) vs CPU Uncompress_Read
decomp_bench: gpu/decomp_bench.cu
	nvcc -O3 -arch=sm_80 -Xcompiler -fopenmp -o gpu/decomp_bench gpu/decomp_bench.cu

# A3b batching-queue unit test (pure pthreads, no CUDA)
batch_queue_test: gpu/batch_queue_test.c gpu/batch_queue.c gpu/batch_queue.h
	$(CC) $(CFLAGS) -Igpu -o gpu/batch_queue_test gpu/batch_queue_test.c gpu/batch_queue.c -lpthread

# A3b genome-resident kernel validation (discover_g / trace_g vs per-task)
genome_resident_test: gpu/genome_resident_test.cu gpu/fastga_gpu.cu gpu/fastga_gpu.h gpu/disc_format.h gpu/trace_format.h
	nvcc -O3 -arch=sm_80 -c gpu/fastga_gpu.cu -o gpu/fastga_gpu.o
	nvcc -O3 -arch=sm_80 -Igpu -o gpu/genome_resident_test gpu/genome_resident_test.cu gpu/fastga_gpu.o
