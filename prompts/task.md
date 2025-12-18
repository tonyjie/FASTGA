# Task
Optimize the performance of FastGA algorithm on multi-threaded CPUs without breaking the existing functinality. 

# Problem
FastGA is a pairwise whole genome aligner. It searches for all local DNA alignments between two high quality genomes. It is a highly optimized algorithm that is designed to be fast and efficient. 
The core assumption is that the genomes are nearly complete involving at most several thousand contigs with a sequence quality of Q40 or better. Based on a novel adaptive seed finding algorithm and the first wave-based local aligner developed for daligner (2012), the tool can for example compare two 2Gbp bat genomes finding almost all regions over 100bp that are 70% or more similar in about 5.0 minutes wall clock time on my MacPro with 8 cores (about 28 CPU minutes). Moreover, it uses a trace point concept to record all the found alignments in a compressed and indexable ONEcode file in a very space-efficient manner, e.g. just 44.5MB for over 635,000 local alignments in our running example. These trace point encodings of the alignments can then be swiftly translated into .psl or .paf format on demand with programs provided here.

Using FastGA can be as simple as calling it with two FASTA files containing genome assemblies where each entry is a scaffold with runs of N's separating and potentially giving the estimated distance between the contigs thereof. By default a PAF file encoding all the local alignments found between the two genomes is streamed to the standard output.

We believe there are still opportunities to optimize the performance of FastGA algorithm, e.g. better utilizing the multi-threaded CPUs, changing the data structures to be more cache-friendly, etc. Your task is to find the best optimization strategy to achieve the best performance without breaking the existing functinality.

# Relavant Code Files
The source code of FastGA is these `.c` and `.h` files located in `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA`. You will try to modify these files to optimize the performance of the FastGA algorithm.

You should carefully read the code to understand the algorithm and the data structures used in the code. Note that there are proprocessing steps including `FAtoGDB` and `GIXmake`, but they are skipped in our testing because we already cached the `.1gdb` and `.gix` files in the `EXAMPLE` directory. They are not the bottleneck of the FastGA algorithm.


# Directories
All experiments and generated files should be saved under: /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_1218_1422


# Tools provided
## Script to run the FastGA algorithm
Run the FastGA algorithm with 32 threads on a sample dataset, save the alignment results, and save the log file including the run time and CPU utilization information. 

Script: `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/scripts/run_fastga.sh`

Usage: `bash run_fastga.sh <output_directory>`, where <output_directory> is the directory to save the alignment results and the log file.

Example: `bash run_fastga.sh /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_1218_1422/iter0`

It will run the FastGA binary with 32 threads, take two FASTA files as input, and save the alignment results in `.1aln` file, and save the output log file including the run time and CPU utilization information. An example log file is provided in `prompts/example_output/fastga.log`.
Note that `.1gdb` and `.gix` files are already generated in the `EXAMPLE` directory, so `FAtoGDB` and `GIXmake` are not needed to be run again in our testing. This is because the `FastGA` step is the most time-consuming step, and we want to focus on the optimization of the `FastGA` step. 

The final run time can be grasped from the `wall time` in the log file. For example, in the example `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/prompts/example_output/fastga.log`, it includes the following information: 
```
/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/FastGA -vk -v -T32 -L:/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_12_28_13_18/fastga.log -P/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_12_28_13_18 -1:/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_12_28_13_18/H1vH2.1aln /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE/HAP1.fasta.gz /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/EXAMPLE/HAP2.fasta.gz

  Using 32 threads

  Adaptive seed merge for G1

  Total seeds = 51082720, ave. len = 28.5, seeds per genome position = 0.6

  Resources for phase:  4.122u  29.586s  4.686w  719.3%

  Seed sort and alignment search, 6 parts

  Total hits over 85bp = 338700, 376250 aln's, 323569 non-redundant aln's of ave len 1953

  Sorting and merging alignments

  Resources for phase:  1:49.883u  2.423s  8.428w  1332.5%

  Total Resources:  1:54.008u  32.026s  13.133w  1111.9%  0MB
```
The last line `Total Resources:  1:54.008u  32.026s  13.133w  1111.9%  0MB` shows 1minute 54.008 seconds user CPU time, 32.026 seconds system CPU time, and 13.133 seconds wall clock time with 1111.9% CPU utilization. The 13.133 seconds wall clock time is the final run time of the FastGA algorithm that we aim to optimize. 

## Script to verify the alignment output
Compare two `.1aln` files by hashing their sorted PAF outputs.
Script: `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/scripts/verify_exact_match.sh`

Usage: `bash verify_exact_match.sh --baseline <file> --candidate <file>`, where <file> is the `.1aln` file to compare. The baseline file is already generated in `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/prompts/example_output/H1vH2.1aln`. The candidate file is the `.1aln` file newly generated to compare with the baseline file.

Example: `bash verify_exact_match.sh --baseline /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/prompts/example_output/H1vH2.1aln --candidate /work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/experiments/exp_1218_1422/iter0/H1vH2.1aln`

The success of the verification will be indicated by the exit code of the script. 0 means the two `.1aln` files are exactly the same, 1 means the two `.1aln` files are different, and 2 means the script failed to run. A success verification will show the following output: 
```
Comparing: H1vH2.1aln vs H1vH2.1aln
OK: PAF (no CIGAR) - 81e8fb2eb608ab276ca7b858df35cedc097eb23a97c33bf6ad2d7a5233effc21
PASS: All checks passed
```

# Workflow

0. **Preparation (Git Setup)**
   - Before starting, create and switch to a new git branch named `agent-optimize-1218_1422`.
   - command: `git checkout -b agent-optimize-1218_1422`
   - If the branch exists, handle the error or use a unique name.

1. **Baseline**
   - Run the script `scripts/run_fastga.sh` to get baseline performance.
   - Create a file `baseline/result.txt` with the timing.

2. **Optimization iterations (5 iterations)**
   - Loop from iter_0 to iter_4:
     - Create a new directory `experiments/exp_1218_1422/iter_X` for each iteration. All the generated files in this iteration should be saved in this directory.
     - **Modify Code**: Apply code changes to relavant C files (`.c`, `.h`) in FastGA repository to optimize performance.
     - **Compile**: Run `make -j 32` in `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/`. Ensure it compiles successfully. If the compilation fails, fix the errors and try again.
     - **Run Modified FastGA**: Run the `scripts/run_fastga.sh` script to run the modified FastGA algorithm.
     - **Verify the alignment output**: Run the script `scripts/verify_exact_match.sh` to verify the alignment output. If the verification fails, revert the changes and try a different optimization strategy.
     - **Record & Compare**: Read the wall clock time from the output log file. Compare with the best previous performance. Document the performance and key code changes in `iter_X/report.md`. Also analysis why the performance improved or worsened in the report. You should also document the trial and error process of this iteration in the report. 
     - **Git Snapshot (Crucial Step)**: Commit the changes regardless of whether performance improved or worsened. The commit message MUST follow this format: `Iteration X: Time=[Time]s - [Brief Description of Change]` Commands: `git add --all` `git commit -m "Iteration X: Time=12s - Improved cache utilization"`
     f. **Decision Making**:
        - If performance **improved**: Keep the changes. Proceed to next iteration based on this code.
        - If performance **worsened**: 
          - **Revert the code** to the previous good state using git before trying the next idea.
          - Command: `git reset --hard HEAD~1` (Or simply undo changes to files).


# Other docuemnts Provided
- README.md: The README file of the FastGA repository, saved in `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/README.md`. 
- PDF version of the FastGA paper: saved in `/work/shared/users/phd/jl4257/Project/genomics-agent/FASTGA/paper_pdf/Myers and Durbin - FastGA Fast Genome Alignment.pdf`.

