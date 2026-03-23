#!/usr/bin/env python3
"""Compress ktab partition files into block-compressed format.

Reads existing .ktab.N files (uncompressed), divides entries into blocks,
compresses each block with zstd, writes .ktab.N files in a new format with
a block offset table.

New compressed ktab format:
  int     KMER            (same)
  int64   nels            (negative => indicates compressed format; abs(nels) = actual count)
  int     block_entries   (entries per block)
  int     num_blocks      (number of compressed blocks)
  int64   block_offsets[num_blocks+1]  (byte offset from file start to each compressed block)
  <compressed block 0>
  <compressed block 1>
  ...

The negative nels acts as a format discriminator: old readers see a negative
entry count and will fail with an error rather than silently misreading data.
"""

import struct
import sys
import os
import zstandard as zstd

STREAM_BLOCK = 0x20000  # 131072 entries per block (matches libfastk's read-ahead)


def compress_ktab(input_path, output_path, block_entries=STREAM_BLOCK, level=9):
    """Compress a single ktab partition file."""
    with open(input_path, 'rb') as f:
        # Read header
        kmer = struct.unpack('i', f.read(4))[0]
        nels = struct.unpack('q', f.read(8))[0]
        data = f.read()

    # Compute entry size from data
    if nels <= 0:
        print(f"  Skipping {input_path}: already compressed or empty (nels={nels})")
        return 0, 0

    entry_size = len(data) // nels
    if len(data) % nels != 0:
        print(f"  ERROR: data size {len(data)} not divisible by nels {nels}")
        return 0, 0

    # Divide into blocks
    num_blocks = (nels + block_entries - 1) // block_entries
    compressor = zstd.ZstdCompressor(level=level)
    compressed_blocks = []

    for b in range(num_blocks):
        start = b * block_entries * entry_size
        end = min((b + 1) * block_entries * entry_size, len(data))
        block_data = data[start:end]
        compressed = compressor.compress(block_data)
        compressed_blocks.append(compressed)

    # Compute block offsets (from start of file)
    header_size = (4 +  # KMER
                   8 +  # nels (negative)
                   4 +  # block_entries
                   4 +  # num_blocks
                   (num_blocks + 1) * 8)  # block_offsets

    offsets = [0] * (num_blocks + 1)
    offsets[0] = header_size
    for b in range(num_blocks):
        offsets[b + 1] = offsets[b] + len(compressed_blocks[b])

    # Write output
    with open(output_path, 'wb') as f:
        f.write(struct.pack('i', kmer))
        f.write(struct.pack('q', -nels))  # negative = compressed flag
        f.write(struct.pack('i', block_entries))
        f.write(struct.pack('i', num_blocks))
        for off in offsets:
            f.write(struct.pack('q', off))
        for block in compressed_blocks:
            f.write(block)

    orig_size = 12 + len(data)  # header + data
    comp_size = offsets[-1]
    return orig_size, comp_size


def main():
    if len(sys.argv) < 2:
        print("Usage: compress_ktab.py <genome_root> [--level N] [--block-size N]")
        print("  Compresses .<genome>.ktab.* files in the current directory")
        print("  Creates .<genome>.ktab.*.zst compressed files")
        sys.exit(1)

    root = sys.argv[1]
    level = 9
    block_size = STREAM_BLOCK

    i = 2
    while i < len(sys.argv):
        if sys.argv[i] == '--level':
            level = int(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--block-size':
            block_size = int(sys.argv[i + 1])
            i += 2
        else:
            i += 1

    # Find all ktab parts
    part = 1
    total_orig = 0
    total_comp = 0
    while True:
        input_path = f".{root}.ktab.{part}"
        if not os.path.exists(input_path):
            break
        output_path = f".{root}.ktab.{part}.zst"
        print(f"Compressing {input_path}...")
        orig, comp = compress_ktab(input_path, output_path, block_size, level)
        if orig > 0:
            ratio = orig / comp if comp > 0 else 0
            print(f"  {orig:,} -> {comp:,} bytes ({ratio:.2f}x, {comp * 100 / orig:.1f}%)")
            total_orig += orig
            total_comp += comp
        part += 1

    if total_orig > 0:
        ratio = total_orig / total_comp
        print(f"\nTotal: {total_orig:,} -> {total_comp:,} bytes ({ratio:.2f}x, {total_comp * 100 / total_orig:.1f}%)")
        print(f"Savings: {total_orig - total_comp:,} bytes ({(total_orig - total_comp) / 1048576:.1f} MB)")

    # Optionally replace originals
    if '--replace' in sys.argv:
        print("\nReplacing original files...")
        for p in range(1, part):
            orig = f".{root}.ktab.{p}"
            comp = f".{root}.ktab.{p}.zst"
            if os.path.exists(comp):
                os.rename(comp, orig)
                print(f"  {comp} -> {orig}")


if __name__ == '__main__':
    main()
