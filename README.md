network-aware-bwa
=================

The Burrows Wheeler Aligner with networking capabilities.  The
unmaintainable multi-threading code was removed from upstream BWA and a
new function 'bwa bam2bam' was added.

This function adds a new workflow that reads a BAM file and writes a BAM
file, in a single invocation and with no intermediate files.  It uses
0MQ (ZeroMQ) to allow multithreading and network transparent
multiprocessing.

Why would anyone want this?  First off, the simple workflow makes life
easy.  Second, if you're aligning lots of data and running on a single
machine becomes impractical, you will be happy about this program, which
distributes work on the fly to a fully dinamic network of cooperating
processes and machines.

