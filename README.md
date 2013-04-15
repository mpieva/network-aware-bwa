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
distributes work on the fly to a fully dynamic network of cooperating
processes and machines.

To use the improved workflow, run something like

  bwa bam2bam input.bam -f output.bam -t 8
  
This would run locally with 8 worker threads, handling everything (both single ended
and paired end reads) automatically.  To make use of the network distribution, run

  bwa bam2bam input.bam -f output.bam -t 0 -p PORT
  
where PORT is a free port number.  Then run many instances of this

  bwa worker -h HOST -p PORT -t N
  
where HOST and PORT denote the host running 'bam2bam' and N is the number of desired 
worker threads.  Bwa doesn't case how these jobs are started, so you can use any load 
balancer or grid scheduler you happen to run.
