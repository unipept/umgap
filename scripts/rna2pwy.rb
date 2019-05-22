#!/usr/bin/env ruby

# Usage: ./rna2pwy.rb [FILE]
#
# Reads FASTA-records from stdin (or FILE if FILE is given) and tries to assign
# a pathway to each sequence using the 'kmer_pwys.tsv' file in the current
# directory. Per six sequences (resulting from six-frame-translation), one is
# selected with the most assigned pathways as the final assigned pathway for 
# the original read. A single FASTA-record will then be written to stdout with 
# the pathway ID and how many kmers matched that pathway.

K = 9
KMER_PWY_FILE = 'kmer_pwys.tsv'

Record = Struct.new :header, :sequences

INPUT = if ARGV[0]
  File.open(ARGV[0])
else
  STDIN
end

def stdin_line
  INPUT.readline.strip
end

def each_record(&block)
  while !INPUT.eof?
    header = stdin_line
    raise 'No header' unless header.start_with?('>')
    seqs = [stdin_line]
    5.times do
      raise 'Header mismatch' if header != stdin_line
      seqs << stdin_line
    end
    yield Record.new(header, seqs)
  end
end

def kmer_hash(kmer_file)
  File.open(kmer_file) do |file|
    _ = file.readline # discard header
    return file.readlines
        .map(&:strip)
        .map{|line| line.split("\t")}
        .map{|kmer, ids| [kmer, ids.split(',').map(&:to_i)]}
        .to_h
  end
end

KMER_HASH = kmer_hash(KMER_PWY_FILE)

# Split into kmers
def split_in_kmers(sequence)
  (sequence.length - K + 1).times.map{|s| sequence[s..(s+K-1)]}
end

# Map a sequence to ids using kmers and return the kmers that occur
# the most frequent and returns
# [max_count, [list of ids occuring max_count times]]
def seq_to_ids(sequence)
  id_counts = split_in_kmers(sequence)
    .map{ |kmer| KMER_HASH[kmer] }
    .reject(&:nil?)
    .group_by(&:itself)
    .transform_values(&:count)
  max_count = id_counts.values.max
  if max_count
    ids = id_counts.select{ |id, count| count == max_count }.keys
    [ids, max_count]
  end
end

# Takes a list of sequences, tries to map each sequence to
# ids ands then takes the sequence wiith the highest amount of hits
def most_likely_id(sequences)
  ids_counts = sequences.map{ |seq| seq_to_ids(seq) }
                        .reject(&:nil?)
                        .to_h
  max_count = ids_counts.values.max
  if max_count
    ids = ids_counts.select{ |ids, count| count == max_count }
                    .keys
                    .flatten
                    .uniq
    [ids, max_count]
  end
end

each_record do |record|
  ids, count = most_likely_id(record.sequences)
  if ids
    puts record.header
    ids.each do |id|
      puts "#{id} (#{count})"
    end
  end
end

