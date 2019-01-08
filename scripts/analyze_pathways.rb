#!/bin/env ruby
require 'byebug'

K = 9

KMER_PWY_FILE = 'kmer_pwys.tsv'
KMER_PROT_FILE = 'kmer_prots.tsv'

Record = Struct.new :header, :sequences

#INPUT = File.open('test.fasta')
INPUT = File.open(ARGV[0])

def stdin_line
  INPUT.readline.strip
end

def each_record(&block)
  while !INPUT.eof?
    header = stdin_line
    raise 'No header' unless header.start_with?('>')
    seqs = [stdin_line]
    5.times do
      raise 'Header mismatch' if header != stdin_line # header should be the same
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
        .map{|kmer, ids| [kmer, ids.to_i]}
        .to_h
  end
end

KMER_HASH = kmer_hash(KMER_PWY_FILE)

# Split into kmers
def split_in_kmers(sequence)
  (sequence.length - K - 1).times.map do |start|
    sequence[start..(start+K-1)]
  end
end

def kmers_to_id(kmers)
  id_count = Hash.new(0)
  kmers.map{ |kmer| KMER_HASH[kmer] }
       .reject(&:nil?)
       .each { |id| id_count[id] += 1 }
  id_count.max_by { |id, count| }
end

def most_likely_id(sequences)
  most_likely = sequences.map{ |seq|   split_in_kmers(seq) }
                         .map{ |kmers| kmers_to_id(kmers)  }
                         .reject(&:nil?)
                         .max_by{ |id, count| count }
  return most_likely.nil? ? 0 : most_likely.first
end

PWY_NAMES = File.open('pwy_name.tsv')
                .readlines
                .map{|line| line.strip.split("\t")}
                .map{|id,name| [id.to_i,name]}
                .to_h

ids_found = Hash.new(0)
each_record do |record|
  id = most_likely_id(record.sequences)
  if id != 0
    ids_found[id] += 1
    puts record.header
    puts id
  end
end

puts "\nReport:\n"
ids_found.to_a
              .sort_by{ |pwy, count| count }
              .each{ |pwy, count| puts "#{PWY_NAMES[pwy]}\t#{pwy}\t#{count}" }



