#!/usr/bin/env ruby

Record = Struct.new :header, :sequences

INPUT = STDIN

def stdin_line
  return nil if INPUT.eof?
  INPUT.readline.strip
end

def each_record(&block)
  line = stdin_line
  while !line.nil?
    header = line
    raise 'No header' unless header.start_with?('>')
    line = stdin_line
    sequences = []
    while !line.nil? && !line.start_with?('>')
      sequences << line
      line = stdin_line
    end
    yield Record.new(header, sequences)
  end
end


PWY_NAMES = File.open('pwy_name.tsv')
                .readlines
                .map{|line| line.strip.split("\t")}
                .map{|id,xid,name| [id.to_i, [xid, name]]}
                .to_h

puts "#MATCHES, PWY_ID, PWY_XID, PWY_NAME"
enum_for(:each_record)
  .map{ |record|
    record.sequences.map{ |seq| seq.split(",").map(&:to_i) }
  }
  .flatten
  .group_by(&:itself)
  .transform_values(&:count)
  .to_a
  .sort_by{ |pwy, count| -count }
  .each{ |pwy, count|
    xid, name = PWY_NAMES[pwy]
    puts [count, pwy, xid, name].join("\t")
  }
