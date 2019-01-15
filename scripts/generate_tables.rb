#!/usr/bin/env ruby
require 'mysql2'
require 'set'

UNIPROT_BIOCYC_REF_FILE = '/home/rien/Unipept/data/biocyc.tsv'

SQL = Mysql2::Client.new host: 'localhost',
                         username: 'biocyc',
                         password: 'biocyc',
                         database: 'biocyc'

def split_in_kmers(sequence, k)
  (sequence.length - k + 1).times.map{|s| sequence[s..(s+k-1)]}
end

def write_tsv(file, headers, values)
  File.open(file, 'w') do |tsv|
    tsv.write headers.join("\t") + "\n"
    values.each do |row|
      raise '!' if row.length != headers.length
      tsv.write row.join("\t") + "\n"
    end
  end
end

def uniprot_sequence_protid
  refs = File.open(UNIPROT_BIOCYC_REF_FILE)
  id_seq = Hash.new
  refs.each_line do |line|
    _, _, db, id, seq = line.strip.split("\t")
    if db == 'MetaCyc'
      id_seq[id] = seq
    end
  end
  return id_seq
end

def write_pwy_names
  pwy_names = SQL.query(%q{
    SELECT PW.WID,
           DBID.XID,
           PW.Name
    FROM Pathway AS PW
    JOIN DBID
      ON DBID.OtherWID = PW.WID;
    })
  write_tsv 'pwy_name.tsv',
            ["PWY_ID", "PWY_XID", "PWY_NAME"],
            pwy_names.map(&:values)
end

def write_prot_rxn_tsv
  # prot -> enzrxn -> rxn
  prot_rxn = SQL.query(%q{
    SELECT DBID.XID,
           P.WID AS PROT,
           E.WID AS ERXN,
           R.WID AS RXN,
           coalesce(R.ECNumber, R.ECNumberProposed) AS EC
    FROM DBID
    JOIN Protein AS P
      ON DBID.OtherWID = P.WID
    JOIN EnzymaticReaction AS E
      ON P.WID = E.ProteinWID
    JOIN Reaction AS R
      ON E.ReactionWID = R.WID;
  })

  write_tsv 'prot_rxn.tsv',
            ["PROT_NAME", "PROT_ID", "ERXN_ID", "RXN_ID", "EC"],
            prot_rxn.map(&:values)
end

def write_prot_pwy_tsv
  # prot -> enzrxn -> rxn -> pwy_rxn -> pwy
  prot_pwy = SQL.query(%q{
    SELECT PID.XID AS PROT_ID,
           P.WID AS PROT,
           E.WID AS ERXN,
           R.WID AS RXN,
           COALESCE(R.ECNumber, R.ECNumberProposed) AS EC,
           PWY.WID AS PWY,
           PWYID.XID AS PWY_ID
    FROM DBID AS PID
    JOIN Protein AS P
      ON PID.OtherWID = P.WID
    JOIN EnzymaticReaction AS E
      ON P.WID = E.ProteinWID
    JOIN Reaction AS R
      ON E.ReactionWID = R.WID
    JOIN PathwayReaction AS PR
      ON PR.ReactionWID = R.WID
    JOIN Pathway AS PWY
      ON PR.PathwayWID = PWY.WID
    JOIN DBID AS PWYID
      ON PWY.WID = PWYID.OtherWID;
  })
  #write_tsv 'prot_pwy.tsv',
  #          ["PROT_NAME", "PROT_ID", "ERXN_ID", "RXN_ID", "EC", "PWY_ID", "PWY_NAME"],
  #          prot_pwy.map(&:values)
  prot_pwy_hash = Hash.new{ |hash, k| hash[k] = [] }
  prot_pwy.each do |row|
    prot_pwy_hash[row["PROT_ID"]] << row["PWY_ID"]
  end
  write_tsv 'prot_pwy.tsv',
            ["PROT_ID", "PWY_IDS"],
            prot_pwy_hash.transform_values{ |pwy_ids| pwy_ids.join(',') }
end

def protein_sequences
  biocyc = uniprot_sequence_protid

  # prot -> enzrxn -> rxn
  prot_name_id = SQL.query(%q{
    SELECT DBID.XID,
           P.WID
    FROM DBID
    JOIN Protein AS P
      ON DBID.OtherWID = P.WID;
  })

  prot_name_id.map do |row|
      name, id = row.values
      if biocyc.include?(name)
        [name, id, biocyc[name]]
      end
    end
    .reject(&:nil?)
end

def write_rxn_pwys
  # prot -> enzrxn -> rxn -> pwy_rxn -> pwy
  prot_rxns = SQL.query(%q{
    SELECT PID.XID AS PROTID,
           P.WID AS PROT,
           E.WID AS ERXN,
           R.WID AS RXN
    FROM DBID AS PID
    JOIN Protein AS P
      ON PID.OtherWID = P.WID
    JOIN EnzymaticReaction AS E
      ON P.WID = E.ProteinWID
    JOIN Reaction AS R
      ON E.ReactionWID = R.WID;
  })
  prot_rxns.each do |row|
    name, id, erxn, rxn = row.values
    rxn_pwys = SQL.query(%Q{
        SELECT DBID.XID,
               PWY.WID
        FROM PathwayReaction as PR
        WHERE PR.ReactionWID = #{id}
        JOIN Pathway AS PWY
          ON PR.PathwayWID = PWY.WID;
      })
  end
end

def kmer_prots(protseq, k)
  kmers = Hash.new{ |hash, k| hash[k] = [] }
  protseq.each do |prot_name, prot_id, sequence|
    split_in_kmers(sequence, k).uniq.each do |kmer|
      kmers[kmer] << prot_id
    end
  end
  return kmers
end

def protein_reactions
  prot_rxn = Hash.new
  SQL.query(%Q{
      SELECT E.ProteinWID as PROT,
             R.WID AS RXN
      FROM EnzymaticReaction as E
      JOIN Reaction AS R
        ON E.ReactionWID = R.WID;
    })
    .each do |row|
      if prot_rxn.include?(row['PROT'])
        prot_rxn[row['PROT']] << row['RXN']
      else
        prot_rxn[row['PROT']] = [row['RXN']]
      end
    end
  return prot_rxn
end

def reaction_pathways
  rxn_pwy = Hash.new
  SQL.query(%Q{
      SELECT PR.ReactionWID AS RXN,
             PW.WID AS PWY
      FROM PathwayReaction as PR
      JOIN Pathway AS PW
        ON PR.PathwayWID = PW.WID;
    })
    .each do |row|
      if rxn_pwy.include?(row['RXN'])
        rxn_pwy[row['RXN']] << row['PWY']
      else
        rxn_pwy[row['RXN']] = [row['PWY']]
      end
    end
  return rxn_pwy
end

def pathway_ancestry
  SQL.query(%Q{
      SELECT SubPathwayWID AS child,
             SuperPathwayWID AS parent
      FROM SuperPathway;
      }).map(&:values)
        .to_h
end

def kmer_reactions(kmer_prots)
  prot_rxn = protein_reactions
  kmer_prots.map do |kmer, prots|
    reactions = Set.new
    prots.each do |prot_id|
      if prot_rxn.include?(prot_id)
        reactions.merge(prot_rxn[prot_id])
      end
    end
    [kmer, reactions.to_a]
  end
end

def kmer_pathways(kmer_rxns)
  rxn_pwy = reaction_pathways
  parent = pathway_ancestry
  kmer_rxns.map do |kmer, rxns|
    pathways = Set.new
    rxns.each do |rxn_id|
      if rxn_pwy.include?(rxn_id)
        pathways.merge(rxn_pwy[rxn_id])
      end
    end

    # Remove superpathways of which a subpathway
    # is already mapped to this kmer
    to_remove = []
    pathways.each do |pw|
      while parent.include?(pw)
        pw = parent
        to_remove << pw
      end
    end
    pathways -= to_remove

    if pathways.count.positive?
      [kmer, pathways.to_a]
    end
  end.reject(&:nil?)
end

# Returns a hash { n => [k1, k2, ..., km]  } with n the amount of values the 
# kmers k1 .. km are matched with.
def matchcount_kmers(kmers_values)
  counts = Hash.new
  kmers_values.each do |kmer, values|
    count = values.length
    if counts.include?(count)
      counts[count] << kmer
    else
      counts[count] = [kmer]
    end
  end
  return counts
end

# Returns an array of [n, m, p] where m is the amount of kmers (p percent of the total kmers)
# that match with n values.
def kmer_matches_count(kmer_values)
  total = kmer_values.values.sum(&:length)
  kmer_values.to_a
    .map{|matches, values| [matches, values.length, values.length.to_f / total.to_f]}
    .sort
end

def kmer_values_joined_sorted(kmers)
  kmers.map{|kmer, prots| [kmer, prots.join(',')]}.sort
end

puts 'Writing pathway names'
write_pwy_names

puts 'Writing protein pathways'
write_prot_pwy_tsv

puts 'Writing prot_seq.tsv'
protseq = protein_sequences
write_tsv 'prot_seq.tsv',
          ['PROT_NAME', 'PROT_ID', 'PROT_SEQ'],
          protseq

# Protein KMERS
k = 9

puts 'Calculating protein kmers'
protkmers = kmer_prots(protseq, k)

puts 'Writing protein reactions'
write_prot_rxn_tsv

#puts 'Calculation & writing kmer matchcount'

#kmer_prot_matchcount = matchcount_kmers(protkmers)
#write_tsv 'kmer_prot_matchcount.tsv',
#          ["Matchcount", "Count (total: #{protkmers.length})", "Percentage"],
#          kmer_matches_count(kmer_prot_matchcount)

puts 'Sorting & writing kmer_prots.tsv'
write_tsv 'kmer_prots.tsv',
          ['KMER', 'PROT_IDS'],
          kmer_values_joined_sorted(protkmers)

# Reaction KMERS
puts 'Calculating reaction kmers'
rxnkmers = kmer_reactions(protkmers).select{|k,prots| prots.any?}

#puts 'Calculation & writing reaction-kmer matchcount'
#kmer_rxn_matchcount = matchcount_kmers(rxnkmers)
#write_tsv 'kmer_rxn_matchcount.tsv',
#          ["Matchcount", "Count (total: #{rxnkmers.length})", "Percentage"],
#          kmer_matches_count(kmer_rxn_matchcount)

#puts 'Sorting & writing kmer_rxns.tsv'
#write_tsv 'kmer_rxns.tsv',
#          ['KMER', 'RXN_IDS'],
#          kmer_values_joined_sorted(rxnkmers)

# Pathway KMERS
puts 'Calculating pathway kmers'
pwykmers = kmer_pathways(rxnkmers)

#puts 'Calculation & writing pathway-kmer matchcount'

#kmer_pwy_matchcount = matchcount_kmers(pwykmers)
#write_tsv 'kmer_pwy_matchcount.tsv',
#          ["Matchcount", "Count (total: #{pwykmers.length})", "Percentage"],
#          kmer_matches_count(kmer_pwy_matchcount)

puts 'Sorting & writing kmer_pwys.tsv'
write_tsv 'kmer_pwys.tsv',
          ['KMER', 'PWY_IDS'],
          kmer_values_joined_sorted(pwykmers)

# Reaction Pathways
puts 'Writing reaction pathways'
rxn_pwy = reaction_pathways
write_tsv 'rxn_pwy.tsv',
          ['RXN_ID', 'PWY_IDS'],
          rxn_pwy.map{|rxn, pwys| [rxn, pwys.join(',')]}

