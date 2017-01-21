/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.Objects;

/**
 *
 * @author Aranka
 */
public class Seed implements Comparable<Seed>{
    private Deque<Kmer> kmers;
    public int start;
    public int end;
    public int frame;
    public String rank;
    public ArrayList<String> rankOrder;
    public static final String[] ranks = new String[]{"no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum","superclass", "class", "subclass", "infraclass", "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species_group", "species_subgroup", "species", "subspecies", "varietas", "forma"};
    
    public Seed(Deque<Kmer> kmers){
        this.kmers = new ArrayDeque<>(kmers);
        this.rank = kmers.peekFirst().taxonRank;
        this.start = kmers.peekFirst().start;
        this.end = kmers.peekLast().start + kmers.peekLast().k - 1;
        this.rankOrder = new ArrayList<>(Arrays.asList(ranks));
    }
    
    public Deque<Kmer> getKmer(){
        return kmers;
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 47 * hash + Objects.hashCode(this.kmers);
        hash = 47 * hash + this.start;
        hash = 47 * hash + this.end;
        return hash;
    }
    
    @Override
    public boolean equals(Object o){
        if(!(o instanceof Seed)){
            return false;
        }
        Seed s = (Seed) o;
        return (this.kmers.equals(s.getKmer()) && this.start == s.start && this.end == s.end);
    }

    @Override
    public int compareTo(Seed o) {
        if (equals(o)){
            return 0;
        }else{
            if(o.kmers.size() == this.kmers.size()){
                return (rankOrder.indexOf(this.rank) < rankOrder.indexOf(o.rank) ? -1 :
                        (rankOrder.indexOf(this.rank) == rankOrder.indexOf(o.rank)? 0 : 1));
            }else{
                return (this.kmers.size() < o.kmers.size() ? -1: 1);
            }
        }
    }
    
}
