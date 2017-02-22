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
    private final Deque<Kmer> kmers;
    public int start;
    public int end;
    public int frame;
    public int taxon;
    public String rank;
    public String taxonName;
    public ArrayList<String> rankOrder;
    public static final String[] ranks = new String[]{"no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum","superclass", "class", "subclass", "infraclass",
                                                        "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma"};
    public final double[] rankScore = new double[]{0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.4,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.8,0.8,0.9,0.9,1,1,1,1};
    
    public Seed(Deque<Kmer> kmers, int frame){
        this.frame = frame;
        this.kmers = new ArrayDeque<>(kmers);
        this.taxon = kmers.peekFirst().taxonID;
        this.rank = kmers.peekFirst().taxonRank;
        this.start = kmers.peekFirst().start;
        this.taxonName = kmers.peekFirst().taxonName;
        this.end = kmers.peekLast().start;
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
//            if(o.kmers.size() == this.kmers.size()){
//                return (rankOrder.indexOf(this.rank) < rankOrder.indexOf(o.rank) ? -1 :
//                        (rankOrder.indexOf(this.rank) == rankOrder.indexOf(o.rank)? 0 : 1));
//            }else{
//                return (this.kmers.size() < o.kmers.size() ? -1: 1);
//            }
//            if(rankOrder.indexOf(this.rank) == rankOrder.indexOf(o.rank)){
//                return(this.kmers.size() < o.kmers.size() ? -1 :
//                        (this.kmers.size() == o.kmers.size() ?  0 : 1));
//            }else{
//                return (rankOrder.indexOf(this.rank) < rankOrder.indexOf(o.rank) ? -1: 1);
//            }
            double score = calculateScore();
            double oscore = o.calculateScore();
            if(score <= oscore){
                return(score < oscore ? -1 : 0);
            }else{
                return 1;
            }
        }
    }
    
    public double calculateScore(){
        double rankS = this.rankScore[rankOrder.indexOf(this.rank)];
        double length = (double) kmers.size();
        return rankS * length;
    }
    
}
