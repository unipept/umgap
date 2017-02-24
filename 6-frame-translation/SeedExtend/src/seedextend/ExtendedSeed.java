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
import java.util.List;

/**
 *
 * @author Aranka
 */
public class ExtendedSeed {
    public Deque<Kmer> kmers;
    public int start;
    public int end;
    public int taxonID;
    public int frameN;
    public int ngaps = 0;
    public double score;
    private double gapPenalty = 0.5;
    private ArrayList<String> rankOrder;
    private static final String[] ranks = new String[]{"no rank", "superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum","superclass", "class", "subclass", "infraclass",
                                                        "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "tribe", "subtribe", "genus", "subgenus", "species group", "species subgroup", "species", "subspecies", "varietas", "forma"};
    private final double[] rankScore = new double[]{0.1,0.2,0.2,0.2,0.3,0.3,0.3,0.4,0.4,0.4,0.4,0.5,0.5,0.5,0.5,0.5,0.6,0.6,0.6,0.7,0.7,0.8,0.8,0.9,0.9,1,1,1,1};
    
    
    public ExtendedSeed(List<Seed> seeds, int taxonID, int frameN){
        kmers = new ArrayDeque<>();
        for (Seed s: seeds){
            kmers.addAll(s.getKmer());
        }
        this.start = kmers.peekFirst().start;
        this.end = kmers.peekLast().start;
        this.taxonID = taxonID;
        this.frameN = frameN;
        this.rankOrder = new ArrayList<>(Arrays.asList(ranks));
    }
    
    public ExtendedSeed(Seed s, int taxonID, int frameN){
        kmers = new ArrayDeque<>();
        kmers.addAll(s.getKmer());
        this.start = s.start;
        this.end = s.end;
        this.taxonID = taxonID;
        this.frameN = frameN;
        this.rankOrder = new ArrayList<>(Arrays.asList(ranks));
    }
    
    public void extraGaps(int n){
        ngaps += n;
    }
    
    public void addLeft(Seed s){
        Deque<Kmer> s_kmers = s.getKmer();
        while (! s_kmers.isEmpty()){
            kmers.addFirst(s_kmers.pollLast());
        }
        this.start = kmers.peekFirst().start;
    }
    
    public void addRight(Seed s){
        Deque<Kmer> s_kmers = s.getKmer();
        while (! s_kmers.isEmpty()){
            kmers.addLast(s_kmers.pollFirst());
        }
        this.end = kmers.peekLast().start;
    }
    
    public void addLeft(Kmer k){
        kmers.addFirst(k);
        this.start = kmers.peekFirst().start;
    }
    
    public void addRight(Kmer k){
        kmers.addLast(k);
        this.end = kmers.peekLast().start;
    }
    
    public int getLength(){
        int l = end - start;
        return l;
    }
    
    public void calculateScore(){
        score = 0;
        for(Kmer k:kmers){
            score += rankScore[rankOrder.indexOf(k.taxonRank)];
        }
//        int length = kmers.size() + ngaps;
        score -= ngaps*gapPenalty;
//        score = score/(double) length;
    }
}
