/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.util.ArrayDeque;
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
    private String lineage;
    public int taxonID;
    public int frameN;
    
    public ExtendedSeed(List<Seed> seeds, String lineage, int taxonID, int frameN){
        kmers = new ArrayDeque<>();
        for (Seed s: seeds){
            kmers.addAll(s.getKmer());
        }
        this.start = kmers.peekFirst().start;
        this.end = kmers.peekLast().start;
        this.lineage = lineage;
        this.taxonID = taxonID;
        this.frameN = frameN;
    }
    
    public ExtendedSeed(Seed s, String lineage, int taxonID, int frameN){
        kmers = new ArrayDeque<>();
        kmers.addAll(s.getKmer());
        this.start = s.start;
        this.end = s.end;
        this.lineage = lineage;
        this.taxonID = taxonID;
        this.frameN = frameN;
    }
    
    public String getLineage(){
        return lineage;
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
}
