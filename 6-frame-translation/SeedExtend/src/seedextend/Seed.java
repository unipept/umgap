/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.util.Deque;

/**
 *
 * @author Aranka
 */
public class Seed {
    private Deque<Kmer> kmers;
    public int start;
    public int end;
    public int frame;
    
    public Seed(Deque<Kmer> kmers){
        this.kmers = kmers;
        this.start = kmers.peekFirst().start;
        this.end = kmers.peekLast().start + kmers.peekLast().k - 1;
    }
    
    public Deque<Kmer> getKmer(){
        return kmers;
    }
}
