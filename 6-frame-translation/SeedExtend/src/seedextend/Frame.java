/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.util.TreeMap;

/**
 *
 * @author Aranka
 */
public class Frame {
    private TreeMap<Integer,Seed> seeds;
    private TreeMap<Integer,Kmer> kmers;
    private int framenr;
    
    public Frame (int framenr){
        seeds = new TreeMap<>();
        kmers = new TreeMap<>();
        this.framenr = framenr;
    }
    
    public int getFramenr(){
        return framenr;
    }
    
    public void fillSeeds(TreeMap<Integer,Seed> seeds){
        this.seeds = seeds;
    }
    
    public void fillKmers(TreeMap<Integer,Kmer> kmers){
        this.kmers = kmers;
    }
    
    public int getTopSeedposition(){
        int maxSeed = seeds.lastKey();
        for(int i:seeds.keySet()){
            if(seeds.get(i).compareTo(seeds.get(maxSeed))==1){
                maxSeed = i;
            }
        }
        return maxSeed;
    }
}
