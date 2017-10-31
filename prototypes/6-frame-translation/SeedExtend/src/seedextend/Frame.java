/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 *
 * @author Aranka
 */
public class Frame {
    private TreeMap<Integer,Seed> seeds;
    private TreeMap<Integer,Kmer> kmers;
    private ArrayList<ExtendedSeed> extSeeds = new ArrayList<>();
    private final int framenr;
    private int k;
    private int maxGap;
    private double gapPenalty;
    
    public Frame (int framenr, int k, int gap, double penalty){
        seeds = new TreeMap<>();
        kmers = new TreeMap<>();
        this.framenr = framenr;
        this.k = k;
        this.maxGap = gap;
        this.gapPenalty = penalty;
    }
    
    public int getFramenr(){
        return framenr;
    }
    
    public ArrayList<ExtendedSeed> getExtendedSeeds(){
        return extSeeds;
    }
    
    public void fillSeeds(TreeMap<Integer,Seed> seeds){
        this.seeds = seeds;
    }
    
    public void fillKmers(TreeMap<Integer,Kmer> kmers){
        this.kmers = kmers;
    }
    
    public TreeMap<Integer,Seed> getSeeds (){
        return seeds;
    }
    
    public ExtendedSeed extendSeed(int start){
        Seed startSeed = seeds.get(start);
        ExtendedSeed extension = new ExtendedSeed(startSeed,startSeed.taxon, framenr,gapPenalty);
        seeds.remove(start);
        for(Kmer km: startSeed.getKmer()){
            kmers.remove(km.start);
        }
        int left = start;
        int right = startSeed.end;
        boolean elongatedL = true;
        boolean elongatedR = true;
        while(elongatedL && ! (seeds.headMap(left).isEmpty() && kmers.headMap(left).isEmpty())){
            SortedMap<Integer,Seed> leftSeeds =  seeds.headMap(left);
            SortedMap<Integer,Kmer> leftKmers = kmers.headMap(left);
            if(! leftSeeds.isEmpty()){
                Seed leftSeed = seeds.get(leftSeeds.lastKey());
                for(int i = 0; i<=maxGap; i++){
                    if(leftSeed.end + 1 + i == left){
                        for(Kmer km: leftSeed.getKmer()){
                            kmers.remove(km.start);
                        }
                        extension.addLeft(leftSeed);
                        extension.extraGaps(i);
                        left = leftSeeds.lastKey();
                        seeds.remove(leftSeeds.lastKey());
                        break;
                    }else{
                        Kmer leftKmer = kmers.get(leftKmers.lastKey());
                        if(leftKmer.start + 1 + i == left){
                            extension.addLeft(leftKmer);
                            extension.extraGaps(i);
                            left = leftKmers.lastKey();
                            kmers.remove(leftKmers.lastKey());
                            break;
                        }else{
                            if(i==maxGap){
                                elongatedL = false;
                            }
                        }
                    }
                }
            }else{
                Kmer leftKmer = kmers.get(leftKmers.lastKey());
                for(int i = 0; i<= maxGap;i++){
                    if(leftKmer.start + 1 + i == left ){
                        extension.addLeft(leftKmer);
                        extension.extraGaps(i);
                        left = leftKmers.lastKey();
                        kmers.remove(leftKmers.lastKey());
                        break;
                    }else{
                        if(i==maxGap){
                           elongatedL = false;
                        }
                    }
                }
            }
        }
        while(elongatedR && ! (seeds.tailMap(right).isEmpty() &&  kmers.tailMap(right).isEmpty())){
            TreeMap<Integer,Seed> rightSeeds =  new TreeMap<>(seeds.tailMap(right));
            TreeMap<Integer,Kmer> rightKmers = new TreeMap<>(kmers.tailMap(right));
            if(!rightSeeds.isEmpty()){
                Seed rightSeed = rightSeeds.firstEntry().getValue();
                for(int i = 0; i<= maxGap; i++){
                    if(rightSeed.start == right + 1 + i){
                        for(Kmer km: rightSeed.getKmer()){
                            kmers.remove(km.start);
                        }
                        extension.addRight(rightSeed);
                        extension.extraGaps(i);
                        right = rightSeed.end;
                        seeds.remove(rightSeeds.firstKey());
                        break;
                    }else{
                        Kmer rightKmer = kmers.get(rightKmers.firstKey());
                        if(rightKmer.start == right + 1 + i ){
                            extension.addRight(rightKmer);
                            extension.extraGaps(i);
                            right = rightKmers.firstKey();
                            kmers.remove(rightKmers.firstKey());
                            break;
                        }else{
                            if(i==maxGap){
                                elongatedR = false;
                            }
                        }
                    }
                }
            }else{
                for(int i = 0; i<= maxGap; i++){
                    Kmer rightKmer = kmers.get(rightKmers.firstKey());
                    if(rightKmer.start == right + 1 + i ){
                        extension.addRight(rightKmer);
                        extension.extraGaps(i);
                        right = rightKmers.firstKey();
                        kmers.remove(rightKmers.firstKey());
                        break;
                    }else{
                        elongatedR = false;
                    }
                }
            }
        }
        return extension;
    }
    
}
