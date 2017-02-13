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
    
    public Frame (int framenr, int k){
        seeds = new TreeMap<>();
        kmers = new TreeMap<>();
        this.framenr = framenr;
        this.k = k;
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
    
    public int getTopSeedposition(TreeMap<Integer,Seed> seeds){
        if(! seeds.isEmpty()){
            int maxSeed = seeds.lastKey();
            for(int i:seeds.keySet()){
                if(seeds.get(i).compareTo(seeds.get(maxSeed))==1){
                    maxSeed = i;
                }
            }
            return maxSeed;
        }else{
            return -1;
        }
    }
    
    public void extendSeeds(){
        TreeMap<Integer, Seed> tempSeedList = new TreeMap(seeds);
        TreeMap<Integer, Kmer> tempKmerList = new TreeMap(kmers);
        while(! tempSeedList.isEmpty()){
            int start = getTopSeedposition(tempSeedList);
            Seed startSeed = tempSeedList.get(start);
            ExtendedSeed extension = new ExtendedSeed(startSeed,getLineage(startSeed.taxon,startSeed.taxonName),startSeed.taxon, framenr);
            tempSeedList.remove(start);
            for(Kmer km: startSeed.getKmer()){
                tempKmerList.remove(km.start);
            }
            int left = start;
            int right = startSeed.end;
            boolean elongatedL = true;
            boolean elongatedR = true;
            while(elongatedL && ! (tempSeedList.headMap(left).isEmpty() && tempKmerList.headMap(left).isEmpty())){
                SortedMap<Integer,Seed> leftSeeds =  tempSeedList.headMap(left);
                SortedMap<Integer,Kmer> leftKmers = tempKmerList.headMap(left);
                if(! leftSeeds.isEmpty()){
                    Seed leftSeed = tempSeedList.get(leftSeeds.lastKey());
                    if(leftSeed.end + 1 == left && (extension.getLineage().contains(leftSeed.taxonName)|| leftSeed.taxonName.equals("root"))){
                        for(Kmer km: leftSeed.getKmer()){
                            tempKmerList.remove(km.start);
                        }
                        extension.addLeft(leftSeed);
                        left = leftSeeds.lastKey();
                        tempSeedList.remove(leftSeeds.lastKey());
                    }else{
                        Kmer leftKmer = tempKmerList.get(leftKmers.lastKey());
                        if(leftKmer.start + 1 == left && (extension.getLineage().contains(leftKmer.taxonName) || leftKmer.taxonName.equals("root"))){
                            extension.addLeft(leftKmer);
                            left = leftKmers.lastKey();
                            tempKmerList.remove(leftKmers.lastKey());
                        }else{
                            elongatedL = false;
                        }
                    }
                }else{
                    Kmer leftKmer = tempKmerList.get(leftKmers.lastKey());
                    if(leftKmer.start + 1 == left && (extension.getLineage().contains(leftKmer.taxonName) || leftKmer.taxonName.equals("root"))){
                        extension.addLeft(leftKmer);
                        left = leftKmers.lastKey();
                        tempKmerList.remove(leftKmers.lastKey());
                    }else{
                        elongatedL = false;
                    }
                }
            }
            while(elongatedR && ! (tempSeedList.tailMap(right).isEmpty() &&  tempKmerList.tailMap(right).isEmpty())){
                TreeMap<Integer,Seed> rightSeeds =  new TreeMap<>(tempSeedList.tailMap(right));
                TreeMap<Integer,Kmer> rightKmers = new TreeMap<>(tempKmerList.tailMap(right));
                if(!rightSeeds.isEmpty()){
                    Seed rightSeed = rightSeeds.firstEntry().getValue();
                    if(rightSeed.start == right + 1 && (extension.getLineage().contains(rightSeed.taxonName)|| rightSeed.taxonName.equals("root"))){
                        for(Kmer km: rightSeed.getKmer()){
                            tempKmerList.remove(km.start);
                        }
                        extension.addRight(rightSeed);
                        right = rightSeed.end;
                        tempSeedList.remove(rightSeeds.firstKey());
                    }else{
                        Kmer rightKmer = tempKmerList.get(rightKmers.firstKey());
                        if(rightKmer.start == right + 1 && (extension.getLineage().contains(rightKmer.taxonName) || rightKmer.taxonName.equals("root"))){
                            extension.addRight(rightKmer);
                            right = rightKmers.firstKey();
                            tempKmerList.remove(rightKmers.firstKey());
                        }else{
                            elongatedR = false;
                        }
                    }
                }else{
                    Kmer rightKmer = tempKmerList.get(rightKmers.firstKey());
                    if(rightKmer.start == right + 1 && (extension.getLineage().contains(rightKmer.taxonName) || rightKmer.taxonName.equals("root"))){
                        extension.addRight(rightKmer);
                        right = rightKmers.firstKey();
                        tempKmerList.remove(rightKmers.firstKey());
                    }else{
                        elongatedR = false;
                    }
                }
            }
            extSeeds.add(extension);
        }
    }
    
    public String getLineage(int taxonID, String TaxonName){
//        try {
//            Runtime rt = Runtime.getRuntime();
//            //Process pr = rt.exec("cmd /c dir");
//            Process pr = rt.exec("/Users/Aranka/edirect/efetch -db taxonomy -id 470 -format xml | /Users/Aranka/edirect/xtract -pattern Taxon -element Lineage");
//            BufferedReader input = new BufferedReader(new InputStreamReader(pr.getInputStream()));
//            String line;
//            while((line=input.readLine()) != null) {
//                System.out.println(line);
//            }
//            int exitVal = pr.waitFor();
//            System.out.println("Exited with error code "+exitVal);
//
//            } catch(IOException | InterruptedException e) {
//                System.out.println(e.toString());
//                e.printStackTrace();
//            }
//        
        String lineage="";
        URL oracle;
        try {
            oracle = new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id="+ taxonID +"&retmode=xml");
            URLConnection yc = oracle.openConnection();
            BufferedReader in = new BufferedReader(new InputStreamReader(yc.getInputStream()));
            String inputLine;
            while ((inputLine = in.readLine()) != null){ 
                if(inputLine.contains("<Lineage>")){
                    lineage = inputLine.trim().substring(9, inputLine.length()-14);
                }
            }
            in.close();
        } catch (MalformedURLException ex) {
            System.out.println(ex.toString());
        } catch (IOException ex) {
            System.out.println(ex.toString());
        }
        return lineage+"; "+TaxonName;
    }
}
