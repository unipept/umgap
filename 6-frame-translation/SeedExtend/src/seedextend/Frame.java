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
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

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
            ex.printStackTrace();
        } catch (IOException ex) {
            System.out.println(ex.toString());
            ex.printStackTrace();
        }
        return lineage+"; "+TaxonName;
    }
}
