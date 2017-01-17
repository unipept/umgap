/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scorereads;

import java.util.Map;
import java.util.TreeMap;
import java.net.*;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author Aranka
 */
public class Kmer {
    public String aminoSeq;
    public String taxonName;
    public int taxonID;
    public String taxonRank;
    public int k;
    private double score;
    private Map<String,Double> taxonomy_score = new TreeMap<>();
    public String[] lineage;
    
    public Kmer(String seq, String taxonName, int taxonID, String taxonRank, int k, Map<String,Double> taxScore){
        this.aminoSeq = seq;
        this.taxonID = taxonID;
        this.taxonName = taxonName;
        this.taxonRank = taxonRank;
        this.k = k;
        this.taxonomy_score = taxScore;
        calculateScore();
    }
    
    public Kmer(String Unipeptinfo, int k, Map<String,Double> taxScore){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.k = k;
        this.taxonomy_score = taxScore;
        calculateScore();
    }
    
    public Kmer(String Unipeptinfo, int k, Map<String,Double> taxScore, String[] lineage){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.k = k;
        this.taxonomy_score = taxScore;
        calculateScore();
        this.lineage=lineage;
    }
    
    private void calculateScore(){
        this.score = 0;
        // Nog te implementeren
    }
    
//      Niet meer nodig, vooraf berekenen via unipept cli    
//    private String[] getTaxonInfo(int taxonID){
//        String[] info = new String[2];
//        String url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + taxonID;
//        try {
//            URL taxonomyPage = new URL(url);
//        } catch (MalformedURLException ex) {
//            System.out.println("Url to retreive taxonomy could nog be followed");
//            System.out.println("Given url: " + url);
//        }
//        return info;
//    } 
}
