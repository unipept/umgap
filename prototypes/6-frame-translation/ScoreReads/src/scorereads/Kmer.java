/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scorereads;

import java.util.Map;

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
    public String[] lineage;
    
    public Kmer(String seq, String taxonName, int taxonID, String taxonRank, int k, Map<String,Double> taxScore){
        this.aminoSeq = seq;
        this.taxonID = taxonID;
        this.taxonName = taxonName;
        this.taxonRank = taxonRank;
        this.k = k;
    }
    
    public Kmer(String Unipeptinfo, int k, Map<String,Double> taxScore){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.k = k;
    }
    
    public Kmer(String Unipeptinfo, int k, Map<String,Double> taxScore, String[] lineage){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.k = k;
        this.lineage=lineage;
    }
    
    
}
