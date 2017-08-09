/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package drawlca4fgs;

import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author Aranka
 */
public class Peptide {
    public String aminoSeq;
    public String taxonName;
    public int taxonID;
    public String taxonRank;
    public int length;
    private double score;
    private Map<String,Double> taxonomy_score = new TreeMap<>();
    public String[] lineage;
    
    public Peptide(String UnipeptOut, Map<String,Double> taxScore){
        String[] info = UnipeptOut.split(",");
        int n = info.length;
        this.taxonRank = info[n-1];
        this.taxonName = info[n-2];
        this.taxonID = Integer.parseInt(info[n-3]);
        this.aminoSeq = info[n-4];
        this.length = aminoSeq.length();
        this.taxonomy_score = taxScore;
        calculateScore();
    }
    public Peptide(String UnipeptOut){
        String[] info = UnipeptOut.split(",");
        int n = info.length;
        this.taxonRank = info[n-1];
        this.taxonName = info[n-2];
        this.taxonID = Integer.parseInt(info[n-3]);
        this.aminoSeq = info[n-4];
        this.length = aminoSeq.length();
    }
    public Peptide(String UnipeptOut, String[] lineage){
        String[] info = UnipeptOut.split(",");
        int n = info.length;
        this.taxonRank = info[n-1];
        this.taxonName = info[n-2];
        this.taxonID = Integer.parseInt(info[n-3]);
        this.aminoSeq = info[n-4];
        this.length = aminoSeq.length();
        this.lineage = lineage;
    }
    public Peptide(String UnipeptOut, Map<String,Double> taxScore, String[] lineage){
        String[] info = UnipeptOut.split(",");
        int n = info.length;
        this.taxonRank = info[n-1];
        this.taxonName = info[n-2];
        this.taxonID = Integer.parseInt(info[n-3]);
        this.aminoSeq = info[n-4];
        this.length = aminoSeq.length();
        this.taxonomy_score = taxScore;
        calculateScore();
        this.lineage = lineage;
    }
    private void calculateScore(){
        double length_score = 0.1;
        if(length>6 && length<=16){
            length_score = 0.09 * (double) length - 0.44;
        }
        if(length >16){
            length_score = 1;
        }
        double taxon_score = taxonomy_score.get(this.taxonRank); 
        score = taxon_score * length_score;
        score = Math.round(score*100.0)/100.0;
    }
    
    public double getScore(){
        return score;
    }

}