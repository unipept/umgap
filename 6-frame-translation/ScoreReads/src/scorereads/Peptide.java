/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package scorereads;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

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
    
    public Peptide(String aminoSeq, String taxonName, int taxonID, String taxonRank){
        this.aminoSeq = aminoSeq;
        this.taxonName = taxonName;
        this.taxonID = taxonID;
        this.taxonRank = taxonRank;
        this.length = aminoSeq.length();
        try {
            fillTaxScore("taxonomy_score.txt");
        } catch (FileNotFoundException ex) {
            System.out.println("Initiating taxonomy scores failed");
        }
        calculateScore();
    }
    public Peptide(String UnipeptOut, Map<String,Double> taxScore){
        String[] info = UnipeptOut.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.length = aminoSeq.length();
        this.taxonomy_score = taxScore;
        calculateScore();
    }
    public Peptide(String UnipeptOut, Map<String,Double> taxScore, String[] lineage){
        String[] info = UnipeptOut.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.length = aminoSeq.length();
        this.taxonomy_score = taxScore;
        calculateScore();
        this.lineage = lineage;
    }
    private void calculateScore(){
        double length_score = 0.1;
        if(length>6 && length<=16){
            length_score = 0.1 * (double) length - 0.6;
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
    
    private void fillTaxScore(String file) throws FileNotFoundException{
        taxonomy_score.put("no rank", 0.1);
        File input = new File(file);
        Scanner sc = new Scanner(input);
        while(sc.hasNextLine()){
            String[] taxa = sc.next().split(",");
            double s = sc.nextDouble();
            for(String taxon : taxa){
                taxonomy_score.put(taxon, s);
            }
        }
    }
}
