/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

/**
 *
 * @author Aranka
 */
public class Kmer {
    public String aminoSeq;
    public String taxonName;
    public int taxonID;
    public String taxonRank;
    public int start;
    public int k;
    
    public Kmer(String Unipeptinfo, int k, int start){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        for(int i=2; i<info.length-1; i++){
            this.taxonName += info[i];
        }
        this.taxonRank = info[info.length-1];
        this.start = start;
        this.k = k;
    }
    public Kmer(String Unipeptinfo, int k){
        String[] info = Unipeptinfo.split(",");
        this.aminoSeq = info[0];
        this.taxonID = Integer.parseInt(info[1]);
        for(int i=2; i<info.length-1; i++){
            this.taxonName += info[i];
        }
        this.taxonRank = info[info.length-1];
        this.k = k;
    }
    
    public void setStart(int start){
        this.start = start;
    }
}
