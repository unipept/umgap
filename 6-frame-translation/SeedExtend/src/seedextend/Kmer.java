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
        this.taxonName = info[2];
        this.taxonRank = info[3];
        this.start = start;
        this.k = k;
    }
}
