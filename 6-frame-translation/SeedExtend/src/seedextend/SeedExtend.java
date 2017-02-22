/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 *
 * @author Aranka
 */
public class SeedExtend {
    private static final int minSeedSize = 3;
    private static final int gapSize = 2;
    private static int k;
    private static final List<ExtendedSeed> extendedSeeds = new ArrayList<>();
    private static final ArrayList<Frame> frames = new ArrayList<>();

    /**
     * @param args the command line arguments
     * 0: =0 voor peptiden =k voor k-meren
     * 1: sixframe-file
     * 2: lca-file
     *              
     */
    public static void main(String[] args) throws IOException {
        BufferedReader sixframe = null;
        try {
            k = Integer.parseInt(args[0]);
            File sixframeF = new File(args[1]);
            sixframe = new BufferedReader(new FileReader(sixframeF));
            Scanner lcaF = new Scanner(new File(args[2]));
            String lcaHeader = lcaF.nextLine();
            while(lcaF.hasNextLine()){
                String printHeader = lcaHeader.substring(0,lcaHeader.indexOf("|"));
                int frameN = 1;
                while(frameN <= 6){
                    Frame f = new Frame(frameN,k,gapSize);
                    String header;
                    // Initialiseer lijst met seeds, index en deque van kmeren en lees de eerste regels van beide files
                    TreeMap<Integer,Seed> frameSeeds = new TreeMap<>();
                    TreeMap<Integer,Kmer> frameKmers = new TreeMap<>();
                    int prevIndex=0;
                    header = sixframe.readLine();
                    String frame = sixframe.readLine();
                    Deque<Kmer> deq = new ArrayDeque<>();
                    while(lcaHeader.equals(header)&&lcaF.hasNextLine()){
                        // voor alle lca's bij deze frame
                        String nextLCA = lcaF.nextLine();
                        if(nextLCA.startsWith(">")){
                            lcaHeader = nextLCA;
                            if(deq.size() >= 3){
                                Seed newSeed = new Seed(deq, frameN);
                                frameSeeds.put(newSeed.start,newSeed);
                            }
                        }else{
                            // maak de nieuwe kmer
                            Kmer kmer = new Kmer(nextLCA,k);
                            int start = frame.indexOf(kmer.aminoSeq, prevIndex);
                            kmer.setStart(start);
                            frameKmers.put(start,kmer);
                            prevIndex = start;
                            if (!deq.isEmpty()){
                                // voorwaarde om bij de huidige seed te mogen horen
                                if (deq.peekLast().start == start - 1 && deq.peekLast().taxonID == kmer.taxonID){
                                    deq.add(kmer);
                                }else{
                                    // minimale grootte voor een seed
                                    if(deq.size() >= 3){
                                        Seed newSeed = new Seed(deq, frameN);
                                        frameSeeds.put(newSeed.start,newSeed);
                                    }
                                    deq.clear();
                                    deq.add(kmer);
                                }
                            }else{
                                deq.add(kmer);
                            }
                        }
                    }
                    f.fillKmers(frameKmers);
                    f.fillSeeds(frameSeeds);
                    frames.add(f);
                    frameN ++;
                }
                getExtendedSeeds();
                for(ExtendedSeed e:extendedSeeds){
                    System.out.println("Frame: " + e.frameN + "," + e.start + "," + (e.end + 9) + ", " + e.taxonID);
                }
//              TODO:
//                Kies de beste extended seed en geef de aanwezige k-meren terug
            }
        } catch (NumberFormatException | IOException ex) {
            System.out.println("One of both files could not be read");
        }
        
    }
    
    private static void getExtendedSeeds(){
//        TODO:
//          evt. strenger stopcriterium
        while(! frames.get(0).getSeeds().isEmpty()
                || ! frames.get(1).getSeeds().isEmpty()
                || ! frames.get(2).getSeeds().isEmpty()
                || ! frames.get(3).getSeeds().isEmpty()
                || ! frames.get(4).getSeeds().isEmpty()
                || ! frames.get(5).getSeeds().isEmpty()){
            double bestScore = 0;
            int[] bestPos = new int[2];
            for(int i = 0; i < 6; i ++){
                TreeMap<Integer,Seed> seeds = frames.get(i).getSeeds();
                for(Integer pos: seeds.keySet()){
                    Seed test = seeds.get(pos);
                    if(bestScore < test.calculateScore()){
                        bestScore = test.calculateScore();
                        bestPos[0] = i;
                        bestPos[1] = pos;
                    }
                }
            }
            extendedSeeds.add(frames.get(bestPos[0]).extendSeed(bestPos[1]));
        }
    }
    
    private static void printBestSeed(String header){
        ExtendedSeed longest = extendedSeeds.get(0);
        for(int i=1; i< extendedSeeds.size();i++){
            if(extendedSeeds.get(i).getLength()>longest.getLength()){
                longest = extendedSeeds.get(i);
            }
        }
        System.out.println(header + "|Frame " + longest.frameN);
        System.out.println(longest.start + "," + longest.end + ": " + longest.taxonID);
    }
}
