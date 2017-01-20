/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package seedextend;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

/**
 *
 * @author Aranka
 */
public class SeedExtend {
    private static final int minSeedSize = 3;
    private static int k;
    private static Map<Integer,List<Seed>> seedsPerFrame;

    /**
     * @param args the command line arguments
     * 0: =0 voor peptiden =k voor k-meren
     * 1: sixframe-file
     * 2: lca-file
     *              
     */
    public static void main(String[] args) {
        BufferedReader sixframe = null;
        try {
            k = Integer.parseInt(args[0]);
            File sixframeF = new File(args[1]);
            sixframe = new BufferedReader(new FileReader(sixframeF));
            Scanner lcaF = new Scanner(new File(args[2]));
            int frameN = 1;
            String header;
            while(frameN <=6){
                header = sixframe.readLine();
                String frame = sixframe.readLine();
                String lcaHeader = lcaF.nextLine();
                Deque<Kmer> deq = new ArrayDeque<>();
                while(lcaHeader.equals(header)&&lcaF.hasNextLine()){
                    String nextLCA = lcaF.nextLine();
                    if(nextLCA.startsWith(">")){
                        lcaHeader = nextLCA;
                    }else{
                        
                    }
                }
            }
        } catch (NumberFormatException | IOException ex) {
            System.out.println("One of both files could not be read");
        }
        
    }
    
}
