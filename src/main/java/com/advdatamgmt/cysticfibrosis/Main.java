/*
 * File: Main.java
 * Name: Dalton Voth
 * Date: 9-03-2019
 * Course: CS-490 Bioinformatics
 * Desc: Takes 2 input files in DNA FASTA format and compares for point mutations
 */
package com.advdatamgmt.cysticfibrosis;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author dvoth
 */
public class Main {
    private static Map<String, AminoAcid> aminoAcidMapping;
            
    public static void main(String[] args) {
        String sequence1FileName, sequence2FileName;
        File sequence1File, sequence2File;
        String sequence1, sequence2, aminoSequence1, aminoSequence2;
        String transcribed1, transcribed2;
        String comparisonChoice = "", templateChoice1 = "", templateChoice2 = "";
        String directionChoice1 = "", directionChoice2 = "";
        Scanner scanner = new Scanner(System.in); 
        
        initAminoAcidMapping();
        
        System.out.print("Sequence 1 filename: ");
        sequence1FileName = getFileNameFromUserInput();
        System.out.print("Sequence 2 filename: ");
        sequence2FileName = getFileNameFromUserInput();
        
        while (!templateChoice1.equals("T") && !templateChoice1.equals("N")) {
            System.out.print("Is first sequence a (T)emplate Strand or (N)on-Template Strand? ");
            templateChoice1 = scanner.next();
        }
        
        while (!templateChoice2.equals("T") && !templateChoice2.equals("N")) {
            System.out.print("Is second sequence a (T)emplate Strand or (N)on-Template Strand? ");
            templateChoice2 = scanner.next();
        }
        
        while (!directionChoice1.equals("5") && !directionChoice1.equals("3")) {
            System.out.print("Is first sequence from (5)-3 or (3)-5?");
            directionChoice1 = scanner.next();
        }
        
        while (!directionChoice2.equals("5") && !directionChoice2.equals("3")) {
            System.out.print("Is second sequence from (5)-3 or (3)-5?");
            directionChoice2 = scanner.next();
        }
        
        sequence1File = new File("sequences/" + sequence1FileName);
        sequence2File = new File("sequences/" + sequence2FileName);
        sequence1 = getSequenceFromFile(sequence1File);
        sequence2 = getSequenceFromFile(sequence2File);
        
        if (sequence1.length() != sequence2.length()) {
            System.out.println("Sequences differ in length; exiting...");
            System.exit(0);
        }
        
        sequence1 = orientSequence(sequence1, templateChoice1, directionChoice1);
        sequence2 = orientSequence(sequence2, templateChoice2, directionChoice2);
        
        transcribed1 = convertToRNA(sequence1);
        aminoSequence1 = translateSequence(transcribed1);
        transcribed2 = convertToRNA(sequence2);
        aminoSequence2 = translateSequence(transcribed2);
        
        System.out.println("Sequence 1");
        System.out.println(sequence1);
        System.out.println("Transcribed 1");
        System.out.println(transcribed1);
        System.out.println("Protein 1");
        System.out.println(aminoSequence1);
        
        System.out.println("Sequence 2");
        System.out.println(sequence2);
        System.out.println("Transcribed 2");
        System.out.println(transcribed2);
        System.out.println("Protein 2");
        System.out.println(aminoSequence2);
        
        while (!comparisonChoice.equals("A") && !comparisonChoice.equals("N")) {
            System.out.print("\n\nCompare (N)ucleotide or (A)Amino Acid Sequences? ");
            comparisonChoice = scanner.next();
        }
        
        if (comparisonChoice.equals("A")) {
            findMutations(aminoSequence1, aminoSequence2);
        } else if (comparisonChoice.equals("N")) {
            findMutations(transcribed1, transcribed2);
        }
    }
    
    public static String getFileNameFromUserInput() {
        Scanner scanner = new Scanner(System.in); 
        String inputFileName = "";
        boolean fileExists = false;

        while (!fileExists) {
            inputFileName = scanner.nextLine().trim();
            File input = new File("sequences/" + inputFileName); 
            if (!input.exists()) {
                System.out.print("File " + input.getName() + " not found, try again: ");
            } else {
                fileExists = true;
            }
        }
        
        
        return inputFileName;
    }
    
    public static String getSequenceFromFile(File file) {
        String sequence = "", fileHeader = "";
        
        try {
            Scanner scanner = new Scanner(file);
            // Get the file header
            fileHeader = scanner.nextLine();
            
            // Read the rest of the file as a single string, remove whitespace, to uppercase
            sequence = scanner.useDelimiter("\\Z").next().replaceAll("\\s+","").toUpperCase();
            
            // Output header and sequence separately
            System.out.println("Header: " + fileHeader + "\n");
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            System.exit(0);
        }
        
        return sequence;
    }
    
    public static String orientSequence(String sequence, String template, String direction) {
        if (template.equals("N") && direction.equals("3")) {
            return reverseSequence(sequence);
        } else if (template.equals("T") && direction.equals("5")) {
            sequence = reverseSequence(sequence);
            return complementSequence(sequence);
        } else if (template.equals("T") && direction.equals("3")) {
            return complementSequence(sequence);
        } else {
            return sequence;
        }
    }
    
    public static void initAminoAcidMapping() {
        aminoAcidMapping = new HashMap<>();
        aminoAcidMapping.put("UUU", new AminoAcid("Phenylalanine", "Phe", "F"));
        aminoAcidMapping.put("UUC", new AminoAcid("Phenylalanine", "Phe", "F"));
        aminoAcidMapping.put("UUA", new AminoAcid("Leucine", "Leu", "L"));
        aminoAcidMapping.put("UUG", new AminoAcid("Leucine", "Leu", "L"));
        
        aminoAcidMapping.put("CUU", new AminoAcid("Leucine", "Leu", "L"));
        aminoAcidMapping.put("CUC", new AminoAcid("Leucine", "Leu", "L"));
        aminoAcidMapping.put("CUA", new AminoAcid("Leucine", "Leu", "L"));
        aminoAcidMapping.put("CUG", new AminoAcid("Leucine", "Leu", "L"));
        
        aminoAcidMapping.put("AUU", new AminoAcid("Isoleucine", "Ile", "I"));
        aminoAcidMapping.put("AUC", new AminoAcid("Isoleucine", "Ile", "I"));
        aminoAcidMapping.put("AUA", new AminoAcid("Isoleucine", "Ile", "I"));
        aminoAcidMapping.put("AUG", new AminoAcid("Methionine", "Met", "M"));
        
        aminoAcidMapping.put("GUU", new AminoAcid("Valine", "Val", "V"));
        aminoAcidMapping.put("GUC", new AminoAcid("Valine", "Val", "V"));
        aminoAcidMapping.put("GUA", new AminoAcid("Valine", "Val", "V"));
        aminoAcidMapping.put("GUG", new AminoAcid("Valine", "Val", "V"));
        
        aminoAcidMapping.put("UCU", new AminoAcid("Serine", "Ser", "S"));
        aminoAcidMapping.put("UCC", new AminoAcid("Serine", "Ser", "S"));
        aminoAcidMapping.put("UCA", new AminoAcid("Serine", "Ser", "S"));
        aminoAcidMapping.put("UCG", new AminoAcid("Serine", "Ser", "S"));
        
        aminoAcidMapping.put("CCU", new AminoAcid("Proline", "Pro", "P"));
        aminoAcidMapping.put("CCC", new AminoAcid("Proline", "Pro", "P"));
        aminoAcidMapping.put("CCA", new AminoAcid("Proline", "Pro", "P"));
        aminoAcidMapping.put("CCG", new AminoAcid("Proline", "Pro", "P"));
        
        aminoAcidMapping.put("ACU", new AminoAcid("Threonine", "Thr", "T"));
        aminoAcidMapping.put("ACC", new AminoAcid("Threonine", "Thr", "T"));
        aminoAcidMapping.put("ACA", new AminoAcid("Threonine", "Thr", "T"));
        aminoAcidMapping.put("ACG", new AminoAcid("Threonine", "Thr", "T"));
        
        aminoAcidMapping.put("GCU", new AminoAcid("Alanine", "Ala", "A"));
        aminoAcidMapping.put("GCC", new AminoAcid("Alanine", "Ala", "A"));
        aminoAcidMapping.put("GCA", new AminoAcid("Alanine", "Ala", "A"));
        aminoAcidMapping.put("GCG", new AminoAcid("Alanine", "Ala", "A"));
        
        aminoAcidMapping.put("UAU", new AminoAcid("Tyrosine", "Tyr", "Y"));
        aminoAcidMapping.put("UAC", new AminoAcid("Tyrosine", "Tyr", "Y"));
        aminoAcidMapping.put("UAA", new AminoAcid("Stop", "", ""));
        aminoAcidMapping.put("UAG", new AminoAcid("Stop", "", ""));
        
        aminoAcidMapping.put("CAU", new AminoAcid("Histidine", "His", "H"));
        aminoAcidMapping.put("CAC", new AminoAcid("Histidine", "His", "H"));
        aminoAcidMapping.put("CAA", new AminoAcid("Glutamine", "Gln", "Q"));
        aminoAcidMapping.put("CAG", new AminoAcid("Glutamine", "Gln", "Q"));
        
        aminoAcidMapping.put("AAU", new AminoAcid("Asparagine", "Asn", "N"));
        aminoAcidMapping.put("AAC", new AminoAcid("Asparagine", "Asn", "N"));
        aminoAcidMapping.put("AAA", new AminoAcid("Lysine", "Lys", "K"));
        aminoAcidMapping.put("AAG", new AminoAcid("Lysine", "Lys", "K"));
        
        aminoAcidMapping.put("GAU", new AminoAcid("Aspartic acid", "Asp", "D"));
        aminoAcidMapping.put("GAC", new AminoAcid("Aspartic acid", "Asp", "D"));
        aminoAcidMapping.put("GAA", new AminoAcid("Glutamic acid", "Glu", "E"));
        aminoAcidMapping.put("GAG", new AminoAcid("Glutamic acid", "Glu", "E"));
        
        aminoAcidMapping.put("UGU", new AminoAcid("Cysteine", "Cys", "C"));
        aminoAcidMapping.put("UGC", new AminoAcid("Cysteine", "Cys", "C"));
        aminoAcidMapping.put("UGA", new AminoAcid("Stop", "", ""));
        aminoAcidMapping.put("UGG", new AminoAcid("Tryptophan", "Trp", "W"));
        
        aminoAcidMapping.put("CGU", new AminoAcid("Arginine", "Arg", "R"));
        aminoAcidMapping.put("CGC", new AminoAcid("Arginine", "Arg", "R"));
        aminoAcidMapping.put("CGA", new AminoAcid("Arginine", "Arg", "R"));
        aminoAcidMapping.put("CGG", new AminoAcid("Arginine", "Arg", "R"));
        
        aminoAcidMapping.put("AGU", new AminoAcid("Serine", "Ser", "S"));        
        aminoAcidMapping.put("AGC", new AminoAcid("Serine", "Ser", "S"));
        aminoAcidMapping.put("AGA", new AminoAcid("Arginine", "Arg", "R"));
        aminoAcidMapping.put("AGG", new AminoAcid("Arginine", "Arg", "R"));
        
        aminoAcidMapping.put("GGU", new AminoAcid("Glycine", "Gly", "G"));
        aminoAcidMapping.put("GGC", new AminoAcid("Glycine", "Gly", "G"));
        aminoAcidMapping.put("GGA", new AminoAcid("Glycine", "Gly", "G"));
        aminoAcidMapping.put("GGG", new AminoAcid("Glycine", "Gly", "G"));
    }
    
    public static String translateSequence(String sequence) {
        String aminoSequence = "";
        String aminoAcid = "";
        int startCodonIndex = findStartCodonIndex(sequence);
        
        if (startCodonIndex != -1) {
            for (int i = startCodonIndex; i < sequence.length() - 2; i+=3) {
                char nuc1 = sequence.charAt(i);        
                char nuc2 = sequence.charAt(i + 1);  
                char nuc3 = sequence.charAt(i + 2);  

                aminoAcid = new StringBuilder()
                    .append(nuc1)
                    .append(nuc2)
                    .append(nuc3)
                    .toString();

                if (aminoAcid.equals("UAA") || aminoAcid.equals("UAG") || aminoAcid.equals("UGA")) {
                    break;
                } else {
                    aminoSequence += aminoAcidMapping.get(aminoAcid).getAbbreviation1();
                }
            }
        } else {
            System.out.println("Start codon not found");
        }
        
        return aminoSequence;
    }
    
    public static int findStartCodonIndex(String sequence) {
        String aminoAcid = "";
        
        for (int i = 0; i < sequence.length() - 2; i++){
            char nuc1 = sequence.charAt(i);        
            char nuc2 = sequence.charAt(i + 1);  
            char nuc3 = sequence.charAt(i + 2);  
            
            aminoAcid = new StringBuilder()
                .append(nuc1)
                .append(nuc2)
                .append(nuc3)
                .toString();
            
            if (aminoAcid.equals("AUG")) {
                return i;
            }
        }
        
        return -1;
    }
    
    public static String complementSequence(String sequence) {
        String complementedSequence = "";
        
        for (int i = 0; i < sequence.length(); i++){
            char nucleotide = sequence.charAt(i);
            
            if (nucleotide == 'A') {
                complementedSequence += "U";
            } else if (nucleotide == 'T') {
                complementedSequence += 'A';
            } else if (nucleotide == 'C') {
                complementedSequence += 'G';
            } else if (nucleotide == 'G') {
                complementedSequence += 'C';
            }
        }
        
        return complementedSequence;
    }
    
    public static String convertToRNA(String sequence) {
        String complementedSequence = "";
        
        for (int i = 0; i < sequence.length(); i++){
            char nucleotide = sequence.charAt(i);
            
            if (nucleotide == 'T') {
                complementedSequence += "U";
            } else {
                complementedSequence += nucleotide;
            }
        }
        
        return complementedSequence;
    }
    
    public static String reverseSequence(String sequence) {
        return new StringBuilder(sequence).reverse().toString();
    }
    
    public static void findMutations(String sequence1, String sequence2) {
        boolean mismatchFound = false;
        for (int i=0; i < sequence1.length(); i++) {
            if (sequence1.charAt(i) != sequence2.charAt(i)) {
                mismatchFound = true;
                System.out.println("Mismatch found: " + sequence1.charAt(i) + (i + 1) + sequence2.charAt(i)  );
            }
        }
        
        if (!mismatchFound) {
            System.out.println("Sequences are identical");
        }
    }
}
