/*
 * File: AminoAcid.java
 * Name: Dalton Voth
 * Date: 9-03-2019
 * Course: CS-490 Bioinformatics
 * Desc: Structure of an amino acid, used to map 3 nucleotide bases to its
 *       corresponding amino acid codon
 */
package com.advdatamgmt.cysticfibrosis;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author dvoth
 */
public class AminoAcid {
    public String fullName;
    public String abbreviation3;
    public String abbreviation1;

    public AminoAcid(String name, String abr3, String abr1) {
        this.fullName = name;
        this.abbreviation3 = abr3;
        this.abbreviation1 = abr1;
    }
    
    public String getFullName() {
        return fullName;
    }

    public void setFullName(String fullName) {
        this.fullName = fullName;
    }

    public String getAbbreviation3() {
        return abbreviation3;
    }

    public void setAbbreviation3(String abbreviation3) {
        this.abbreviation3 = abbreviation3;
    }

    public String getAbbreviation1() {
        return abbreviation1;
    }

    public void setAbbreviation1(String abbreviation1) {
        this.abbreviation1 = abbreviation1;
    }
}
