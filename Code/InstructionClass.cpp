//
//  InstructionClass.cpp
//  AssemblyByNumbers
//
//  Created by Avraam Tapinos.
//  Copyright © 2017 Avraam Tapinos. All rights reserved.
//

#include "InstructionClass.hpp"
InstructionClass::InstructionClass(bool c){
    
    std::cout << "AssemblyByNumbers de novo project!!!!\n";
    if(c == true){
        std::cout << "Wrong parameter provided.\n";
    }
    else{
        std::cout << "Help Instructions!!!\n";
    }
    std::cout << "AssemblyByNumbers de novo project usage:\n";
    std::cout << "Parameters:\n";
    std::cout << "-fq [string]                                 -- (Required parameter) fastq file directory and file name\n\n";
    std::cout << "-fo [string]                                 -- (Optional parameter) output directory to save the consensus contigs fasta file. If not provided the fasta file is save in the same as the fastq file directory with the name Contigs.fasta.\n";
    std::cout << "-so [string]                                 -- (Optional parameter) output directory to save the sam file with the sequences overlap alignment sam file. If not provided the sam file is save in the same as the fastq file directory with the name Alignment.sam\n";
    std::cout << "-ovlp [integer]                              -- (Optional parameter) overlap size to use for the analysis [default: 100]\n\n";
    std::cout << "-rep [string]                                -- (Optional parameter) Provide the representation method to use for converting nucleotide sequences to numerical sequences [default: Tetrahedron]\n";
    std::cout << "Options:\n";
    std::cout << "01) Atomic_Numbers\n";
    std::cout << "02) Complex_Numbers\n";
    std::cout << "03) Dna_Walk\n";
    std::cout << "04) EIIP_numbers\n";
    std::cout << "05) Integer_numbers\n";
    std::cout << "06) Pair_numbers\n";
    std::cout << "07) Real_numbers\n";
    std::cout << "08) Tetrahedron             -- default option\n";
    std::cout << "09) Voss_indicators";
    std::cout << "10) Z_curve\n\n";
    std::cout << "-tra [string]                                -- (Optional parameter) Provide the transformation method to use for compressing data to a lower dimensional space [default: DWT]\n";
    std::cout << "Options:\n";
    std::cout << "01) DFT \n";
    std::cout << "02) DWT                     -- default option\n";
    std::cout << "03) PAA\n\n";
    std::cout << "-clvl [integer]                              -- (Optional parameter) compression level to use for the analysis [default: 4]\n\n";
    std::cout << "-knn [integer]                               -- (Optional parameter)the number of knn neighbour data to use for the KNN search [default: 5]\n\n";
    std::cout << "-bin [true/false]                            -- (Optional parameter) Preprocess data and bin them to groups. Reads in each bin is assemble to contigs and then all contigs are re assembled. Proposed for big dataset. [default: false]\n";
    std::cout << "-binmet [int (1 or 2)]                       -- (Optional parameter) if bin is true then binning method will be selected. 1 indicate a faster binning approach and 2 a more accurate clustering approach. [default: 1]\n\n";
    std::cout << "-sen [double]                                -- (Optional parameter) Sensitivity level for matching overlaps based on their to use when aligning reads with smith waterman alignment score. Value must be between -2 and 2. A low value will result in higher sensitivity and higher specificity matching. [default: 1.0]";
    std::cout << "-ssi [true/false]                            -- (Optional parameter) Incrementally make matching process more sensitive [default: false]\n";
    std::cout << "-siv [double]                                -- (Optional parameter) Value to use to increment matching sensitivity level [default: 0.5]\n";
    std::cout << "Example:";
    std::cout << "./AssemblyByumbers -fq ./metagenomic.fastq\n";
    std::cout << "./AssemblyByumbers -fq ./metagenomic.fastq -so ./Results/OVLASEMBLY.sam -fo ./Results/CONTIGS.fasta -ovlp 300 -rep Voss_indicators -tra DFT -clvl 6 -knn 3\n";
    std::cout << "./AssemblyByumbers -fq ./metagenomic.fastq -so ./Results/OVLASEMBLY.sam -fo ./Results/CONTIGS.fasta -ovlp 300 -rep Voss_indicators -tra DFT -clvl 6 -knn 3 -p true\n";
    
};
