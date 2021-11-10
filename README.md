# extractUTRs

This code is used to find 5'UTRs and 3'UTRs using PomBase (**S. pombe**) annotations and a transcriptome assembled with Stringtie2 or RATTLE (https://github.com/comprna/RATTLE).


# Usage

extractUTRs.py [-h] [-t ASSEMBLEDGTF] [-c ANNOTATIONSGTF]
                    [-f RATTLE_FASTQ] [-a ASSEMBLER] [-o OUTPUT]

optional arguments:
  -h, --help            
  
  show this help message and exit
  
  -t , --assembledGTF 
  
                        assembled transcripts in GTF format
                        
  -c , --annotationsGTF 
  
                        CDS annotations in GTF format
                        
  -f , --RATTLE_FASTQ 
  
                        RATTLE output 'transcriptome.fq' in FASTQ format
                        
  -a , --assembler 
  
                        program used to assemble transcripts <Stringtie2> or
                        <RATTLE>
                        
  -o , --output 
  
                        output directory (default: current folder)
                        
