# BioSimSpace Examples

This repository contains example `.pdb` files that successful pass through tleap ( [Ambertools](http://ambermd.org/AmberTools.php) ) providing full input and output for a command line use of tleap. These serve as test cases and examples that can be used for setting up systems with BioSimSpace. 

Directory structure:    

    protein   
      |---input    
            |----input files needed for tleap   
      |---output   
            |---- output generated for tleap
        
In the `input` directory contains all necessary information needed for running tleap and the `output` shows the example output generated by running the input commands. Command line arguments should be place in a file called `command.txt` and the input `input.txt`.  
  

## Single molecule system examples 

All the examples below should consist of a single molecule. 

### Protein systems
HSP90 (2JJC)  
MDM2 (1Z1M)      
*Thrombin*   
*Cyclophilin A*   
*Lysozyme*   


### DNA/RNA system
*DNA Hairpin (1D66)*

### Other molecules with parameters (other_molecules)
*water*   
*lipids*   
*ions*   
*co-factors*

### Other molecules without parameters (small-molecules)
benzene   
FXR_79
cyclosporine
HSP90_82p2
HSP90_82p2-sqmfailure
HSP90_82p7
HSP90_84p1-sqmfailure

### Mapping examples (small-molecules at the moment) 
benzene~chlorobenzene*
*benzene~chlorobenzene-spatiallyclosest*
*HSP90_82p2~HSP90_84p1*
*FXR_79~benzene*
*FXR_79~HSP90_82p2* 

